import enum
from abc import ABC, abstractmethod
from collections import namedtuple
from math import sin, cos, asin, copysign
from typing import Callable, List, Tuple

import numpy as np
from scipy.integrate import solve_ivp, cumtrapz
from scipy.misc import derivative

from doublebeam.models import VelocityModel3D, LinearVelocityLayer
from doublebeam.raytracing.ray import Ray3D
from doublebeam.utils import Index


# TODO Use __all__ to specify which symbols can be imported from this module


def cartesian_to_ray_s(x, z, xm, _theta):
    """Conversion of 2D cartesian coordinates x, z to ray coordinate
    s, see eq. 14a in Hill1990 Gaussian beam migration"""
    return (x - xm) * sin(_theta) - z * cos(_theta)


def cartesian_to_ray_n(x, z, xm, theta):
    """Conversion of 2D cartesian coordinates x, z to ray coordinate
    n, see eq. 14b in Hill1990 Gaussian beam migration"""
    return (x - xm) * cos(theta) + z * sin(theta)


def horizontal_slowness(v, theta):
    return sin(theta) / v


def vertical_slowness(v, theta):
    return cos(theta) / v


def critical_angle(v1: float, v2: float) -> float:
    """
    Use Snells law to calculate the critical angle at an interface
    :param v1: velocity before interface in m/s
    :param v2: velocity after interface in m/s
    :return: critical angle in rad
    """
    if v1 < v2:
        return asin(v1/v2)
    return np.pi


def snells_law(p: np.ndarray, v_before: float, v_after: float,
               wave_type: str = "T") -> np.ndarray:
    """
    Calculate initial slowness of a wave transmitted across an interface.
    Eq. (2.4.70) in Cerveny - Seismic ray theory (2001).
    :param p: Slowness vector before interface
    :param v_before: Velocity before the interface
    :param v_after: Velocity after the interface
    :param wave_type: Select for which wave type (reflected/transmitted) the
    slowness should be computed. "T" for transmitted, "R" for reflected.
    """
    # TODO implement checking for critical angle, which should abort ray tracing
    if wave_type == "R":
        # TODO quick exit not valid for non 1D models
        p[Index.Z] = -p[Index.Z]
        return p
    else:
        minus_plus = -1
    # n is normal vector of the interface. n should be oriented to the side
    # the transmitted wave propagates, else the minus/plus relation for
    # transmitted/reflected waves isn't valid.
    n = np.array((0, 0, copysign(1, p[Index.Z])))
    pn = p.dot(n)
    eps = copysign(1, pn)
    return p - (pn + minus_plus * eps
                * (v_after**-2 - v_before**-2 + pn**2)**0.5) * n


class _IVPResultStatus(enum.IntEnum):
    """Enum to make the int status returned from solve_ivp more readable"""
    FAILED = -1
    END_REACHED = 0
    TERMINATION_EVENT = 1


# Tuple containing state of ODE system for kinematic ray tracing
_ODEStateKinematic3D = namedtuple("ODEStateKinematic3D", ["x", "y", "z", "px", "py", "pz", "T"])

# Tuple containing state of ODE system for dynamic ray tracing
_ODEStateDynamic3D = namedtuple("ODEStateDynamic3D", "x, y, z, px, py, pz, T")

_IVPEventFunction = Callable[[float, _ODEStateKinematic3D], float]


class RayTracerBase(ABC):

    def __init__(self, velocity_model: VelocityModel3D):
        self.model = velocity_model
        self.layer = None

    # TODO change signature of _velocity to accept np array instead of
    #  three floats
    @staticmethod
    def _velocity(layer: LinearVelocityLayer, x: float, y: float,
                  z: float) -> float:
        """
        Evaluate velocity at a given depth in a layer with linearly varying
        velocity
        :param layer: The layer to evaluate
        """
        return layer["intercept"] + layer["gradient"] * z

    @abstractmethod
    def _trace_layer(self, ray: Ray3D, slowness: np.ndarray,
                     max_step_s: float) -> None:
        pass

    def trace_stack(self, ray: Ray3D, ray_code: str = None,
                    max_step: float = 1) -> None:
        """
        Trace ray through a stack of layers. The ray type at an interface is
        chosen by the ray code.
        :param ray: Ray to trace through the model
        :param ray_code: Specifies which ray (Transmitted/Reflected) to follow
        at an interface in the model. "T" stands for the transmitted ray, "R"
        for the reflected ray. If not given or empty, the ray will only be
        traced through the layer in which its starting point resides.
        :param max_step: Max step s for the integration.
        """
        top, bottom = self.model.vertical_boundaries()
        if not top <= ray.start[Index.Z] <= bottom:
            raise ValueError(f"Ray {ray} starts outside of model")

        index = self.model.layer_index(ray.start[Index.Z])
        self.layer = self.model[index]
        initial_slownesses = ray.last_slowness
        self._trace_layer(ray, initial_slownesses, max_step)
        if not ray_code:
            return
        for wave_type in ray_code:
            v_top, v_bottom = self.model.interface_velocities(ray.last_point[Index.Z])
            if wave_type == "T":
                # reflected waves stay in the layer and dont change the index
                index += -1 if ray.direction == "up" else 1
            self.layer = self.model[index]
            if ray.direction == "down":
                new_p = snells_law(ray.last_slowness, v_top, v_bottom, wave_type)
            else:
                new_p = snells_law(ray.last_slowness, v_bottom, v_top, wave_type)
            self._trace_layer(ray, new_p, max_step)


class KinematicRayTracer3D(RayTracerBase):

    def _trace(self, s: float, y: _ODEStateKinematic3D) -> _ODEStateKinematic3D:
        """
        Implement ray tracing system (3.1.10) from Cerveny - Seismic ray theory
        (2001).
        :param s: Arclength along ray
        :param y: Previous state
        :return: New state
        """
        x, y, z, px, py, pz, T = y
        v = self._velocity(self.layer, x, y, z)
        dxds = px * v
        dyds = py * v
        dzds = pz * v
        # TODO creating a lambda every iteration is probably not conducive to
        #  performance.
        dpxds = derivative((lambda x_: 1 / self._velocity(self.layer, x_, y, z)),
                           x, dx=0.0001)
        dpyds = derivative((lambda y_: 1 / self._velocity(self.layer, x, y_, z)),
                           y, dx=0.0001)
        dpzds = derivative((lambda z_: 1 / self._velocity(self.layer, x, y, z_)),
                           z, dx=0.0001)
        dTds = 1. / v
        return _ODEStateKinematic3D(dxds, dyds, dzds, dpxds, dpyds, dpzds, dTds)

    def _trace_layer(self, ray: Ray3D, initial_slowness: np.ndarray,
                     max_step_s: float) -> None:
        """
        Trace a ray through a single layer of the model and set the parameters
        (path, slowness, traveltime) of the ray.
        :param ray: Ray to trace
        :param initial_slowness: Initial value of px, py, pz (slowness along
        the corresponding axis) of the ray in s/m.
        :param max_step_s: Max step the solver takes for the integration
        variable s.
        """
        if self.layer["gradient"] == 0:
            # constant velocity layer, use analytic ray tracing
            # this uses (A.1) and (4) from Analytical ray tracing system:
            # Introducing art, a C-library designed for seismic applications
            # (Miqueles et al., 2013)
            c = self.layer["intercept"]
            p0 = initial_slowness
            x0 = ray.last_point
            z0 = ray.last_point[Index.Z]
            z = self.layer["top_depth"] if p0[Index.Z] < 0 else self.layer["bot_depth"]
            s_end = (z - z0) / (c * p0[Index.Z])
            num_steps = int(s_end / max_step_s)
            s = np.linspace(0, s_end, num_steps)
            path = x0[None, :] + c * s[:, None] * p0[None, :]
            time = ray.last_time + s / c
            ray.path.append(path)
            ray.slowness.append(p0.reshape(1, 3))
            ray.travel_time.append(time)
            return
        # workaround: make events active only after a step is taken because
        # sometimes events trigger directly after tracing starts and the
        # algorithm thinks an interface was reached. Keep event function
        # artificially away from zero to avoid that
        upper_layer_event: _IVPEventFunction = lambda s, y_: y_[2] - self.layer["top_depth"] if s > max_step_s else 1
        lower_layer_event: _IVPEventFunction = lambda s, y_: y_[2] - self.layer["bot_depth"] if s > max_step_s else -1
        for func in (upper_layer_event, lower_layer_event):
            func.terminal = True
        initial_state = _ODEStateKinematic3D(*ray.last_point, *initial_slowness,
                                             ray.last_time)
        result = solve_ivp(self._trace, (0, np.inf), initial_state,
                           max_step=max_step_s, events=(upper_layer_event,
                                                        lower_layer_event))
        x, y, z, px, py, pz, t = result.y
        ray.path.append(np.vstack((x, y, z)).T)
        ray.slowness.append(np.vstack((px, py, pz)).T)
        ray.travel_time.append(t)


class DynamicRayTracer3D(RayTracerBase):

    def _trace(self, s: float, y: _ODEStateDynamic3D) -> _ODEStateDynamic3D:
        # TODO same as kinematic ray tracing, maybe move to base class
        x, y, z, px, py, pz, T = y
        v = self._velocity(self.layer, x, y, z)
        dxds = px * v
        dyds = py * v
        dzds = pz * v
        dpxds = derivative((lambda x_: 1 / self._velocity(self.layer, x_, y, z)),
                           x, dx=0.0001)
        dpyds = derivative((lambda y_: 1 / self._velocity(self.layer, x, y_, z)),
                           y, dx=0.0001)
        dpzds = derivative((lambda z_: 1 / self._velocity(self.layer, x, y, z_)),
                           z, dx=0.0001)
        dTds = 1. / v
        return _ODEStateDynamic3D(dxds, dyds, dzds, dpxds, dpyds, dpzds, dTds)

    def _trace_layer(self, ray: Ray3D, initial_slowness: np.ndarray,
                     max_step_s: float) -> None:
        # TODO factor out the common code between dynamic and kinematic ray tracing
        #  for creating events by putting them into a function in the base class
        upper_layer_event: _IVPEventFunction = lambda s, y_: y_[2] - self.layer["top_depth"] if s > max_step_s else 1.
        lower_layer_event: _IVPEventFunction = lambda s, y_: y_[2] - self.layer["bot_depth"] if s > max_step_s else -1.
        for func in (upper_layer_event, lower_layer_event):
            func.terminal = True
        V0 = self.model.eval_at(*ray.last_point)
        # for a layer with constant gradient of velocity, P is constant
        P = np.array([1j/V0, 0, 0, 1j/V0]).reshape(2, 2)
        beam_width_m = 10
        beam_frequency_Hz = 40
        Q = np.array([beam_frequency_Hz*beam_width_m**2 / V0, 0,
                      0, beam_frequency_Hz*beam_width_m**2 / V0],
                     dtype=np.complex128).reshape(2, 2)
        initial_state = _ODEStateDynamic3D(*ray.last_point, *initial_slowness,
                                           ray.last_time)
        result = solve_ivp(self._trace, (0, np.inf), initial_state,
                           max_step=max_step_s, events=(upper_layer_event,
                                                        lower_layer_event))
        x, y, z, px, py, pz, t = result.y
        # make P multi dimensional according to the number of steps by adding an
        # empty last axis and repeating the array along it
        num_steps = len(t)
        P = np.repeat(P[..., None], num_steps, axis=-1)
        # this make the last axis the number of points
        # move number of points as first axis
        P = np.moveaxis(P, -1, 0)
        # use special case for layer with constant gradient of velocity
        # see Cerveny2001, section 4.8.3
        V = np.array([self.model.eval_at(*point) for point in zip(x, y, z)])
        sigma = cumtrapz(V**2, t, initial=0)
        Q = Q[np.newaxis, ...]
        Q = Q + sigma[..., np.newaxis, np.newaxis] * P
        self.P.append(P)
        self.Q.append(Q)
        ray.path.append(np.vstack((x, y, z)).T)
        ray.slowness.append(np.vstack((px, py, pz)).T)
        ray.travel_time.append(t)

    def trace_stack(self, ray: Ray3D, ray_code: str = None,
                    max_step: float = 1) -> Tuple[np.ndarray, np.ndarray]:
        self.P: List[np.ndarray] = []
        self.Q: List[np.ndarray] = []
        super().trace_stack(ray, ray_code, max_step)
        # TODO mutable list is not cleaned after exit
        return self.P, self.Q


