import enum
from abc import ABC, abstractmethod
from collections import namedtuple
from math import sin, cos, asin, copysign
from typing import Callable, List, Tuple, Optional

import numpy as np
from scipy.integrate import solve_ivp, cumtrapz
from scipy.misc import derivative

from doublebeam.models import VelocityModel3D, LinearVelocityLayer
from doublebeam.raytracing.ray import Ray3D
from doublebeam.utils import Index, DIGITS_PRECISION, angle


# TODO Use __all__ to specify which symbols can be imported from this module.
#  Even better would be to create the ability for a top level import from doublebeam package


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
        self.index: int = None

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

        self.index = self.model.layer_index(ray.start[Index.Z])
        self.layer = self.model[self.index]
        initial_slownesses = ray.last_slowness
        self._trace_layer(ray, initial_slownesses, max_step)
        if not ray_code:
            return
        for wave_type in ray_code:
            new_p = self.continue_ray_across_interface(ray, wave_type)
            self._trace_layer(ray, new_p, max_step)
    
    def continue_ray_across_interface(self, ray: Ray3D, wave_type: str):
        v_top, v_bottom = self.model.interface_velocities(ray.last_point[Index.Z])
        if wave_type == "T":
            # reflected waves stay in the layer and dont change the index
            self.index += -1 if ray.direction == "up" else 1
        self.layer = self.model[self.index]
        if ray.direction == "down":
            new_p = snells_law(ray.last_slowness, v_top, v_bottom, wave_type)
        else:
            new_p = snells_law(ray.last_slowness, v_bottom, v_top, wave_type)
        return new_p


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
            z0 = x0[Index.Z]
            z = self.layer["top_depth"] if p0[Index.Z] < 0 else self.layer["bot_depth"]
            s_end = (z - z0) / (c * p0[Index.Z])
            num_steps = int(s_end / max_step_s)
            s = np.linspace(0, s_end, num_steps)
            path = x0[None, :] + c * s[:, None] * p0[None, :]
            time = ray.last_time + s / c
            ray.path.append(path)
            ray.slowness.append(p0.reshape(1, 3))
            ray.travel_time.append(time)
            # workaround numerical issue where target is sometimes "overshot"
            # TODO deduplicate this for constant velocity and linear gradient
            if (path.T[Index.Z][-1] > self.layer["bot_depth"]
                    or path.T[Index.Z][-1] < self.layer["top_depth"]):
                path.T[Index.Z][-1] = round(path.T[Index.Z][-1], ndigits=DIGITS_PRECISION)
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
        # workaround numerical issue where target is sometimes "overshot"
        # TODO deduplicate this for constant velocity and linear gradient
        if z[-1] > self.layer["bot_depth"] or z[-1] < self.layer["top_depth"]:
            z[-1] = round(z[-1], ndigits=DIGITS_PRECISION)
        ray.path.append(np.vstack((x, y, z)).T)
        ray.slowness.append(np.vstack((px, py, pz)).T)
        ray.travel_time.append(t)


class InterfacePropagatorMatrix:
    """
    This class implements the required methods to transform matrices P and Q,
    which are calculated during dynamic ray tracing, across an interface.
    The main equation used is eq. 4.4.75 from Cerveny2001.
    Values corresponding to a point behind the interface are marked with a tilde.

    Table of symbols
    ----------------------
    x_i: general Cartesian coordinate system.
    z_i: coordinates of the local Cartesian coordinate system at Q.
    i_n^(z): basis vectors of the local Cartesian coordinate system, with its
        origin at point Q, defined by eq. 4.4.21 Cerveny2001.
        Its unit vector i_3^(z) coincides with the unit vector normal to the
        interface at Q. The other two vectors are specified arbitrarily in the
        tangent plane so that the system is mutually orthogonal and right-handed.
        For a velocity model containing only horizontal interfaces, the i_3^(z)
        axis coincides with the z-axis of the general Cartesian coordinate
        system.
    q_i: Coordinates of the ray centered coordinate system.
    e_i: Basis vectors of the ray centered coordinate system of the incident
        wave at Q.
        e_3 is the unit tangent to the ray. e_1 and e_2 are perpendicular to the
        ray. General definition can be found in chapter 4.1.1, Cerveny2001.
    kappa: angle between e_2 and i_2^(z) (p. 293 Cerveny2001)
        0 <= kappa <= 2*pi
        cos(kappa) = e_2 * i_2^(z)
        sin(kappa) = e_1 * i_2^(z)
    For the special case of V = V(z) and horizontal interfaces:
        Rays will be planar, e_2 will be constant and can be chosen to be
        orthogonal to the ray plane. Using 4.4.21 for i_2 means that it also
        will be orthogonal to the ray plane. This means kappa is 0.
    # TODO find a good compromise between different style of member functions and attributes
    """

    def __init__(self):
        # cache these in class to save array creation and only change one value
        self.G_parallel = np.identity(2)
        self.G_parallel_tilde = np.identity(2)
        self.G_orthogonal: np.ndarray
        self.G_orthogonal_tilde: np.ndarray


    def u(self, wave_type: str, V: float, V_tilde: float, i_S: float,
          i_R: float, epsilon: float) -> float:
        """
        Eq. 4.4.51 from Cerveny2001
        :param wave_type: String specifying wave type, valid values are "T" for
        transmitted and "R" for reflected
        :param V: Velocity before the interface, in m/s
        :param V_tilde: Velocity after the interface, in m/s
        :param i_S: Acute angle of incidence, 0 <= i_S <= pi/2
        :param i_R: Acute angle of reflection/transmission
        """

        minus_plus = -1 if wave_type == "T" else 1
        return epsilon * (1. / V * cos(i_S) + minus_plus * 1. / V_tilde * cos(i_R))

    def G(self, epsilon: float, i_S: float, i_R: float, wave_type: str) -> np.ndarray:
        """
        Left eq. from 4.4.48, Cerveny2001
        :return:
        """
        self.G_parallel[0][0] = epsilon * cos(i_S)
        # upper sign probably corresponds to transmitted wave
        plus_minus = 1 if wave_type == "T" else -1
        self.G_orthogonal[0][0] = plus_minus * epsilon * cos(i_R)
        return self.G_parallel @ self.G_orthogonal

    def G_tilde(self, kappa: float) -> np.ndarray:
        """
        Right equation from 4.4.48, Cerveny2001
        :return:
        """
        cos_kappa, sin_kappa = cos(kappa), sin(kappa)
        self.G_orthogonal = np.array(((cos_kappa, -sin_kappa),
                                      (sin_kappa, cos_kappa)))
        self.G_orthogonal_tilde = self.G_orthogonal
        return self.G_parallel_tilde @ self.G_orthogonal_tilde

    def E(self, V, i_S, i_R, epsilon, old_gradient) -> np.ndarray:
        """
        Eq. 4.4.53 from Cerveny2001
        :return:
        """
        # TODO modify this to work with a more general velocity model
        # dV_dzi means the derivative of the velocity after the z_i coordinate
        # for V=V(z)
        dV_dz1 = 0
        dV_dz2 = 0
        dV_dz3 = old_gradient
        E11 = -sin(i_S) * V**-2 * ((1 + cos(i_S)**2) * dV_dz1 - epsilon * cos(i_S) * sin(i_S) * dV_dz3)
        E12 = -sin(i_S) * V**-2 * dV_dz2
        E22 = 0
        return np.array(((E11, E12),
                         (E12, E22)))

    def E_tilde(self, wave_type: str, V_tilde, i_R, i_S, epsilon, new_gradient):
        """
        Eq. 4.4.54 from Cerveny2001
        :return:
        """
        dV_tilde_dz1 = 0
        dV_tilde_dz2 = 0
        dV_tilde_dz3 = new_gradient
        minus_plus = -1 if wave_type == "R" else 1
        E11 = -sin(i_R) * V_tilde**-2 * ((1 + cos(i_R)**2) * dV_tilde_dz1 + minus_plus * epsilon * cos(i_R) * sin(i_R) * dV_tilde_dz3)
        E12 = -sin(i_R) * V_tilde**-2 * dV_tilde_dz2
        E22 = 0
        return np.array(((E11, E12),
                         (E12, E22)))

    def D(self) -> np.ndarray:
        """
        Eq. 4.4.15 from Cerveny2001
        For the currently implemented velocity layer with horizontal interfaces
        only this function is zero everywhere since it contains the second
        derivative of the interface function Sigma in the numerator.
        Sigma = Sigma(z3) for horizontal interfaces.
        :return:
        """
        return np.zeros((2, 2))

    def do(self, P: np.ndarray, Q: np.ndarray, i_S: float, i_R: float,
                 wave_type: str, V_before: float, V_after: float,
                 epsilon:float, old_gradient: float, new_gradient: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Transform P, Q with interface propagator matrix.
        :param P:
        :param Q:
        :return: Tuple of transformed matrices P, Q
        """
        # TODO fix order of operations where G_tilde has to be called before G
        # TODO this is only valid for the simple velocity model V = V(z) and horizontal interfaces
        kappa = 0.
        G_tilde = self.G_tilde(kappa)
        G = self.G(epsilon, i_S, i_R, wave_type)
        G_inverted = np.linalg.inv(G)
        G_tilde_inverted = np.linalg.inv(G_tilde)
        E = self.E(V_before, i_S, i_R, epsilon, old_gradient)
        E_tilde = self.E_tilde(wave_type, V_after, i_R, i_S, epsilon, new_gradient)
        u = self.u(wave_type, V_before, V_after, i_S, i_R, epsilon)
        D = self.D()
        # eq. 4.4.67
        P_tilde = G_tilde_inverted @ (G @ P + (E - E_tilde - u * D) @ G_inverted.T @ Q)
        # eq. 4.4.64
        Q_tilde = G_tilde.T @ G_inverted.T @ Q
        return P_tilde, Q_tilde


class DynamicRayTracer3D(KinematicRayTracer3D):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.P: List[np.ndarray] = []
        self.Q: List[np.ndarray] = []
        self.ic = InterfacePropagatorMatrix()
        
    def _trace_layer(self, ray: Ray3D, initial_slowness: np.ndarray,
                     max_step_s: float) -> None:
        super()._trace_layer(ray, initial_slowness, max_step_s)
        # make P multi dimensional according to the number of steps by adding an
        # empty last axis and repeating the array along it
        num_steps = len(ray.travel_time[-1])
        P0 = np.repeat(self.P0[..., None], num_steps, axis=-1)
        # this make the last axis the number of points
        # move number of points as first axis
        P0 = np.moveaxis(P0, -1, 0)
        # use special case for layer with constant gradient of velocity
        # see Cerveny2001, section 4.8.3
        V = np.array([self.model.eval_at(*point) for point in ray.path[-1]])
        sigma = cumtrapz(V**2, ray.travel_time[-1], initial=0)
        Q0 = self.Q0[np.newaxis, ...]
        Q0 = Q0 + sigma[..., np.newaxis, np.newaxis] * P0
        self.P.append(P0)
        self.Q.append(Q0)

    def continue_ray_across_interface(self, ray: Ray3D, wave_type: str):
        # TODO modify unit vector of interface for a more general velocity model
        interface_unit_vector = np.array((0, 0, 1.))
        V_top, V_bottom = self.model.interface_velocities(ray.last_point[Index.Z])
        i_S = angle(ray.last_slowness, interface_unit_vector)
        direction = ray.direction
        old_gradient = self.model.gradients[self.index]
        new_slowness = super().continue_ray_across_interface(ray, wave_type)
        new_gradient = self.model.gradients[self.index]
        i_R = angle(new_slowness, interface_unit_vector) if wave_type == "T" else i_S
        # epsilon is introduced by eq. 2.4.71, Cerveny2001
        epsilon = copysign(1., ray.last_slowness @ interface_unit_vector)
        if direction == "down":
            V_before, V_after = V_top, V_bottom
        else:
            V_before, V_after = V_bottom, V_top
        if wave_type == "R":
            V_after = V_before

        self.P0, self.Q0 = self.ic.do(self.last_P, self.last_Q, i_S, i_R, wave_type,
                                    V_before, V_after, epsilon, old_gradient,
                                    new_gradient)

        return new_slowness

    # TODO change function signature to take GaussBeam class instead which is an extended Ray
    def trace_stack(self, ray: Ray3D, ray_code: str = None,
                    max_step: float = 1) -> Tuple[List[np.ndarray],
                                                  List[np.ndarray]]:
        V0 = self.model.eval_at(*ray.last_point)
        beam_width_m = 10
        beam_frequency_Hz = 40
        # for a layer with constant gradient of velocity, P is constant
        self.P0 = np.array([1j/V0, 0, 0, 1j/V0]).reshape(2, 2)
        self.Q0 = np.array([beam_frequency_Hz*beam_width_m**2 / V0, 0,
                      0, beam_frequency_Hz*beam_width_m**2 / V0],
                      dtype=np.complex128).reshape(2, 2)
        super().trace_stack(ray, ray_code, max_step)
        P, Q = self.P, self.Q
        self.P, self.Q = [], []
        return P, Q

    @property
    def last_P(self) -> Optional[np.ndarray]:
        try:
            return self.P[-1][-1]
        except IndexError:
            return None

    @property
    def last_Q(self) -> Optional[np.ndarray]:
        try:
            return self.Q[-1][-1]
        except IndexError:
            return None