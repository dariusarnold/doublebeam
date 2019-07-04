import enum
from collections import namedtuple
from math import sin, cos, asin, acos, radians, copysign
from typing import Callable, List, Union

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.misc import derivative

from doublebeam.core.models import VelocityModel3D, LinearVelocityLayer


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


class Ray3D:

    def __init__(self, start_x: float, start_y: float, start_z: float,
                 theta: float, phi: float):
        """
        :param start_x: x coordinate of start point of ray in m
        :param start_y: y coordinate of start point of ray in m
        :param start_z: z coordinate of start point of ray in m
        :param theta: Angle against downgoing vertical axis (z) at start
        point in rad, increasing upwards. 0 <= theta <= pi
        :param phi: Angle against x axis at start point in rad, with increasing
        angle towards the y axis
        0 <= phi <= 2*pi
        """
        self.start = np.array((start_x, start_y, start_z))
        self.theta = theta
        self.phi = phi
        # These will be set after the ray is traced
        self.path: List[np.ndarray] = []
        self.slowness: List[np.ndarray] = []
        self.travel_time: List[np.ndarray] = []

    @property
    def last_point(self) -> np.ndarray:
        """
        Return last point of the rays path. If the ray hasn't been traced yet,
        return the starting point of the ray.
        """
        try:
            return self.path[-1][-1]
        except IndexError:
            return self.start

    @property
    def last_slowness(self) -> Union[np.ndarray, None]:
        """
        Return last slowness value. If the ray hasn't been traced yet, return
        None.
        """
        try:
            return self.slowness[-1][-1]
        except IndexError:
            return None

    @property
    def last_time(self) -> float:
        """
        Return last travel time of the ray.
        :return:
        """
        try:
            return self.travel_time[-1][-1]
        except IndexError:
            return 0.

    @property
    def direction(self) -> str:
        """
        Return direction of the ray along the z-axis. A downgoing ray will
        return "down", while an upgoing ray will return "up". The special case
        of a horizontal ray (vertical slowness pz = 0) will return "horizontal".
        """
        try:
            pz = self.last_slowness[Index.Z]
        except TypeError:
            # if None is returned (ray hasn't been traced yet), use starting
            # angle against vertical axis
            pz = np.pi/2 - self.theta
        if pz < 0:
            return "up"
        elif pz > 0:
            return "down"
        return "horizontal"


def initial_slowness3D(ray: Ray3D, v0: float) -> np.ndarray:
    """
    Calculate initial vertical and horizontal slowness for a ray.
    For geometric definitions see chapter 3.2.1 in Cerveny - Seismic ray theory
    :param ray: Ray instance
    :param v0: Velocity at ray starting point
    :return: Tuple of slowness values
    """
    px = 1/v0 * sin(ray.theta) * cos(ray.phi)
    py = 1/v0 * sin(ray.theta) * sin(ray.phi)
    pz = 1/v0 * cos(ray.theta)
    return np.array((px, py, pz))


def length(vector: np.ndarray) -> float:
    """Calculate length of vector"""
    return (vector @ vector)**0.5


def angle(vector1: np.ndarray, vector2: np.ndarray) -> float:
    """Calculate angle between two vectors in rad"""
    return acos((vector1 @ vector2) / (length(vector1) * length(vector2)))


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
        return p * np.array((1, 1, -1))
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


def plot_ray_in_model(ray: Ray3D, velocity_model: VelocityModel3D):
    fig, ax = plt.subplots()
    ax.invert_yaxis()
    ax.set_ylabel("Depth (m)")
    ax.set_xlabel("Offset along x (m)")
    for depth in velocity_model.interface_depths:
        ax.axhline(depth, color="k")
    for segment in ray.path:
        x, z = segment.T[Index.X], segment.T[Index.Z]
        ax.plot(x, z)
    plt.show()


class IVPResultStatus(enum.IntEnum):
    """Enum to make the int status returned from solve_ivp more readable"""
    FAILED = -1
    END_REACHED = 0
    TERMINATION_EVENT = 1


ODEState3D = namedtuple("ODEState3D", ["x", "y", "z", "px", "py", "pz", "T"])

IVPEventFunction = Callable[[float, ODEState3D], float]


class Index(enum.IntEnum):
    X = 0
    Y = 1
    Z = 2


class NumericRayTracer3D:

    def __init__(self, velocity_model: VelocityModel3D):
        self.model = velocity_model
        self.layer = None

    @staticmethod
    def _velocity(layer: LinearVelocityLayer, x: float, y: float,
                  z: float) -> float:
        """
        Evaluate velocity at a given depth in a layer with linearly varying
        velocity
        :param layer: The layer to evaluate
        """
        return layer["intercept"] + layer["gradient"] * z

    def _trace(self, s: float, y: ODEState3D) -> ODEState3D:
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
        dpxds = derivative((lambda x_: 1 / self._velocity(self.layer, x_, y, z)),
                           x, dx=0.0001)
        dpyds = derivative((lambda y_: 1 / self._velocity(self.layer, x, y_, z)),
                           y, dx=0.0001)
        dpzds = derivative((lambda z_: 1 / self._velocity(self.layer, x, y, z_)),
                           z, dx=0.0001)
        dTds = 1. / v
        return ODEState3D(dxds, dyds, dzds, dpxds, dpyds, dpzds, dTds)

    def _trace_layer(self, ray: Ray3D, initial_slownesses: np.ndarray,
                     max_step_s: float) -> Ray3D:
        """
        Trace a ray through a single layer of the model
        :param ray: Ray to trace
        :param initial_slownesses: Initial value of px, py, pz (slowness along
        the corresponding axis) of the ray in s/m.
        :param max_step_s: Max step the solver takes for the integration
        variable s.
        :return: Ray with its parameters set
        """
        # workaround: events only activate after a step is taken because
        # sometimes events trigger directly after tracing starts and the
        # depth doesn't change enough. Keep event function artificially away
        # from zero to avoid that
        upper_layer_event: IVPEventFunction = lambda s, y_: y_[2] - self.layer["top_depth"] if s > max_step_s else 1
        lower_layer_event: IVPEventFunction = lambda s, y_: y_[2] - self.layer["bot_depth"] if s > max_step_s else -1
        for func in (upper_layer_event, lower_layer_event):
            func.terminal = True
        # TODO change signature of _velocity to accept np array instead of
        #  three floats
        initial_state = ODEState3D(*ray.last_point, *initial_slownesses,
                                   ray.last_time)
        result = solve_ivp(self._trace, (0, np.inf), initial_state,
                           max_step=max_step_s, events=(upper_layer_event,
                                                        lower_layer_event))
        x, y, z, px, py, pz, t = result.y
        ray.path.append(np.vstack((x, y, z)).T)
        ray.slowness.append(np.vstack((px, py, pz)).T)
        ray.travel_time.append(t)
        return ray

    def trace_stack(self, ray: Ray3D, ray_code: str = None,
                    max_step: float = 1) -> Ray3D:
        """
        Trace ray through a stack of layers. The ray type at an interface is
        chosen by the ray code.
        :param ray: Ray to trace through the model
        :param ray_code: Specifies which ray (Transmitted/Reflected) to follow
        at an interface in the model. "T" stands for the transmitted ray, "R"
        for the reflected ray. If not given, the ray will only be traced through
        the layer in which its starting point resides.
        :param max_step: Max step s for the integration.
        :return: Ray with parameters set
        """
        index = self.model.layer_index(ray.start[Index.Z])
        self.layer = self.model[index]
        initial_slownesses = initial_slowness3D(ray, self._velocity(self.layer, *ray.last_point))
        ray = self._trace_layer(ray, initial_slownesses, max_step)
        if ray_code is None:
            return ray
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
            ray = self._trace_layer(ray, new_p, max_step)
        return ray


def main():
    layers = [(0, 100, 1800, 4), (100, 200, 2400, 0), (200, 300, 2400, 1),
              (300, 400, 2700, 0), (400, 500, 2250, 1.5)]
    vm = VelocityModel3D(layers)
    ray = Ray3D(0, 0, 0, radians(20), radians(0))
    nrt = NumericRayTracer3D(vm)
    ray = nrt.trace_stack(ray, "TTTTRTTTT")
    plot_ray_in_model(ray, vm)


if __name__ == '__main__':
    main()
