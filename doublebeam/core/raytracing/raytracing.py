import enum
import time
from collections import namedtuple
from math import sin, cos, asin, acos, radians, degrees, copysign, sqrt
from typing import Tuple, Callable, List

import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate
import scipy.misc
from scipy.integrate import solve_ivp

from doublebeam.core.models import VelocityModel1D


def cartesian_to_ray_s(x, z, xm, _theta):
    """Conversion of 2D cartesian coordinates x, z to ray coordinate
    s, see eq. 14a in Hill1990 Gaussian beam migration"""
    return (x - xm) * sin(_theta) - z * cos(_theta)


def cartesian_to_ray_n(x, z, xm, theta):
    """Conversion of 2D cartesian coordinates x, z to ray coordinate
    n, see eq. 14b in Hill1990 Gaussian beam migration"""
    return (x - xm) * cos(theta) + z * sin(theta)


def dvx():
    """Derivative of velocity after x. 0 For 1D model"""
    return 0.


def dvz(velocity_model, z, delta=0.0001):
    """Derivative of velocity after z"""
    if z < delta:
        # special case: evaluating derivative would eval model at negative depth
        return 0.
    return scipy.misc.derivative(velocity_model.eval_at, z, delta)


def horizontal_slowness(v, theta):
    return sin(theta) / v


def vertical_slowness(v, theta):
    return cos(theta) / v


class Ray2D:

    def __init__(self, start_x: float, start_z: float, theta: float):
        """
        :param start_x: x coordinate of start point of ray in m
        :param start_z: z coordinate of start point of ray in m
        :param theta: angle of ray against vertical at the surface in rad
        """
        self.x0 = start_x
        self.z0 = start_z
        self.theta = theta
        # Stores x, y coordinates of ray path
        self.path = None  # type: Tuple[np.array, np.array]
        self.layer_boundaries_crossed_depths = [-99]

    def get_last_boundary_crossed_depth(self):
        try:
            return self.layer_boundaries_crossed_depths[-1]
        except IndexError:
            return None


class Ray3D:

    def __init__(self, start_x: float, start_y: float, start_z: float, theta: float, phi: float):
        """
        :param start_x: x coordinate of start point of ray in m
        :param start_y: y coordinate of start point of ray in m
        :param start_z: z coordinate of start point of ray in m
        :param theta: Angle against vertical at start point in rad
        0 <= theta <= pi
        :param phi: Angle against x axis at start point in rad, with increasing
        angle towards the y axis
        0 <= phi <= 2*pi
        """
        self.x0 = start_x
        self.y0 = start_y
        self.z0 = start_z
        self.theta = theta
        self.phi = phi


def calc_initial_slowness3D(ray: Ray3D, v0: float) -> Tuple[float, float, float]:
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
    return px, py, pz


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


def snells_law(px: float, pz: float, z: float, z_prev: float, velocity_model: VelocityModel1D, wave_type: str = "T"):
    """
    Compute new slowness vector using Snells law
    :param px: x component of slowness vector before interface
    :param pz: z component of slowness vector before interface
    :param z: depth after interface
    :param z_prev: depth before interface
    :param wave_type: Select for which wave type (reflected/transmitted) the
    slowness should be computed. "T" for transmitted, "R" for reflected.
    :return: slowness vector (x, z) after interface
    """
    # create slowness vector p
    # TODO keep p as a np.ndarray everywhere
    p = np.array((px, pz))
    # normal vector in 1D model will always be vertical
    # n should be oriented to the side where the transmitted wave propagates
    # else the minus/plus relation for transmitted/reflected waves isn't valid
    n = np.array((0, copysign(1, pz)))
    # look at angle to determine critical incidence
    angle_in = angle(p, n)
    print("Angle in", degrees(angle_in))
    eps = copysign(1, np.dot(p, n))
    v = velocity_model.eval_at(z)
    v_prev = velocity_model.eval_at(z_prev)
    angle_crit = critical_angle(v_prev, v)
    print("Critical angle", degrees(angle_crit))
    if angle_in > angle_crit:
        print("Total reflection")
        p[1] *= -1
        return p
    #angle_out = asin(sin(angle_in) * v/v_prev)
    #print("Angle out", degrees(angle_out))
    pn = np.dot(p, n)
    minusplus = -1 if wave_type == "T" else 1
    # TODO not stable when using reflected wave
    # calculate slowness vector after interface using eq. 2.4.70 from
    # Cerveny - Seismic ray theory
    p_new = p - n * (pn + minusplus * eps * sqrt(v**-2 - v_prev**-2 + pn**2))
    return p_new


def plot_ray(ray: Ray2D):
    x, z = ray.path
    plt.plot(x, z, label="Ray path")
    ax = plt.gca()
    # invert y axis so positive depth values are shown downwards
    ax.invert_yaxis()
    # set aspect ratio to equal so angles stay true
    ax.set_aspect("equal")
    plt.xlabel("x (km)")
    plt.ylabel("z (km)")
    plt.legend()
    plt.show()





class IVPResultStatus(enum.IntEnum):
    """Enum class to make the int status returned from scipys solve_ivp more readable"""
    FAILED = -1
    END_REACHED = 0
    TERMINATION_EVENT = 1


ODEState = namedtuple("ODEState", ["x", "z", "px", "pz"])

IVPEventFunction = Callable[[float, ODEState], float]

class RayTracer2D:
    """
    Class for ray tracing in a 2D velocity model.
    """

    def __init__(self, velocity_model: VelocityModel1D):
        """
        :param velocity_model: Velocity model to use for ray tracing
        """
        self._velocity_model = velocity_model
        # Generate a function which has a zero crossing at the boundary depth
        # for all boundary depths in the models to apply Snells law
        # skip first interface depth since its the surface at z = 0 and skip
        # last interface since model stops there
        #  TODO add condition that stops integration once model bottom is reached
        crossings : List[IVPEventFunction] = [lambda s, y: y[1] - depth
            for depth in velocity_model.interface_depths[1:-1]]
        # set to True to stop integration at the boundary so we can apply Snells law
        for f in crossings:
            f.terminal = True
        # return z coordinate which has a natural zero crossing at the surface
        surfaced: IVPEventFunction = lambda s, y: y[1]
        # set True to stop integration once the ray reaches the surface
        surfaced.terminal = True
        self._events = [*crossings, surfaced]

    def _trace(self, s: float, y: ODEState) -> ODEState:
        """
        Standard raypath equations in 2D from Hill1990 Gaussian beam migration eq. 2a-2d
        :param s: Current value of integration variable s (arclength along ray)
        :param y: List of values for x, y coordinate, horizontal and vertical slowness
        :return:
        """
        x, z, px, pz = y
        v = self._velocity_model.eval_at(z)
        dxds = v * px
        dzds = v * pz
        # TODO simplify dvx by replacing it with zero for 2D case
        dpxds = -1 * v ** -2 * dvx()
        dpzds = -1 * v ** -2 * dvz(self._velocity_model, z)
        dydt = ODEState(dxds, dzds, dpxds, dpzds)
        return dydt

    def ray_trace(self, ray: Ray2D, s_end: float, ds: float=0.01) -> Ray2D:
        v0 = self._velocity_model.eval_at(ray.z0)
        px0 = horizontal_slowness(v0, ray.theta)
        pz0 = vertical_slowness(v0, ray.theta)
        min_float_step = np.finfo(float).eps
        initial_state = ODEState(ray.x0, ray.z0+min_float_step, px0, pz0)
        result: scipy.integrate._ivp.ivp.OdeResult = solve_ivp(self._trace, [0, s_end], initial_state, max_step=ds,
                           events=self._events)
        x_values, z_values = [], []
        while result.status == IVPResultStatus.TERMINATION_EVENT:
            # events are returned in the order as passed to solve_ivp
            crossing_events, surface_events = result.t_events
            # stop integration when surface was reached
            if surface_events.size > 0:
                break
            s_event = crossing_events[0]
            x_, z_, px_, pz_ = result.y
            x_values.append(x_)
            z_values.append(z_)
            px, pz = snells_law(px_[-1], pz_[-1], z_[-1], z_[-2], self._velocity_model)
            min_float_step = copysign(min_float_step, pz_[-1])
            # move z behind the interface just passed by the ray so the event wont
            # trigger again. Increase z (positive step) when ray goes down, decrease
            # z (negative step) when ray goes up
            result = solve_ivp(self._trace, [s_event, s_end],
                               ODEState(x_[-1], z_[-1]+min_float_step, px, pz),
                               max_step=ds, events=self._events)
        x_values.append(result.y[0])
        z_values.append(result.y[1])
        x_values = np.concatenate(x_values)
        z_values = np.concatenate(z_values)
        ray.path = (x_values, z_values)
        return ray


if __name__ == '__main__':
    start_x, start_z = 0, 0
    initial_angle_degrees = 20
    ray = Ray2D(start_x, start_z, radians(initial_angle_degrees))
    vm = VelocityModel1D.from_string("0, 1, 3, 4, 0, 0, 1, 1\n1, 101, 6, 156, 0, 0, 1, 1")
    s_end = 20
    ray_tracer = RayTracer2D(vm)
    a = time.time()
    ray = ray_tracer.ray_trace(ray, s_end)
    b = time.time()
    print(b-a)
    plot_ray(ray)
