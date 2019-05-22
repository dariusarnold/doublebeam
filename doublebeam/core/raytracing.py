import enum
import itertools
from math import sin, cos, asin, acos, radians, degrees, copysign, sqrt
from typing import Tuple, Sequence

import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.integrate
import scipy.misc

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
    return 0


def dvz(velocity_model, z, delta=0.0001):
    """Derivative of velocity after z"""
    if z < delta:
        # special case: evaluating derivative would eval model at negative depth
        return 0
    return scipy.misc.derivative(velocity_model.eval_at, z, delta)


def calc_px(v, theta):
    """Calculate horizontal slowness"""
    return sin(theta) / v


def calc_pz(v, theta):
    """Calculate vertical slowness"""
    return cos(theta) / v


class Ray2D:

    def __init__(self, pierce_point_xm: float, start_z: float, theta: float):
        """
        :param pierce_point_xm: x coordinate of start point of ray in km
        :param theta: angle of ray against vertical at the surface in rad
        """
        self.x0 = pierce_point_xm
        self.z0 = start_z
        self.theta = theta
        self.layer_boundaries_crossed_depths = [-99]

    def get_last_boundary_crossed_depth(self):
        try:
            return self.layer_boundaries_crossed_depths[-1]
        except IndexError:
            return None


class Ray3D:

    def __init__(self, pierce_point_xm: float, pierce_point_ym: float, theta: float, phi: float):
        """
        :param pierce_point_xm: x coordinate of start point of ray in km
        :param pierce_point_ym: y coordinate of start point of ray in km
        :param theta: Angle against vertical at start point in rad
        0 <= theta <= pi
        :param phi: Angle against x axis at start point in rad, with increasing
        angle towards the y axis
        0 <= phi <= 2*pi
        """
        self.xm = pierce_point_xm
        self.ym = pierce_point_ym
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


def length(vector) -> float:
    """Calculate length of vector"""
    return np.dot(vector, vector)**0.5


def angle(vector1, vector2) -> float:
    """Calculate angle between two vectors"""
    return acos(np.dot(vector1, vector2) / (length(vector1) * length(vector2)))


def critical_angle(v1, v2):
    """
    Use Snells law to calculate the critical angle at an interface
    :param v1: velocity before interface
    :param v2: velocity after interface
    :return: critical angle in rad
    """
    if v1 < v2:
        return asin(v1/v2)
    return np.pi


def snells_law(px: float, pz: float, z: float, z_prev: float, velocity_model, wave_type="T"):
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


def plot_ray(points: Sequence[Tuple[float, float]]):
    x1, z1 = zip(*points)
    plt.plot(x1, z1, label="Ray path")
    ax = plt.gca()
    # invert y axis so positive depth values are shown downwards
    ax.invert_yaxis()
    # set aspect ratio to equal so angles stay true
    ax.set_aspect("equal")
    plt.xlabel("x (km)")
    plt.ylabel("z (km)")
    plt.legend()
    plt.show()


def trace(t, y):
    """
    Standard raypath equations in 2D from Hill1990 Gaussian beam migration eq. 2a-2d
    :param t:
    :param y:
    :return:
    """
    x, z, px, pz = y
    v = vm.eval_at(z)
    dxds = v * px
    dzds = v * pz
    dpxds = -1 * v**-2 * dvx()
    dpzds = -1 * v**-2 * dvz(vm, z)
    dydt = [dxds, dzds, dpxds, dpzds]
    return dydt


class IVPResultStatus(enum.IntEnum):
    """Enum class to make the int status returned from scipys solve_ivp more readable"""
    FAILED = -1
    END_REACHED = 0
    TERMINATION_EVENT = 1


def ray_trace_scipy(ray: Ray2D, velocity_model, s_end: float, ds: float = 0.01) -> Sequence[Tuple[float, float]]:
    # Generate a function which has a zero crossing at the boundary depth
    # for all boundary depths in the models to apply Snells law
    crossings = [lambda t, y: y[1] - interface_depth for interface_depth in velocity_model.interface_depths[1:-1]]  # skip first interface depth since its the surface at z = 0 and skip last interface since model stops there # TODO add condition that stops integration once model bottom is reached
    # set to True to stop integration at the boundary so we can apply Snells law
    for function in crossings:
        function.terminal = True

    # return z coordinate which has a natural zero crossing at the surface
    surfaced = lambda t, y: y[1]
    # set True to stop integration once the ray reaches the surface
    surfaced.terminal = True

    # This is a workaround so that trace can access the velocity model
    # TODO find better solution
    global vm
    vm = velocity_model

    z0 = ray.z0
    v0 = velocity_model.eval_at(z0)
    px0 = calc_px(v0, ray.theta)
    pz0 = calc_pz(v0, ray.theta)

    min_float_step = np.finfo(float).eps
    x_values = []
    z_values = []
    # move z slightly below surface so event wont trigger immediately
    result = sp.integrate.solve_ivp(trace, [0, s_end], [ray.x0, z0+min_float_step, px0, pz0],
                                    max_step=ds, events=[*crossings, surfaced])  # type: scipy.integrate._ivp.ivp.OdeResult
    while result.status == IVPResultStatus.TERMINATION_EVENT:
        # events are returned in the order as passed to solve_ivp
        crossing_events, surface_events = result.t_events
        # stop integration when surface was reached
        if surface_events.size > 0:
            break
        s_event = crossing_events[0]
        _x, _z, _px, _pz = result.y
        x_values.append(_x)
        z_values.append(_z)
        px, pz = snells_law(_px[-1], _pz[-1], _z[-1], _z[-2], velocity_model)
        # move z behind the interface just passed by the ray so the event wont
        # trigger again. Increase z (positive step) when ray goes down, decrease
        # z (negative step) when ray goes up
        min_float_step = copysign(min_float_step, _pz[-1])
        result = scipy.integrate.solve_ivp(trace, [s_event, s_end], [_x[-1], _z[-1]+min_float_step, px, pz], max_step=ds, events=[*crossings, surfaced])
    x_values.append(result.y[0])
    z_values.append(result.y[1])
    # scipy_x and scipy_z are lists of lists. Every sublist contains part of a
    # ray path between interfaces. To plot them, unpack the inner lists
    x_values = list(itertools.chain.from_iterable(x_values))
    z_values = list(itertools.chain.from_iterable(z_values))
    # create a list of (x, z) tuples where every tuple is one point of the ray
    return list(zip(x_values, z_values))


if __name__ == '__main__':
    start_x, start_z = 0, 0
    initial_angle_degrees = 20
    ray = Ray2D(start_x, start_z, radians(initial_angle_degrees))
    vm = VelocityModel1D.from_string("0, 1, 3, 4, 0, 0, 1, 1\n1, 101, 6, 156, 0, 0, 1, 1")
    s_end = 20
    points = ray_trace_scipy(ray, vm, s_end)
    plot_ray(points)
