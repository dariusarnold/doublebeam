from math import sin, cos, asin, acos, radians, degrees, isclose, copysign, sqrt
from cmath import sqrt as csqrt
import itertools
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import scipy.misc
import scipy.integrate


def cartesian_to_ray_s(x, z, xm, _theta):
    """Conversion of 2D cartesian coordinates x, z to ray coordinate
    s, see eq. 14a in Hill1990 Gaussian beam migration"""
    return (x - xm) * sin(_theta) - z * cos(_theta)


def cartesian_to_ray_n(x, z, xm, theta):
    """Conversion of 2D cartesian coordinates x, z to ray coordinate
    n, see eq. 14b in Hill1990 Gaussian beam migration"""
    return (x - xm) * cos(theta) + z * sin(theta)


class MockVelocityModel1D:

    def __init__(self, boundary_depth_km):
        self.boundary_depth_km = boundary_depth_km

    def eval_at(self, z):
        """Mock 1D velocity model. Get velocity as a function of depth.
        Can be replaced with VelocityModel1D when layers with linear velocity
        change are implemented"""
        if z < self.boundary_depth_km:
            v0 = 3
            slope = 1
        else:
            v0 = 4.5  # km/s
            slope = 1.5  # km/s per km depth
        return v0 + z * slope

    def at_boundary(self, z, abs_tol=0.01):
        return isclose(z, self.boundary_depth_km, abs_tol=abs_tol)

    def boundary_crossed(self, z1, z2):
        """Return true if a boundary is between the depths z1, z2"""
        # sort limits so that even when going upwards (where z1 > z2) a crossing is detected
        _z1 = min(z1, z2)
        _z2 = max(z1, z2)
        return _z1 <= self.boundary_depth_km <= _z2


def dvx():
    """Derivative of velocity after x. 0 For 1D model"""
    return 0


def dvz(velocity_model, z, delta=0.0001):
    """Derivative of velocity after z"""
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
    pn = np.dot(p, n)
    minusplus = -1 if wave_type == "T" else 1
    # TODO not stable when using reflected wave
    # calculate slowness vector after interface using eq. 2.4.70 from
    # Cerveny - Seismic ray theory
    p_new = p - n * (pn + minusplus * eps * sqrt(v**-2 - v_prev**-2 + pn**2))
    return p_new


def plot_rays(points1, points2):
    x1, z1 = zip(*points1)
    x2, z2 = zip(*points2)
    plt.plot(x1, z1, label="My integration")
    plt.plot(x2, z2, label="Scipy")
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
    x, z, px, pz = y
    v = vm1.eval_at(z)
    dxds = v * px
    dzds = v * pz
    dpxds = -1 * v**-2 * dvx()
    dpzds = -1 * v**-2 * dvz(vm1, z)
    dydt = [dxds, dzds, dpxds, dpzds]
    return dydt


def ray_trace_euler(ray, velocity_model, s_end, ds=0.01):
    s = 0
    x = ray.x0
    z = ray.z0
    v0 = velocity_model.eval_at(z)
    px = calc_px(v0, ray.theta)
    pz = calc_pz(v0, ray.theta)
    points = [(x, z)]
    while s < s_end and z >= 0:
        v = velocity_model.eval_at(z)
        x += v * px * ds
        z += v * pz * ds
        if velocity_model.boundary_crossed(points[-1][1], z):
            z_prev = points[-1][1]
            px, pz = snells_law(px, pz, z, z_prev, velocity_model)
        else:
            px -= (1 / v**2) * dvx() * ds
            pz -= (1 / v**2) * dvz(velocity_model, z) * ds
        s += ds
        points.append((x, z))
    return points


def ray_trace_scipy(ray, velocity_model, s_end, ds=0.01):
    # Function which has a zero crossing at the boundary depth
    crossed = lambda t, y: y[1] - velocity_model.boundary_depth_km
    # set to True to stop integration at the boundary so we can apply Snells law
    crossed.terminal = True

    # return z coordinate which has a natural zero crossing at the surface
    surfaced = lambda t, y: y[1]
    # set True to stop integration once the ray reaches the surface
    surfaced.terminal = True

    global vm1
    vm1 = velocity_model

    start_z = ray.z0
    V0 = velocity_model.eval_at(start_z)
    px0 = calc_px(V0, ray.theta)
    pz0 = calc_pz(V0, ray.theta)
    ds = 0.01

    min_float_step = np.finfo(float).eps
    scipy_x = []
    scipy_z = []
    # move z slightly below surface so event wont trigger immediately
    result = scipy.integrate.solve_ivp(trace, [0, s_end], [ray.x0, start_z+min_float_step, px0, pz0], max_step=ds, events=[crossed, surfaced])
    while result.status == 1:
        # stop integration when surface was reached
        if result.t_events[1].size > 0:
            break
        s_event = result.t_events[0][0]
        _x, _z, _px, _pz = result.y
        scipy_x.append(_x)
        scipy_z.append(_z)
        px, pz = snells_law(_px[-1], _pz[-1], _z[-1], _z[-2], velocity_model)
        # move z behind the interface just passed by the ray so the event wont
        # trigger again. Increase z (positive step) when ray goes down, decrease
        # z (negative step) when ray goes up
        min_float_step = copysign(min_float_step, _pz[-1])
        result = scipy.integrate.solve_ivp(trace, [s_event, s_end], [_x[-1], _z[-1]+min_float_step, px, pz], max_step=ds, events=[crossed, surfaced])
    scipy_x.append(result.y[0])
    scipy_z.append(result.y[1])
    # scipy_x and scipy_z are lists of lists. Every sublist contains part of a
    # ray path between interfaces. To plot them, unpack the inner lists
    scipy_x = list(itertools.chain.from_iterable(scipy_x))
    scipy_z = list(itertools.chain.from_iterable(scipy_z))
    # for plotting, create
    scipy_points = list(zip(scipy_x, scipy_z))
    return scipy_points



if __name__ == '__main__':
    start_x, start_z = 0, 0
    initial_angle_degrees = 20
    ray = Ray2D(start_x, start_z, radians(initial_angle_degrees))
    boundary_depth_km = 1
    vm = MockVelocityModel1D(boundary_depth_km)
    s_end = 20

    points = ray_trace_euler(ray, vm, s_end)
    scipy_points = ray_trace_scipy(ray, vm, s_end)

    plot_rays(points, scipy_points)
