from math import sin, cos, radians, isclose, copysign
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import scipy.misc


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

    def __init__(self, pierce_point_xm: float, theta: float):
        """
        :param pierce_point_xm: x coordinate of start point of ray in km
        :param theta: angle of ray against vertical at the surface in rad
        """
        self.xm = pierce_point_xm
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


def snells_law(px: float, pz: float, z: float, z_prev: float, wave_type="T"):
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
    n = np.array((0, 1))
    eps = copysign(1, np.dot(p, n))
    # calculate slowness vector after interface using eq. 2.4.70 from
    # Cerveny - Seismic ray theory
    v = vm.eval_at(z)
    v_prev = vm.eval_at(z_prev)
    pn = np.dot(p, n)
    minusplus = -1 if wave_type == "T" else 1
    #TODO not stable when using reflected wave
    p_new = p - n * (pn + minusplus * eps * (v**-2 - v_prev**-2 + pn**2)**0.5)
    return p_new


if __name__ == '__main__':
    start_x, start_z = 0, 0
    #TODO stability of transmission not given for angles > 30°
    initial_angle_degrees = 29
    ray = Ray2D(start_x, radians(initial_angle_degrees))
    vm = MockVelocityModel1D(1)
    V0 = vm.eval_at(start_z)
    px0 = calc_px(V0, ray.theta)
    pz0 = calc_pz(V0, ray.theta)

    ds = 0.01

    x = ray.xm
    z = start_z
    px = px0
    pz = pz0
    s = 0
    points = []
    points.append((start_x, start_z))
    while s < 9 and z >= 0:
        x += vm.eval_at(z) * px * ds
        z += vm.eval_at(z) * pz * ds
        if vm.boundary_crossed(points[-1][1], z):
            z_prev = points[-1][1]
            px, pz = snells_law(px, pz, z, z_prev)
        else:
            px -= (1 / vm.eval_at(z)**2) * dvx() * ds
            pz -= (1 / vm.eval_at(z)**2) * dvz(vm, z) * ds
        s += ds
        points.append((x, z))

    x, z = zip(*points)
    plt.plot(x, z)
    ax = plt.gca()
    # invert y axis so positive depth values are shown downwards
    ax.invert_yaxis()
    # set aspect ratio to equal so angles stay true
    ax.set_aspect("equal")
    plt.xlabel("x (km)")
    plt.ylabel("z (km)")
    plt.show()