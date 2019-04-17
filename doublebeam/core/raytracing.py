from math import sin, cos, radians, isclose
from typing import Tuple

import matplotlib.pyplot as plt
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


if __name__ == '__main__':
    start_x, start_z = 0, 0
    initial_angle_degrees = 45
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
    while s < 1:
        x += vm.eval_at(z) * px * ds
        z += vm.eval_at(z) * pz * ds
        px -= (1 / vm.eval_at(z)**2) * dvx() * ds
        pz -= (1 / vm.eval_at(z)**2) * dvz(vm, z) * ds
        s += ds
        points.append((x, z))

    x, z = zip(*points)
    plt.scatter(x, z)
    ax = plt.gca()
    # invert y axis so positive depth values are shown downwards
    ax.invert_yaxis()
    # set aspect ratio to equal so angles stay true
    ax.set_aspect("equal")
    plt.xlabel("x (km)")
    plt.ylabel("z (km)")
    plt.show()