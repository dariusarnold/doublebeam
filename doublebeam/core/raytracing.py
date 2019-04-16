from math import sin, cos, radians

import matplotlib.pyplot as plt


def cartesian_to_ray_s(x, z, xm, _theta):
    """Conversion of 2D cartesian coordinates x, z to ray coordinate
    s, see eq. 14a in Hill1990 Gaussian beam migration"""
    return (x - xm) * sin(_theta) - z * cos(_theta)


def cartesian_to_ray_n(x, z, xm, theta):
    """Conversion of 2D cartesian coordinates x, z to ray coordinate
    n, see eq. 14b in Hill1990 Gaussian beam migration"""
    return (x - xm) * cos(theta) + z * sin(theta)


def velocity(z):
    """Mock 1D velocity model. Get velocity as a function of depth.
    Can be replaced with VelocityModel1D when layers with linear velocity
    change are implemented"""
    boundary_depth_km = 0.4
    if z < boundary_depth_km:
        v0 = 3
        slope = 0.75
    else:
        v0 = 4.5 # km/s
        slope = 0.5 # km/s per km depth
    return v0 + z * slope

def dvx(x, z, delta=0.0001):
    """Derivative of velocity after x"""
    return (velocity(z) - velocity(z)) / delta


def dvz(x, z, delta=0.0001):
    return (velocity(z+delta) - velocity(z-delta)) / delta

def diff(func, x, z, delta):
    # broken
    return (func)


def calc_px(v, theta):
    return sin(theta) / v


def calc_pz(v, theta):
    return cos(theta) / v


class Ray2D:

    def __init__(self, pierce_point_xm: float, theta: float):
        """
        :param pierce_point_xm: x coordinate of pierce point of ray at surface in km
        :param theta: angle of ray against vertical at the surface in rad
        """
        self.xm = pierce_point_xm
        self.theta = theta


if __name__ == '__main__':
    start_x, start_z = 0, 0
    ray = Ray2D(start_x, radians(45))
    V0 = velocity(start_z)
    px0 = calc_px(V0, ray.theta)
    pz0 = calc_pz(V0, ray.theta)

    ds = 0.01

    x = ray.xm
    z = 0
    px = px0
    pz = pz0
    theta = ray.theta
    s = 0
    points = []
    s = cartesian_to_ray_s(x, z, ray.xm, ray.theta)
    while s < 1:
        x += velocity(z) * px * ds
        z += velocity(z) * pz * ds
        px -= (1 / velocity(z)**2) * dvx(x, z) * ds
        px -= (1 / velocity(z)**2) * dvz(x, z) * ds
        s += ds
        points.append((x, z))

    x, z = zip(*points)
    plt.scatter(x, z)
    plt.show()