from math import sin, cos, radians

import matplotlib.pyplot as plt


def cartesion_to_ray_s(x, z, xm,_theta):
    return (x - xm) * sin(_theta) - z * cos(_theta)

def cartesion_to_ray_n(x, z, xm, theta):
    return (x - xm) * cos(theta) + z * sin(theta)

def velocity(x, z):
    v0 = 4.5 # km/s
    slope = 0.5 # km/s per km depth
    return v0 + z * slope

def dvx(x, z, delta=0.0001):
    """Derivative of velocity after x"""
    return (velocity(x+delta, z) - velocity(x-delta, z)) / delta

def dvz(x, z, delta=0.0001):
    return (velocity(x, z+delta) - velocity(x, z-delta)) / delta

def diff(func, x, z, delta):
    # broken
    return (func)

def calc_px(v, theta):
    return sin(theta) / v

def calc_pz(v, theta):
    return cos(theta) / v

class Ray2D:

    def __init__(self, pierce_point_xm: float, theta: float):
        # x coordinate of pierce point of ray at surface in km
        self.xm = pierce_point_xm
        # angle of ray against vertical at the surface in rad
        self.theta = theta

if __name__ == '__main__':
    ray = Ray2D(0, radians(45))
    V0 = velocity(ray.xm, 0)
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
    s = cartesion_to_ray_s(x, z, ray.xm, ray.theta)
    while s < 1:
        x += velocity(x, z) * px * ds
        z += velocity(x, z) * pz * ds
        px -= (1 / velocity(x, z)**2) * dvx(x, z) * ds
        px -= (1 / velocity(x, z)**2) * dvz(x, z) * ds
        s += ds
        points.append((x, z))

    x, z = zip(*points)
    plt.scatter(x, z)
    plt.show()