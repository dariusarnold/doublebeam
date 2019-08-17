import enum
from math import acos, cos, sin, atan2, radians
from typing import Tuple

import numpy as np

"""
Common and general functions and classes 
"""


# constant specifying to how many digits a number should be rounded to remedy
# numeric issues where a float falls slightly outside of its valid range.
# Lower this because 13 still had issues
DIGITS_PRECISION = 12


class Index(enum.IntEnum):
    """
    Class that improves indexing readability by replacing "magic" 0,1,2 by names
    """
    X = 0
    Y = 1
    Z = 2


def length(vector: np.ndarray) -> float:
    """Calculate length of vector"""
    return (vector @ vector)**0.5


def angle(vector1: np.ndarray, vector2: np.ndarray, acute: bool = True) -> float:
    """
    Calculate angle between two vectors in rad.
    :param vector1: N dimensional vector
    :param vector2: N dimensional vector
    :param acute: If true, return the acute angle, else return the obtuse angle.
    This is not interpreted to be the reflex angle.
    :return: Angle between vector1 and vector2 in rad
    """
    angle = acos(np.clip((vector1 @ vector2) / (length(vector1) * length(vector2)), -1, 1))
    if acute:
        return angle
    else:
        return np.pi - angle


def angle_clockwise(vector1: np.ndarray, vector2: np.ndarray) -> float:
    """
    Calculate the clockwise angle from the first vector to the second vector in
    the x-y horizontal plane.
    :return: Clockwise angle from vector1 to vector2 in radians
    """
    angle1 = atan2(*vector1[0:2])
    angle2 = atan2(*vector2[0:2])
    difference = angle2 - angle1
    if difference < 0:
        # special case in the bottom right quadrant
        return 2 * np.pi + difference
    return difference


def horizontal_distance(point1: np.ndarray, point2: np.ndarray) -> float:
    """
    Calculate distance between two points projected on a horizontal plane.
    """
    vector = point1 - point2
    vector[Index.Z] = 0.
    return length(vector)


def safe_divide(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    Divide a by b but return 0 on places where b is zero without dividing
    """
    return np.divide(a, b, out=np.zeros_like(a), where=b != 0)


def slowness_3D(theta: float, phi: float, velocity: float) -> np.ndarray:
    """
    Calculate initial vertical and horizontal slowness for a ray.
    For geometric definitions see chapter 3.2.1 in Cerveny - Seismic ray theory
    :param theta: Angle against downgoing vertical axis (z) at start
    point in rad, increasing upwards. 0 <= theta <= pi
    :param phi: Angle against x axis at start point in rad, with increasing
    angle towards the y axis
    0 <= phi <= 2*pi
    :param velocity: Velocity in m/s at start point of the ray, used to
    calculate slowness
    :return: Tuple of slowness values
    """
    px = 1/velocity * sin(theta) * cos(phi)
    py = 1/velocity * sin(theta) * sin(phi)
    pz = 1/velocity * cos(theta)
    return np.array((px, py, pz))


def generate_grid_coordinates(target_depth: float, x_extent: Tuple[float, float],
                              y_extent: Tuple[float, float], num_x: int,
                              num_y: int) -> np.ndarray:
    """
    Generate coordinates of a evenly spaced grid.
    :param x_extent: x_start, x_end tuple.
    :param y_extent: y_start, y_end tuple.
    :param target_depth: Depth of fractures in m.
    :param num_x: Number of targets along x-axis.
    :param num_y: Number of targets along y-axis.
    :return: Array containing all target locations with shape (num_x, num_y, 3)
    where the first axis is the x axis, the second axis the y axis and the third
    axis contains the z value.
    """
    x = np.linspace(*x_extent, num_x)
    y = np.linspace(*y_extent, num_y)
    xx, yy = np.meshgrid(x, y)
    return np.ascontiguousarray(np.stack((xx, yy, np.full_like(xx, target_depth))).T)


def unit_vector(v: np.ndarray) -> np.ndarray:
    """
    Create unit vector.
    :param v: Vector specifying direction.
    :return: Same vector with length 1.
    """
    return v / length(v)


def pairwise(iterable):
    """
    Create pairwise iterator from iterable which returns pairs of values.
    "s -> (s0, s1), (s2, s3), (s4, s5), ..."
    :param iterable:
    """
    a = iter(iterable)
    return zip(a, a)


def generate_vector_arc(num_values: int, central_direction: np.ndarray):
    """
    Generate num_values of unit vectors representing fracture orientation.
    The range of vectors returned will have an angle from 0° to 180° increasing
    in clockwise direction. The first vector will have a clockwise angle of -90°
    or 270° degrees of the central direction, the last vector will have a
    clockwise angle of 90°. All vectors will be in the horizontal plane, and the
    axis will be projected on the horizontal plane as well.
    :param num_values: How many vectors to return.
    :param central_direction: Direction of the center of the vectors.
    :return: Array of shape (num_values, 3) containing the vectors.
    """
    x_axis = np.array((1., 0, 0))
    angles = np.linspace(0, radians(180), num_values)
    # -90° to rotate to the left
    angle_against_axis = angle_clockwise(x_axis, central_direction) - radians(90)
    angles += angle_against_axis
    vectors = np.zeros((num_values, 3))
    vectors[:, Index.X] = np.cos(angles)
    vectors[:, Index.Y] = -np.sin(angles)
    return vectors
