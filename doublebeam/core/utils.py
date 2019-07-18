import enum
from math import acos

import numpy as np

"""
Common and general functions and classes 
"""


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
    angle = acos((vector1 @ vector2) / (length(vector1) * length(vector2)))
    if acute:
        return angle
    else:
        return np.pi - angle


def horizontal_distance(point1: np.ndarray, point2: np.ndarray) -> float:
    """
    Calculate distance between two points projected on a horizontal plane.
    """
    vector = point1 - point2
    vector[Index.Z] = 0.
    return length(vector)
