import enum

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