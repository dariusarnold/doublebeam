from typing import List, Union, Tuple

import numpy as np

from doublebeam.utils import Index, slowness_3D


class Ray3D:

    def __init__(self, start_coordinates: Union[np.ndarray, Tuple[float, float, float]],
                 slowness: np.ndarray):
        """
        :param start_coordinates: x, y, z coordinate triple of start point of
        ray in m
        :param slowness: Array of three floats: px, py, pz
        """
        self.start = np.asarray(start_coordinates)
        # These will be set after the ray is traced
        self.path: List[np.ndarray] = []
        slowness = np.reshape(slowness, (1, 3))
        self.slowness: List[np.ndarray] = [slowness]
        self.travel_time: List[np.ndarray] = []

    @classmethod
    def from_angle(cls, start_coordinates: Union[np.ndarray, Tuple[float, float, float]],
                   theta: float, phi: float, velocity: float) -> "Ray3D":
        """
        Create a ray from angles and calculate slowness values.
        :param start_coordinates: x, y, z coordinate triple of start point of
        ray in m
        :param theta: Angle against downgoing vertical axis (z) at start
        point in rad, increasing upwards. 0 <= theta <= pi
        :param phi: Angle against x axis at start point in rad, with increasing
        angle towards the y axis
        0 <= phi <= 2*pi
        :param velocity: Velocity in m/s at start point of the ray, used to
        calculate slowness
        :return:
        """
        p = slowness_3D(theta, phi, velocity)
        return cls(start_coordinates, p)


    @property
    def last_point(self) -> np.ndarray:
        """
        Return last point of the rays path. If the ray hasn't been traced yet,
        return the starting point of the ray.
        """
        try:
            return self.path[-1][-1]
        except IndexError:
            return self.start

    @property
    def last_slowness(self) -> Union[np.ndarray, None]:
        """
        Return last slowness value. If the ray hasn't been traced yet, return
        None.
        """
        try:
            return self.slowness[-1][-1]
        except IndexError:
            return None

    @property
    def last_time(self) -> float:
        """
        Return last travel time of the ray.
        :return:
        """
        try:
            return self.travel_time[-1][-1]
        except IndexError:
            return 0.

    @property
    def direction(self) -> str:
        """
        Return direction of the ray along the z-axis. A downgoing ray will
        return "down", while an upgoing ray will return "up". The special case
        of a horizontal ray (vertical slowness pz = 0) will return "horizontal".
        """
        pz = self.last_slowness[Index.Z]
        if abs(pz) < np.finfo(type(pz)).eps:
            return "horizontal"
        if pz < 0:
            return "up"
        elif pz > 0:
            return "down"
