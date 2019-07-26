from typing import List, Union

import numpy as np

from doublebeam.core.utils import Index


class Ray3D:

    def __init__(self, start_x: float, start_y: float, start_z: float,
                 theta: float, phi: float):
        """
        :param start_x: x coordinate of start point of ray in m
        :param start_y: y coordinate of start point of ray in m
        :param start_z: z coordinate of start point of ray in m
        :param theta: Angle against downgoing vertical axis (z) at start
        point in rad, increasing upwards. 0 <= theta <= pi
        :param phi: Angle against x axis at start point in rad, with increasing
        angle towards the y axis
        0 <= phi <= 2*pi
        """
        self.start = np.array((start_x, start_y, start_z))
        self.theta = theta
        self.phi = phi
        # These will be set after the ray is traced
        self.path: List[np.ndarray] = []
        self.slowness: List[np.ndarray] = []
        self.travel_time: List[np.ndarray] = []

    @classmethod
    # TODO refactor this and init to take tuple or array
    def from_slowness(cls, start_x, start_y, start_z, slowness_x: float,
                      slowness_y: float) -> "Ray3D":
        raise NotImplementedError

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
        try:
            pz = self.last_slowness[Index.Z]
        except TypeError:
            # if None is returned (ray hasn't been traced yet), use starting
            # angle against vertical axis
            pz = np.pi/2 - self.theta
        if pz < 0:
            return "up"
        elif pz > 0:
            return "down"
        return "horizontal"
