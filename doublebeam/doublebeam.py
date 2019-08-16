from typing import Optional

import numpy as np

from doublebeam.models import VelocityModel3D
from doublebeam.raytracing.twopoint import TwoPointRayTracing
from doublebeam.utils import generate_grid_coordinates, Index, unit_vector


class FractureParameters:

    def __init__(self, fracture_depth: float, fracture_spacing: float,
                 fracture_orientation: np.ndarray):
        """
        Aggregate class holding fracture information for vertical fracture
        planes.
        :param fracture_depth: Depth of top of fractures in m.
        :param fracture_spacing: Distance between fracture planes
        :param fracture_orientation: Vector orthogonal to fracture planes.
        """
        self.depth = fracture_depth
        self.spacing = fracture_spacing
        self.orientation = unit_vector(fracture_orientation)


class DoubleBeam:

    def __init__(self, model: VelocityModel3D, window_length: float,
                 fracture_parameters: FractureParameters):
        self.model = model
        self.fracture_info = fracture_parameters
        self.target_locations = None
        self.source_beam_centers: Optional[np.ndarray] = None
        self.window_length = window_length
        self.twopoint = TwoPointRayTracing(model)

    def generate_targets(self, num_x: int, num_y: int) -> None:
        self.target_locations = generate_grid_coordinates(self.fracture_info.depth, (0, self.model.x_width),
                                                          (0, self.model.y_width), num_x, num_y)

    def generate_source_beam_centers(self, num_x: int, num_y: int) -> None:
        self.source_beam_centers = generate_grid_coordinates(self.model.vertical_boundaries()[0], (0, self.model.x_width),
                                                             (0, self.model.y_width), num_x, num_y)

    def scattered_slowness(self, slowness: np.ndarray, phi_hat: np.ndarray,
                           fracture_spacing: float, frequency: float) -> np.ndarray:
        """
        Create new slowness vector for wave scattered from fractures.
        :param slowness: Incoming slowness vector.
        :param phi_hat: unit vector orthogonal to fracture planes.
        :param fracture_spacing: Distance between fracture planes in m.
        :param frequency: Frequency of seismic wave in Hz.
        :return: Modified slowness vector.
        """
        ps = slowness[:Index.Z]
        return ps - np.copysign(1, ps @ phi_hat[:Index.Z]) / (fracture_spacing * frequency) * phi_hat

