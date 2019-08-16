from typing import Optional
import numpy as np
from doublebeam.models import VelocityModel3D
from doublebeam.raytracing.twopoint import TwoPointRayTracing
from doublebeam.utils import unit_vector


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
