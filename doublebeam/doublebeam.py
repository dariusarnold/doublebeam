
import numpy as np

from doublebeam.models import VelocityModel3D
from doublebeam.raytracing.twopoint import TwoPointRayTracing
from doublebeam.utils import generate_grid_coordinates, Index, unit_vector


def scattered_slowness(slowness: np.ndarray, phi_hat: np.ndarray,
                       fracture_spacing: float, frequency: float) -> np.ndarray:
    """
    Create new slowness vector for wave scattered from fractures.
    :param slowness: Incoming slowness vector.
    :param phi_hat: unit vector orthogonal to fracture planes.
    :param fracture_spacing: Distance between fracture planes in m.
    :param frequency: Frequency of seismic wave in Hz.
    :return: Modified slowness vector.
    """
    slowness = np.array(slowness, dtype=np.float64)
    ps = slowness[:Index.Z]
    ps -= np.copysign(1, ps @ phi_hat[:Index.Z]) / (fracture_spacing * frequency) * phi_hat[:Index.Z]
    return slowness


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
        if model.x_width == 0 or model.y_width == 0:
            raise ValueError("Model width can't be 0!")
        self.model = model
        self.fracture_info = fracture_parameters
        target_num_x, target_num_y = 10, 10
        self.target_locations = generate_grid_coordinates(self.fracture_info.depth, (0, self.model.x_width),
                                                          (0, self.model.y_width), target_num_x, target_num_y)
        source_num_x, source_num_y = 20, 20
        self.source_beam_centers = generate_grid_coordinates(self.model.vertical_boundaries()[0],
                                                             (0, self.model.x_width), (0, self.model.y_width),
                                                             source_num_x, source_num_y)
        self.window_length = window_length
        self.twopoint = TwoPointRayTracing(model)



