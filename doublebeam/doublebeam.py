import itertools as it
from typing import Tuple

import numpy as np

from doublebeam import GaussBeam
from doublebeam.models import VelocityModel3D
from doublebeam.raytracing.initial_value import DynamicRayTracer3D
from doublebeam.raytracing.twopoint import TwoPointRayTracing
from doublebeam.utils import grid_coordinates, Index, unit_vector, generate_vector_arc


def scattered_slowness(slowness: np.ndarray, phi_hat: np.ndarray,
                       fracture_spacing: float, frequency: float) -> np.ndarray:
    """
    Create new slowness vector for wave scattered from fractures.
    :param fracture_spacing: Distance between fracture planes in m.
    :param slowness: Incoming slowness vector.
    :param phi_hat: unit vector orthogonal to fracture planes. Can also compute
    scattered direction for multiple given fracture normal vectors in the shape
    (N, 3), where N is the number of normal vectors given.
    :param frequency: Frequency of seismic wave in Hz.
    :return: Modified slowness vector.
    """
    slowness = np.array(slowness, dtype=np.float64)
    slowness_z = slowness[Index.Z]
    phi_hat = np.atleast_2d(phi_hat)
    signs = np.expand_dims(np.copysign(1, slowness[:2] @ phi_hat[:, :2].T), -1)
    slowness = slowness - signs / (fracture_spacing * frequency) * phi_hat
    slowness[:, 2] = slowness_z
    return np.squeeze(slowness)


class FractureParameters:

    def __init__(self, fracture_depth: float, fracture_orientation: np.ndarray,
                 frac_spacing_min: float, frac_spacing_max: float, num_fracture_spacings: int,
                 num_fracture_orientations: int):
        """
        :param fracture_depth: Depth of top of fractures in m.
        :param fracture_orientation: Vector orthogonal to fracture planes.
        :param frac_spacing_min: Minimum distance between fracture planes in m.
        :param frac_spacing_max: Maximum distance between fracture planes in m.
        :param num_fracture_spacings: Number of samples between min and max fracture
        spacing.
        :param num_fracture_orientations: Number of samples for fracture
        orientations to scan in a 180° arc around the given fracture direction.
        """
        self.depth = fracture_depth
        self.orientation = unit_vector(fracture_orientation)
        self.spacings: np.ndarray = np.linspace(frac_spacing_min, frac_spacing_max,
                                                num_fracture_spacings)
        self.num_fracture_orientations = num_fracture_orientations


class GridGeometry:

    def __init__(self, offset_x: Tuple[float, float], offset_y: Tuple[float, float],
                 num_points_x: int, num_points_y: int):
        """
        Class describing geometry of a cartesian grid.
        :param offset_x: Tuple of (grid start coordinate along x axis, grid end
        coordinate along x axis).
        :param offset_y: Tuple of (grid start coordinate along y axis, grid end
        coordinate along y axis).
        :param num_points_x: Number of points along x axis.
        :param num_points_y: Number of points along y axis.
        """
        self.offset_x = offset_x
        self.offset_y = offset_y
        self.num_points_x = num_points_x
        self.num_points_y = num_points_y


class DoubleBeamParameters:

    def __init__(self, source_geometry: GridGeometry, target_geometry: GridGeometry,
                 fracture_info: FractureParameters, window_length: float):
        """
        Aggregate class holding information required for the double beam
        algorithm.
        :param source_geometry: Info about geometry of source beam centers.
        :param target_geometry: Info about geometry of target locations.
        :param fracture_info: Info about fracture parameters.
        :param window_length: Time in seconds for windowing seismic data.
        """
        self.source = source_geometry
        self.target = target_geometry
        self.fractures = fracture_info
        self.window_length = window_length


class DoubleBeam:

    def __init__(self, model: VelocityModel3D, params: DoubleBeamParameters):
        if model.x_width == 0 or model.y_width == 0:
            raise ValueError("Model width can't be 0!")

        self.model = model
        self.params = params
        self.target_locations = grid_coordinates(self.params.fractures.depth, self.params.target.offset_x,
                                                 self.params.target.offset_y, self.params.target.num_points_x,
                                                 self.params.target.num_points_y)
        self.source_beam_centers = grid_coordinates(self.model.vertical_boundaries()[0], self.params.source.offset_x,
                                                    self.params.source.offset_y, self.params.source.num_points_x,
                                                    self.params.source.num_points_y)
        self.twopoint = TwoPointRayTracing(model)
        self.dynamic_rt = DynamicRayTracer3D(model)

    def _direct_ray_code(self, source: np.ndarray, receiver: np.ndarray,
                         model: VelocityModel3D) -> str:
        """
        Generate a ray code that represents the direct ray connecting source and receiver
        :param source: x, y, z coordinates of source
        :param receiver: x, y, z coordinates of receiver
        :param model: Velocity model containing interfaces
        """
        return model.num_of_interfaces_between(source, receiver) * "T"

    def algorithm(self, width: float, frequency: float):
        # divide surface shots into source beam centers
        for target, src_beam_center in it.product(self.target_locations.reshape(-1, 3),
                                                  self.source_beam_centers.reshape(-1, 3)):
            # compute source beam b_s between source beam center and target
            slowness = self.twopoint.trace(src_beam_center, target)
            source_beam = GaussBeam(src_beam_center, slowness, width, frequency)
            self.dynamic_rt.trace_stack(source_beam, self._direct_ray_code(src_beam_center, target, self.model))
            # scan all possible fracture spacings and fracture orientations
            fracture_normals = generate_vector_arc(self.params.fractures.num_fracture_orientations,
                                                   self.params.fractures.orientation)
            # find scattered beam direction and construct receiver beam b_g
            slowness_scattered = scattered_slowness(slowness, fracture_normals, self.params.fractures.spacings,
                                                    frequency)
