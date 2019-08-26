from pathlib import Path
from typing import Tuple, List, Union

import numpy as np

from doublebeam.utils import Index

# Layer which has a linear velocity gradient over its depth
# v(z) = intercept + z * gradient
LinearVelocityLayer = np.dtype([
    ('top_depth', np.float64),
    ('bot_depth', np.float64),
    ('intercept', np.float64),
    ('gradient', np.float64)
])


def evaluate_at_linear(layer: LinearVelocityLayer, depth: float) -> float:
    """
    Evaluate velocity at a given depth in a layer with linearly varying velocity
    :param layer: The layer to evaluate
    :param depth: Depth in m at which to evaluate the property
    """
    return layer["intercept"] + layer["gradient"] * depth


class VelocityModel3D:
    """Velocity model storing wave velocity for horizontal layers."""

    def __init__(self, layers: Union[List[Tuple[float, float, float, float]], np.ndarray],
                 x_width: float = 0., y_width: float = 0.):
        """
        :param layers: A layer is specified by the depth value of its top
        coordinate, the depth of its bottom, the velocity intercept and the
        velocity gradient. The values are given as a list of tuples, where every
        tuple represents one layer. The order of values in the tuple is as
        described above.
        # TODO Will it be necessary to implement checking these in eval functions?
        :param x_width: Width of model in meters along x axis.
        :param y_width: Width of model in meters along y axis.
        """
        self.x_width = x_width
        self.y_width = y_width
        self.layers = np.asarray(layers, dtype=LinearVelocityLayer)
        depths = self.layers["top_depth"]
        depths = np.append(depths, self.layers["bot_depth"][-1])
        self.interface_depths = depths
        self.layer_heights = self.layers["bot_depth"] - self.layers["top_depth"]
        self.gradients = self.layers["gradient"]
        self.intercepts = self.layers["intercept"]
        self.velocities_top = self.layers["intercept"] + self.layers["gradient"] * self.layers["top_depth"]
        self.velocities_bot = self.layers["intercept"] + self.layers["gradient"] * self.layers["bot_depth"]

    def __getitem__(self, index: int) -> np.ndarray:
        """
        Implement indexing to get layers
        :param index: Number of layer, 0 for top most, increasing downwards.
        """
        return self.layers[index]

    @classmethod
    def convert_to_gradient_intercept(cls, layers: Union[List[Tuple[float, float, float, float]],
                                                         np.ndarray]) -> np.ndarray:
        """
        Convert layers specified by their top and bottom depth and their top and
        bottom velocity to intercept and gradient for the velocity
        :param layers: List of layers, where every layer is specified by four
        values in the listed order:
        Depth top (m), depth bottom (m), velocity top (m/s), velocity bottom (m/s)
        """
        layers = np.asarray(layers).T
        if layers.shape[0] != 4:
            raise ValueError(f"Wrong input shape, expected (N, 4), "
                             f"got {layers.T.shape}.")
        heights = layers[1] - layers[0]
        gradients = (layers[3] - layers[2]) / heights
        intercepts = layers[2] - gradients * layers[0]
        return np.array([(dt, db, i, g) for dt, db, i, g in
                         zip(layers[0], layers[1], intercepts, gradients)],
                        dtype=LinearVelocityLayer)

    @classmethod
    def from_file(cls, filepath: Union[Path, str], *args, **kwargs) -> "VelocityModel3D":
        """
        Load model from file. Formatting and content of file is described in
        README.md.
        Use args or kwargs to pass additional parameters, such as the width, to
        the model constructor.
        :param filepath: Path to to file containing model data
        """
        # TODO update README to reflect changes: Only LinearVelocityLayer is kept
        try:
            raw_data = np.loadtxt(str(filepath), delimiter=",")
        except IndexError:
            msg = f"Error parsing velocity model file {str(filepath)}"
            raise ValueError(msg)
        raw_data = cls.convert_to_gradient_intercept(raw_data)
        return cls(raw_data, *args, **kwargs)

    @classmethod
    def from_string(cls, model: str, *args, **kwargs) -> "VelocityModel3D":
        """
        Create model from string description. The same rules as for model files
        apply for creation from string.
        Use args or kwargs to pass additional parameters, such as the width, to
        the model constructor.
        :param model: Velocity model string.
        """
        # np.loadtxt also takes generators, so create one
        try:
            raw_data = np.loadtxt((line for line in model.split("\n")), delimiter=",")
        except IndexError:
            msg = f"Error parsing velocity model string {model}"
            raise ValueError(msg)
        raw_data = cls.convert_to_gradient_intercept(raw_data)
        return cls(raw_data, *args, **kwargs)

    def layer_index(self, points: Union[float, np.ndarray]) -> Union[int, np.ndarray]:
        """
        Find the layer within the model that contains the point and return
        its index in self.layers. While the top of a layer belongs to the layer,
        the bottom depth belongs to the layer below, except for the bottom most
        layer which includes its bottom depth.
        :param points: Single float value specifying depth or array of shape
        (N, 3) where every entry is a x, y, z coordinate triple.
        """
        try:
            depths = points.T[Index.Z]
        except (AttributeError, IndexError):
            depths = np.array(points)
        top, bottom = self.vertical_boundaries()
        if np.any(np.logical_or(depths < top, depths > bottom)):
            raise LookupError(f"Depth {depths} contains values out of model "
                              f"range: {top}, {bottom}")
        # wrap in asarray since searchsorted returns scalar for single point
        # which then doesn't support assignment
        indices = np.asarray(np.searchsorted(self.interface_depths, depths,
                                             side="right") - 1)
        indices[indices >= len(self)] -= 1
        return indices

    def eval_at(self, points: Union[float, np.ndarray]) -> float:
        """
        Return velocity at a point in the model
        :param points: Array of points of shape (N, 3) where every entry is a
        tuple of three x, y, z coordinates
        :return: Velocity in m/s
        """
        # TODO add exception when an array of the wrong shape is passed
        try:
            z = points.T[Index.Z]
        except (AttributeError, IndexError):
            z = np.array(points)
        top_depth, bottom_depth = self.vertical_boundaries()
        if np.any(np.logical_or(z < top_depth,  z > bottom_depth)):
            raise LookupError(f"Can't evaluate model at negative depth {z}")
        layer = self.layers[self.layer_index(points)]
        return layer["intercept"] + layer["gradient"] * z

    def interface_crossed(self, z1: float, z2: float) -> bool:
        """
        Return true if at least one interface is between the depths z1 and z2.
        The depths do not have to be ordered.
        :param z1: depth in m
        :param z2: depth in m
        """
        z_upper = min(z1, z2)  # higher point
        z_lower = max(z1, z2)  # lower point
        # get mask: true for all interfaces which lay between the two points
        interfaces_between = np.logical_and(self.interface_depths < z_lower,
                                            self.interface_depths > z_upper)
        return np.any(interfaces_between)

    def interface_velocities(self, z: float) -> Tuple[float, float]:
        """
        Return velocities above and below the closest interface.
        0 is returned for the velocity outside of the model, e.g. when the
        interface between the top layer and the one above is requested.
        :param z: Depth in m
        :return: Tuple of: (velocity above interface, velocity below interface)
        where both velocities are in m/s
        """
        index = self.layer_index(z)
        midpoint_height = self.layers[index]["top_depth"] + self.layer_heights[index] / 2
        if index == 0 and z < midpoint_height:
            # above half depth of top layer
            return 0., self.velocities_top[index]
        if index == len(self) - 1 and z > midpoint_height:
            # below half depth of bottom layer
            return self.velocities_bot[index], 0.
        if z > midpoint_height:
            # eval interface with layer below
            return (self.velocities_bot[index],
                    self.velocities_top[index+1])
        else:
            # eval interface with layer above
            return (self.velocities_bot[index-1],
                    self.velocities_top[index])

    def num_of_interfaces_between(self, point1: np.ndarray, point2: np.ndarray) -> int:
        """
        Return the number of interfaces that lie between two points.
        """
        l1 = self.layer_index(point1[Index.Z])
        l2 = self.layer_index(point2[Index.Z])
        return abs(l1 - l2)

    def vertical_boundaries(self) -> Tuple[float, float]:
        """
        Return tuple of upper and lower boundary of the model.
        """
        return self.interface_depths[0], self.interface_depths[-1]

    def __len__(self) -> int:
        """
        Return number of layers in the model.
        """
        return len(self.layers)
