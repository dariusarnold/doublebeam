import numpy as np
from pathlib import Path


# Layer which has a constant velocity over its depth
ConstantVelocityLayer = np.dtype([
    ('top_depth', np.float_),
    ('bot_depth', np.float_),
    ('top_p_velocity', np.float_),
    ('top_s_velocity', np.float_),
    ('top_density', np.float_),
])

# Layer which has a linear velocity gradient over its depth
LinearVelocityLayer = np.dtype([
    ('top_depth', np.float_),
    ('bot_depth', np.float_),
    ('top_p_velocity', np.float_),
    ('bot_p_velocity', np.float_),
    ('top_s_velocity', np.float_),
    ('bot_s_velocity', np.float_),
    ('top_density', np.float_),
    ('bot_density', np.float_),
])


def evaluate_at(layer: ConstantVelocityLayer, prop: str) -> float:
    """
    Evaluate material properties of a velocity layer
    :param layer: The layer to evaluate
    :param prop: Which property to evaluate.
    "p" = P wave velocity (km/s)
    "s" = S wave velocity (km/s)
    "r" or "d" = density (gm/cm^3)
    :return: The requested layer property
    """

    if prop == "p":
        key = "top_density"
    elif prop == "s":
        key = "top_s_velocity"
    elif prop in "rd":
        key = "top_density"
    return layer[key]


def evaluate_at_linear(layer: LinearVelocityLayer, depth: float, prop: str) -> float:
    """
    Evaluate material properties at a given depth in a layer with linearly
    varying properties
    :param layer: The layer to evaluate
    :param depth: Depth in km at which to evaluate the property
    :param prop: Which property to evaluate
    "p" = P wave velocity (km/s)
    "s" = S wave velocity (km/s)
    "r" or "d" = density (gm/cm^3)
    :return: The requested layer property
    """

    def eval_linear_gradient(x_left: float, x_right: float, y_left: float,
                             y_right: float, x_position: float) -> float:
        """
        Used for interpolation between the layers property at top and bottom.
        The gradient is given by two points, which have an x value and a
        corresponding y value. Linear interpolation between the two points at the
        position given by value is used for the result
        :return: Value of linear gradient at x_position
        """
        slope = (y_right - y_left) / (x_right - x_left)
        return y_left + slope * x_position

    if prop == "p":
        top_key, bottom_key = "top_p_velocity", "bot_p_velocity"
    elif prop == "s":
        top_key, bottom_key = "top_s_velocity", "bot_s_velocity"
    elif prop in "rd":
        top_key, bottom_key = "top_density", "bot_density"

    top_value, bottom_value = layer[top_key], layer[bottom_key]
    return eval_linear_gradient(layer["top_depth"], layer["bot_depth"],
                                top_value, bottom_value, depth)



class VelocityModel1D:
    """Simple vertical velocity model storing P and S wave velocity as well
    as density"""

    def __init__(self, layers):
        self.layers = layers

    @classmethod
    def from_file(cls, filepath: Path):
        try:
            # this fails if there are less than exactly 8 values, so next try
            # to convert to ConstantVelocityLayer
            raw_data = np.loadtxt(str(filepath), dtype=LinearVelocityLayer, delimiter=",")
        except IndexError:
            try:
                raw_data = np.loadtxt(str(filepath), dtype=ConstantVelocityLayer, delimiter=",")
            except ValueError:
                msg = f"Error parsing velocity model file {str(filepath)}"
                raise ValueError(msg)
        return cls(raw_data)

    def layer_number(self, depth_km: float) -> int:
        """
        Find the layer within the model that contains the depth and return
        its index in self.layers
        """
        # Get mask: True for all layers that are above the given depth
        layer = self.layers["top_depth"] <= depth_km
        # get indices of all True values, meaning indices of all layers above
        # the given depth
        layer = np.where(layer)[0]
        if len(layer):
            # return the layer which top depth is directly above depth_km
            return layer[-1]
        else:
            raise LookupError(f"No layer found in model at depth {depth_km}")

    def eval_at(self, depth_km: float, prop: str) -> float:
        layer = self.layers[self.layer_number(depth_km)]
        if layer.dtype == ConstantVelocityLayer:
            return evaluate_at(layer, prop)
        elif layer.dtype == LinearVelocityLayer:
            return evaluate_at_linear(layer, depth_km, prop)

    def __len__(self):
        return len(self.layers)
