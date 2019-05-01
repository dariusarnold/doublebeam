from pathlib import Path

import numpy as np

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
        key = "top_p_velocity"
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
        return y_left + slope * (x_position - x_left)

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
        depths = layers["top_depth"]
        depths = np.append(depths, layers["bot_depth"][-1])
        self.interface_depths = depths

    @classmethod
    def from_file(cls, filepath: Path):
        # Load model from file. Formatting and content of file is described
        # in README..md
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

    @classmethod
    def from_string(cls, model: str):
        """
        Create model from string description. The same rules as for model files
        apply for creation from string
        """
        # np.loadtxt also takes generators, so create one
        try:
            raw_data = np.loadtxt((line for line in model.split("\n")), dtype=LinearVelocityLayer, delimiter=",")
        except IndexError:
            try:
                raw_data = np.loadtxt((line for line in model.split("\n")), dtype=ConstantVelocityLayer, delimiter=",")
            except ValueError:
                msg =f"Error parsing velocity model string {model}"
                raise ValueError(msg)
        return cls(raw_data)

    def _layer_number(self, depth_km: float) -> int:
        """
        Find the layer within the model that contains the depth and return
        its index in self.layers
        """
        # Get mask: True for the layer in which the depth is
        layer = np.logical_and(self.layers['top_depth'] <= depth_km,
                               depth_km < self.layers['bot_depth'])
        # get indices of the True value. meaning index of the wanted layer
        layer = np.where(layer)[0]
        if len(layer):
            # return the layer which top depth is directly above depth_km
            return layer[0]
        else:
            # depth must have been outside of VelocityModel since all indices were false
            raise LookupError(f"No layer found in model at depth {depth_km}")

    def eval_at(self, depth_km: float, prop: str) -> float:
        layer = self.layers[self._layer_number(depth_km)]
        if layer.dtype == ConstantVelocityLayer:
            return evaluate_at(layer, prop)
        elif layer.dtype == LinearVelocityLayer:
            return evaluate_at_linear(layer, depth_km, prop)

    def interface_crossed(self, z1: float, z2: float) -> bool:
        """Return true if at least one interface is between the depths z1 and z2."""
        z_upper = min(z1, z2)  # higher point
        z_lower = max(z1, z2)  # lower point
        # get mask: true for all interfaces which lay between the two points
        interfaces_between = np.logical_and(self.interface_depths < z_lower, self.interface_depths > z_upper )
        return np.any(interfaces_between)

    def __len__(self):
        return len(self.layers)
