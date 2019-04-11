import numpy as np


VelocityLayer = np.dtype([
    ('top_depth', np.float_),
    #('bot_depth', np.float_),
    ('top_p_velocity', np.float_),
    #('bot_p_velocity', np.float_),
    ('top_s_velocity', np.float_),
    #('bot_s_velocity', np.float_),
    ('top_density', np.float_)
    #('bot_density', np.float_),
])


def evaluate_at(layer: VelocityLayer, depth: float, prop: str) -> float:
    """
    Evaluate material properties at a given depth in a velocity layer
    :param layer: The layer to evaluate
    :param depth: At which depth to evaluate the layer
    :param prop: Which property to evaluate.
    "p" = P wave velocity (km/s)
    "s" = S wave velocity (km/s)
    "r" or "d" = density (gm/cm³)
    :return: The requested layer property
    """
    if prop == "p":
        return layer["top_density"]
    elif prop == "s":
        return layer["top_s_velocity"]
    elif prop in "rd":
        return layer["top_density"]


class VelocityModel1D:
    """Simple vertical velocity model storing P and S wave velocity as well
    as density"""

    def __init__(self, layers):
        self.layers = layers

    @classmethod
    def from_file(cls, filepath):
        raw_data = np.loadtxt(filepath, dtype=VelocityLayer, delimiter=",")
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
            # return the layer which top depths is directly above depth_km
            return layer[-1]
        else:
            raise LookupError(f"No layer found at depth {depth_km}")


    def eval_at(self, depth_km: float, prop: str) -> float:
        layer = self.layers[self.layer_number(depth_km)]
        return evaluate_at(layer, depth_km, prop)


    def __len__(self):
        return len(self.layers)
