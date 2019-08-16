import matplotlib.pyplot as plt
import numpy as np

from doublebeam.models import VelocityModel3D
from doublebeam.raytracing.ray import Ray3D
from doublebeam.utils import Index
from mpl_toolkits.mplot3d import Axes3D


def plot_seismogram_gather(seismograms: np.ndarray) -> None:
    """Plot all seismograms from a shot to recreate fig 5b) from Hu2018a"""
    plt.pcolormesh(seismograms.T, cmap="RdGy", antialiased=True)
    cb = plt.colorbar()
    # invert y axis so origin is in top left
    plt.ylim(plt.ylim()[::-1])
    # ugly and overly specific way to limit the plotting to 1.5-3 secs
    # this is valid for 4 secs trace length with dt = 0.004 s
    plt.ylim(ymin=750, ymax=375)
    # label y axis with seconds
    plt.yticks(np.linspace(375, 750, 4), [f"{x:.2f}" for x in np.linspace(1.5, 3, 4)])
    cb.set_label("Amplitude")
    plt.xlabel("Trace #")
    plt.ylabel("Time (s)")
    plt.show()


def plot_ray_in_model_2D(ray: Ray3D, velocity_model: VelocityModel3D) -> None:
    """
    Plot a ray path projected on the x-z plane and draw horizontal lines
    for interfaces.
    """
    fig, ax = plt.subplots()
    ax.invert_yaxis()
    ax.set_ylabel("Depth (m)")
    ax.set_xlabel("Offset along x (m)")
    for depth in velocity_model.interface_depths:
        ax.axhline(depth, color="k")
    for segment in ray.path:
        x, z = segment.T[Index.X], segment.T[Index.Z]
        ax.plot(x, z)
    plt.show()


def plot_ray_in_model_3D(ray: Ray3D, velocity_model: VelocityModel3D) -> None:
    """
    Create and show 3D plot of ray path
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.invert_zaxis()
    for segment in ray.path:
        x, y, z = segment.T
        ax.plot(x, y, z)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_zlabel("z (m)")
    # use extent of axes and use them as limits for the planes marking interface
    # boundaries.
    x_limits = ax.get_xlim()
    y_limits = ax.get_ylim()
    x, y = np.meshgrid(x_limits, y_limits)
    for depth in velocity_model.interface_depths:
        z = np.full_like(x, depth)
        ax.plot_surface(x, y, z, color="lightgrey", alpha=.25)
    plt.show()


def plot_velocity_model_diagram(model: VelocityModel3D) -> None:
    """
    Create plot like fig. 3 from Fang2019 for a velocity model where v = v(z).
    """
    top, bottom = model.vertical_boundaries()
    number_of_evaluation_points = 1000
    points = np.zeros((number_of_evaluation_points, 3))
    points.T[Index.Z] = np.linspace(top, bottom, 1000)
    velocities = model.eval_at(points)
    plt.plot(velocities, points.T[Index.Z])
    interface_depths = model.interface_depths
    plt.hlines(interface_depths, *plt.xlim(), linestyle="--", alpha=0.4)
    plt.xlabel("Velocity (m/s)")
    plt.ylabel("Depth (m)")
    ax = plt.gca()
    ax.invert_yaxis()
    plt.title("Velocity model")
    plt.show()

