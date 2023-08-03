#!/usr/bin/env python3
"""
Basic visualization tools for SeisFlows to visualize waveforms, models, etc.
"""
import os
import numpy as np
import matplotlib.pyplot as plt

from seisflows.tools import unix


def plot_waveforms(tr_obs, tr_syn, tr_adj=None, fid_out=None, **kwargs):
    """
    Very simple plotting routine to show waveforms and adjoint sources
    manipulated by the Default preprocessing module. Plots are simple and are
    provided in a default style that can be adjusted via keyword arguments.

    :type tr_obs: obspy.core.stream.Stream
    :param tr_obs: observed seismogram, data
    :type tr_syn: obspy.core.stream.Stream
    :param tr_syn: synthetic seismogram
    :type tr_adj: obspy.core.stream.Stream
    :param tr_adj: optional adjoint source. if not given, none will be plotted
    :type fid_out: str
    :param fid_out: name and path to save output file. If none given, output
        file will be saved to current working directory and named based on the
        trace ID of the obs data
    """
    dpi = kwargs.get("dpi", 100)
    figsize = kwargs.get("figsize", (800 / dpi, 300 / dpi))
    lw = kwargs.get("linewidth", 1)
    obs_color = kwargs.get("obs_color", "k")
    syn_color = kwargs.get("syn_color", "r")
    adj_color = kwargs.get("adj_color", "g")

    f, ax = plt.subplots(figsize=figsize, dpi=dpi)

    lines = []  # for legend
    lines += ax.plot(tr_obs.times(), tr_obs.data, c=obs_color, lw=lw, 
                     label="obs", zorder=6)
    lines += ax.plot(tr_syn.times(), tr_syn.data, c=syn_color, lw=lw, 
                     label="syn", zorder=6)

    if tr_adj is not None:
        twax = ax.twinx()
        lines += twax.plot(tr_adj.times(), tr_adj.data, c=adj_color, lw=lw,
                           label="adj", ls="--", alpha=0.75, zorder=5)
        twax.set_ylabel("Adj. Amplitude")

    plt.title(tr_syn.id)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Amplitude")

    labels = [l.get_label() for l in lines]
    ax.legend(lines, labels, loc="upper right")

    if not fid_out:
        fid_out = f"./{tr_syn.id.replace('.', '_')}.png"

    # Overwrite existing figures
    if os.path.exists(fid_out):
        unix.rm(fid_out)

    plt.tight_layout()
    plt.savefig(fid_out)
    plt.close()


def plot_2d_contour(x, z, data, cmap="viridis", zero_midpoint=False):
    """
    Plots values of a SPECEFM2D model/gradient on an unstructured grid

    :type x: np.array
    :param x: x values of GLL mesh
    :type z: np.array
    :param z: z values of GLL mesh
    :type data: np.array
    :param data: D
    :type cmap: str
    :param cmap: matplotlib colormap to be applied to the contour plot. Defaults
        to 'viridis'
    :type zero_midpoint: bool
    :param zero_midpoint: set 0 as the midpoint for the colorbar. Useful for
        diverging colorscales (e.g., for gradients), where the neutral color
        (e.g., white) is set at value=0
    """
    # Figure out aspect ratio of the figure
    r = (max(x) - min(x))/(max(z) - min(z))
    rx = r/np.sqrt(1 + r**2)
    ry = 1/np.sqrt(1 + r**2)

    # Assign zero as the midpoint for things like gradients
    if zero_midpoint:
        abs_max_val = max(abs(data))
        vmin = -1 * abs_max_val
        vmax = abs_max_val
    else:
        vmin, vmax = None, None

    f = plt.figure(figsize=(10 * rx, 10 * ry))
    p = plt.tricontourf(x, z, data, levels=125, cmap=cmap, vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(p, shrink=0.8, pad=0.025) # , format="%.2f")
    plt.axis("image")

    return f, p, cbar


def plot_2d_image(x, z, data, cmap="viridis", zero_midpoint=False,
                  resX=1000, resZ=1000):
    """
    Plots values of a SPECEFM2D model/gradient by interpolating onto a regular grid

    :type x: np.array
    :param x: x values of GLL mesh
    :type z: np.array
    :param z: z values of GLL mesh
    :type data: np.array
    :param data: D
    :type cmap: str
    :param cmap: matplotlib colormap to be applied to the contour plot. Defaults
        to 'viridis'
    :type zero_midpoint: bool
    :param zero_midpoint: set 0 as the midpoint for the colorbar. Useful for
        diverging colorscales (e.g., for gradients), where the neutral color
        (e.g., white) is set at value=0
    :type resX: int
    :param resX: number of points for the interpolation in x- direction (default=1000)
    :type resZ: int
    :param resZ: number of points for the interpolation in z- direction (default=1000)
    """
    # Figure out aspect ratio of the figure
    r = (max(x) - min(x))/(max(z) - min(z))
    rx = r/np.sqrt(1 + r**2)
    ry = 1/np.sqrt(1 + r**2)

    # Assign zero as the midpoint for things like gradients
    if zero_midpoint:
        abs_max_val = max(abs(data))
        vmin = -1 * abs_max_val
        vmax = abs_max_val
    else:
        vmin, vmax = None, None

    f = plt.figure(figsize=(10 * rx, 10 * ry))
    from scipy.interpolate import griddata

    # trick interpolation using the maximum values of z in case of concave topography.
    # nan values helps interpolation act expectedly.
    # Can be tested using the default specfem2D model: simple_topography_and_also_a_simple_fluid_layer
    x = np.append(x, [min(x), max(x)])
    z = np.append(z, [max(z), max(z)])
    data = np.append(data, [np.nan, np.nan])

    xi = np.linspace(min(x), max(x), resX)
    zi = np.linspace(min(z), max(z), resZ)
    X, Z = np.meshgrid(xi, zi)
    V = griddata((x, z), data, (X, Z), method='linear')
    im = plt.imshow(V, vmax=vmax, vmin=vmin,
                    extent=[x.min(), x.max(), z.min(), z.max()],
                    cmap=cmap,
                    origin='lower')

    cbar = plt.colorbar(im, shrink=0.8, pad=0.025)
    plt.axis("image")

    return f, im, cbar

