#!/usr/bin/env python3
"""
Basic visualization tools for SeisFlows to visualize waveforms, models, etc.
Includes utility functions for manipulating image files such as .png and .pdfs.
Makes use use of the Python Pillow package if working with .png files, and
the PyPDF2 package if manipulating pdf files. Internal imports for all functions
to remove Pyatoa-wide dependencies on these packages for short functions.
"""
import os
import numpy as np
import matplotlib.pyplot as plt

from PIL import Image
from pypdf import PdfMerger
from seisflows.tools import unix


def merge_pdfs(fids, fid_out, remove_fids=False):
    """
    Merge a list of pdfs into a single output pdf using the PyPDF2 package.
    Any desired order to the pdfs should be set in the list of input fids.

    :type fids: list
    :param fids: list of paths to .pdf files
    :type fid_out: str
    :param fid_out: path and name of the resulting output .pdf file
    :type remove_fids: bool
    :param remove_fids: remove the original input PDF files if the output PDF
        was successfully created. Defaults to False, original files are kept.
    """
    merger = PdfMerger()
    for fid in fids:
        merger.append(fid)

    merger.write(fid_out)
    merger.close()

    if remove_fids:
        assert(os.path.exists(fid_out)), (
            f"PDF merging failed, cannot remove original input files"
        )
        for fid in fids:
            os.remove(fid)


def imgs_to_pdf(fids, fid_out, remove_fids=False):
    """
    Combine a list of .png files into a single PDF document

    :type fids: list
    :param fids: list of file ids with full pathnames to be combined
    :type fid_out: str
    :param fid_out: the name of the file to be saved with full pathname
    :type remove_fids: bool
    :param remove_fids: remove the original input PNG files if the output PDF
        was successfully created. Defaults to False, original files are kept.
    """
    images = []
    for fid in fids:
        # PNGs need to be converted to RGB to get alpha to play nice
        images.append(Image.open(fid).convert("RGB"))

    image_main = images[0]
    images = images[1:]

    image_main.save(fp=fid_out, format="PDF", resolution=100., save_all=True,
                    append_images=images)

    if remove_fids:
        assert(os.path.exists(fid_out)), (
            f"PNG merging failed, cannot remove original input files"
        )
        for fid in fids:
            os.remove(fid)


def tile_imgs(fids, fid_out):
    """
    Combine a list of images into a single, horizontally tiled image.

    :type fids: list
    :param fids: list of file ids with full pathnames to be tiled
    :type fid_out: str
    :param fid_out: the name of the file to be saved with full pathname
    """
    # .png files require conversion to properly get the alpha layer
    images = []
    for fid in fids:
        images.append(Image.open(fid).convert("RGBA"))

    widths, heights = zip(*(i.size for i in images))
    total_width = sum(widths)
    max_height = max(heights)

    # Create the new image that will be returned
    im_out = Image.new(mode="RGBA", size=(total_width, max_height))
    x_offset = 0
    for im in images:
        im_out.paste(im=im, box=(x_offset, 0))
        x_offset += im.size[0]

    im_out.save(fid_out)


def tif_to_array(fid):
    """
    Convert GeoTiff images (e.g., ETOPO1 topography) to a numpy array for
    conversion and processing

    :type fid: str
    :param fid: .tif(f) file
    :rtype: np.array
    :return: array of data contained within the tiff file
    """
    try:
        im = Image.open(fid)
    except Image.DecompressionBombError as e:
        # If the image is too large, it will throw an error, we'll just need
        # to adjust the acceptable size of the image. This may be bad if you
        # don't trust the image! But let's be reckless...
        error_str = str(e)
        # Assuming the error message looks like:
        # 'Image size (N pixels) exceeds...'  Trying to get value N
        pixels = int(error_str.split()[2][1:])
        print(f"setting Image max pixel count to {pixels}")
        Image.MAX_IMAGE_PIXELS = pixels
        im = Image.open(fid)

    imarray = np.array(im)
    return imarray


def plot_waveforms(tr_obs, tr_syn, tr_adj=None, fid_out=None, title=None,
                   **kwargs):
    """
    Very simple plotting routine to show waveforms and adjoint sources
    manipulated by the Default preprocessing module. Plots are simple and are
    provided in a default style that can be adjusted via keyword arguments.

    Cuts the x-axis (time) to the length of the synthetic seismogram as data
    may be longer which tends to obscure signal of interest.

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

    # Plot observed and synthetic seismograms
    lines = []  # for legend
    lines += ax.plot(tr_obs.times(), tr_obs.data, c=obs_color, lw=lw, 
                     label="obs", zorder=6)
    lines += ax.plot(tr_syn.times(), tr_syn.data, c=syn_color, lw=lw, 
                     label="syn", zorder=6)

    # Plot adjoint source if provided
    if tr_adj is not None:
        twax = ax.twinx()
        lines += twax.plot(tr_adj.times(), tr_adj.data, c=adj_color, lw=lw,
                           label="adj", ls="--", alpha=0.75, zorder=5)
        twax.set_ylabel("Adj. Amplitude")
        twax_ylim = max(abs(tr_adj.data.min()), abs(tr_adj.data.max()))
        twax.set_ylim([-1 * twax_ylim, twax_ylim])

    # Set the x-axis limits based on syn, which is what we are interested in
    ax.set_xlim([tr_syn.times()[0], tr_syn.times()[-1]])
    ax_ylim = max(abs(tr_syn.data.min()), abs(tr_syn.data.max()))
    ax.set_ylim([-1 * ax_ylim, ax_ylim])

    # Plot attributes and labels
    if title is None:
        title = tr_syn.id
    plt.title(title)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Amplitude")

    labels = [l.get_label() for l in lines]
    ax.legend(lines, labels, loc="upper right")

    # Determine where to save the figure
    if not fid_out:
        fid_out = f"./{tr_syn.id.replace('.', '_')}.png"

    # Overwrite existing figures
    if os.path.exists(fid_out):
        unix.rm(fid_out)

    plt.tight_layout()
    plt.savefig(fid_out)
    plt.close()

def plot_optim_stats(fid="output_optim.txt", path_out="./"):
    """
    Line plot of optimization stats, which are written about by 
    Optimize.write_stats(). Intrinsically tied to the format of input file.

    :type fid: str
    :param fid: path to the optimization stats file to plot
    :type path_out: str
    :param path_out: full path (no filename) to save figures. filenames will
        be determined by the header values
    """
    header = open(fid).readlines()[0].strip().split(",")
    stats = np.loadtxt(fid, delimiter=",", skiprows=1).T

    for h, s in zip(header, stats):
        plt.plot(s, "ko-")
        plt.xlabel("Iteration")
        plt.ylabel(h)
        plt.savefig(os.path.join(path_out, f"{h}.png"))
        plt.close("all")

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
    Plots values of a SPECEFM2D model/gradient by interpolating onto a regular 
    grid

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
    :param resX: number of points for the interpolation in x- direction 
        (default=1000)
    :type resZ: int
    :param resZ: number of points for the interpolation in z- direction 
    (default=1000)
    """
    from scipy.interpolate import griddata

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

    # trick interpolation using the maximum values of z in case of concave 
    # topography. nan values helps interpolation act expectedly.
    # Can be tested using the default specfem2D model: 
    # simple_topography_and_also_a_simple_fluid_layer
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

