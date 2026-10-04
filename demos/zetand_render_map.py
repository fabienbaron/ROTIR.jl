#!/usr/bin/env python3
"""Render the extracted ζ And map the same way the original figure draws it, for comparison.

    python3 demos/zetand_render_map.py

Writes `demos/zetAnd/zetAnd2019A_map_extracted.png` and a side-by-side
`zetAnd2019A_map_compare.png` against the raster pulled out of the PDF.

Matched to the original on purpose: Mollweide, the colormap read back OUT of the figure's own
colourbar (it is `gist_heat` — matplotlib's version of IDL's Red Temperature — to 0.0125 rms,
but using the bar itself makes the two panels identical by construction), the same 0.8025-0.9977
colourbar limits, the same unlabelled 30-degree graticule, and the latitude convention of the
FIGURE — `zetand_extract_map.py` flips latitude on the way in (that is the parity ζ And's
visibilities select), so it is flipped back on the way out. If the pipeline is right, the two
panels show the same star.

HOW CLOSE THE TWO PANELS ACTUALLY ARE, measured per pixel over the disc rather than left to
the eye: bias +0.0008, rms 0.0226 — 11.6 % of the 0.195 colourbar range — correlation 0.85.
Splitting on whether a pixel sits at a facet boundary (29 % of the disc do) gives rms 16.2 %
there against 8.9 % in facet interiors, where the correlation is 0.91. So the residual is
resampling: her tessellation is not ours, facet centres do not coincide, and every boundary
pixel in the original is a blend of two colours. Resampling at finer HEALPix levels improves it
smoothly (12.5, 11.6, 9.8, 8.4 % of range at nside_exp 2, 3, 4, 5) without converging, which is
what a grid mismatch looks like and not what a bad colour inversion looks like. The test that
the extraction is RIGHT is not this one — it is that the map fits the visibilities, in
`zetand_published_map.jl`.

The map is drawn by INVERSE projection — every output pixel is carried back to (lon, lat) and
takes the value of the nearest tessel centre — rather than by filling projected quadrilaterals.
Same result at this resolution, and it needs no polygon clipping at the 180-degree meridian,
where a facet's corners otherwise wrap and smear across the whole plot.
"""
import os
import numpy as np
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
ZET  = os.path.join(HERE, "zetAnd")
VMIN, VMAX = 0.8025, 0.9977
NX, NY = 900, 450

def figure_cmap():
    """The colormap sampled from the figure's own colourbar, so the panels cannot differ."""
    raster = os.path.join(ZET, "zetAnd2019A_orig.png")
    bar = np.asarray(Image.open(raster).convert("RGB")).astype(float)[470, 220:618]/255.0
    return colors.ListedColormap(bar)

def inverse_mollweide(xn, yn):
    """Unit box -> (lon, lat) in degrees; NaN outside the ellipse."""
    with np.errstate(invalid="ignore"):
        th = np.arcsin(np.clip(yn, -1, 1))
        lat = np.arcsin(np.clip((2*th + np.sin(2*th))/np.pi, -1, 1))
        lon = np.pi*xn/np.maximum(np.cos(th), 1e-12)
    bad = (xn**2 + yn**2 > 1.0) | (np.abs(lon) > np.pi)
    lon = np.where(bad, np.nan, lon); lat = np.where(bad, np.nan, lat)
    return np.degrees(lon), np.degrees(lat)

def render(values, lonlat, ax, title):
    lon_t, lat_t = lonlat[:, 0], lonlat[:, 1]
    # tessel centres as unit vectors, with latitude FLIPPED back to the figure's convention
    cl, sl = np.cos(np.radians(lon_t)), np.sin(np.radians(lon_t))
    cb, sb = np.cos(np.radians(-lat_t)), np.sin(np.radians(-lat_t))
    tree = cKDTree(np.column_stack([cb*cl, cb*sl, sb]))
    xn, yn = np.meshgrid(np.linspace(-1, 1, NX), np.linspace(1, -1, NY))
    lon, lat = inverse_mollweide(xn, yn)
    good = np.isfinite(lon)
    cl2, sl2 = np.cos(np.radians(lon[good])), np.sin(np.radians(lon[good]))
    cb2, sb2 = np.cos(np.radians(lat[good])), np.sin(np.radians(lat[good]))
    _, idx = tree.query(np.column_stack([cb2*cl2, cb2*sl2, sb2]))
    img = np.full(lon.shape, np.nan)
    img[good] = values[idx]
    # x spans 2 units against y's 1: Mollweide's box is 2:1, and an equal-aspect [-1,1]^2
    # extent draws the ellipse as a circle, which looks like a different projection.
    im = ax.imshow(img, origin="upper", extent=[-2, 2, -1, 1], cmap=figure_cmap(),
                   vmin=VMIN, vmax=VMAX, interpolation="nearest")
    # the same furniture as the original: outline, solid equator and central meridian,
    # dotted graticule every 30 degrees, no labels
    t = np.linspace(0, 2*np.pi, 512)
    ax.plot(2*np.cos(t), np.sin(t), "k-", lw=1.0)
    for lat0 in range(-60, 90, 30):
        ll = np.linspace(-180, 180, 361)
        th = np.arcsin(np.clip(np.sin(np.radians(lat0)), -1, 1))
        for _ in range(40):
            th -= (2*th + np.sin(2*th) - np.pi*np.sin(np.radians(lat0)))/(2 + 2*np.cos(2*th))
        ax.plot(2*np.radians(ll)*np.cos(th)/np.pi, np.full_like(ll, np.sin(th)),
                "k:", lw=0.6, dashes=(1, 3))
    for lon0 in range(-150, 180, 30):
        la = np.linspace(-89.9, 89.9, 361)
        th = np.radians(la).copy()
        for _ in range(40):
            th -= (2*th + np.sin(2*th) - np.pi*np.sin(np.radians(la)))/(2 + 2*np.cos(2*th))
        style = dict(ls="-", lw=1.0) if lon0 == 0 else dict(ls=":", lw=0.6, dashes=(1, 3))
        ax.plot(2*np.radians(lon0)*np.cos(th)/np.pi, np.sin(th), color="k", **style)
    ax.plot([-2, 2], [0, 0], "k-", lw=1.0)
    ax.set_xlim(-2.1, 2.1); ax.set_ylim(-1.1, 1.1); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, fontsize=11)
    return im

def main():
    values = np.loadtxt(os.path.join(ZET, "zetAnd2019A_map.txt"))
    lonlat = np.loadtxt(os.path.join(ZET, "tessel_lonlat_3.txt"))
    print(f"{len(values)} tessels, {np.nanmin(values):.4f}..{np.nanmax(values):.4f}")

    fig, ax = plt.subplots(figsize=(8, 5))
    im = render(values, lonlat, ax, "zet And 2019A — extracted, HEALPix nside_exp 3")
    cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.05, pad=0.02,
                      ticks=[VMIN, VMAX])
    cb.ax.set_xticklabels([f"{VMIN:.4f}", f"{VMAX:.4f}"])
    out = os.path.join(ZET, "zetAnd2019A_map_extracted.png")
    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white"); plt.close(fig)
    print("wrote", out)

    orig = os.path.join(ZET, "zetAnd2019A_orig.png")
    if os.path.exists(orig):
        fig, axes = plt.subplots(1, 2, figsize=(15, 5))
        axes[0].imshow(np.asarray(Image.open(orig))); axes[0].axis("off")
        axes[0].set_title("original figure (from the PDF)", fontsize=11)
        im = render(values, lonlat, axes[1], "extracted and re-rendered")
        out2 = os.path.join(ZET, "zetAnd2019A_map_compare.png")
        fig.savefig(out2, dpi=150, bbox_inches="tight", facecolor="white"); plt.close(fig)
        print("wrote", out2)

if __name__ == "__main__":
    main()
