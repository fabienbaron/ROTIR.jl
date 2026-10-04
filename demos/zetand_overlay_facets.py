#!/usr/bin/env python3
"""Overplot the HEALPix facet edges, and the points the extractor samples, on ζ And's figure.

    python3 demos/zetand_overlay_facets.py

Writes into `demos/zetAnd/`:

    zetAnd2019A_facets.png       the whole figure with facet outlines (blue) and sample points
                                 (green) drawn on it, 3x upscaled
    zetAnd2019A_facets_tip.png   a zoom on the left tip at the equator, which is where a
                                 projection or registration error shows up first
    zetAnd2019A_tip_compare.png  that same tip, original beside our re-render

This is the check that cannot be done numerically: whether the facet grid we believe is in the
figure IS the one in the figure. Everything else — colour inversion, histograms, chi2 — can
agree while the grid is offset, and an offset grid is what produced the first version of this
extraction (it sampled `unit_xyz[:, 1, :]`, a tessel CORNER, 1.7 facet widths from the centre).

Facets that straddle the +-180 deg seam are drawn as both branches, since one of the two is
always off the plot.
"""
import os
import numpy as np
from PIL import Image, ImageDraw
import importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
ZET  = os.path.join(HERE, "zetAnd")
SCALE = 3

spec = importlib.util.spec_from_file_location("ex", os.path.join(HERE, "zetand_extract_map.py"))
ex = importlib.util.module_from_spec(spec); spec.loader.exec_module(ex)

def facets(nside_exp=3):
    """Facet corners and centres as (npix, 5, 2) in degrees, from ROTIR's own tessellation."""
    path = os.path.join(ZET, f"facets_{nside_exp}.txt")
    if not os.path.exists(path):
        raise SystemExit(f"missing {path}; write it from Julia with\n"
                         "  using ROTIR, DelimitedFiles\n"
                         "  t = tessellation_healpix(3)\n"
                         "  rows = [(i, k, atand(t.unit_xyz[i,k,2], t.unit_xyz[i,k,1]),\n"
                         "           asind(clamp(t.unit_xyz[i,k,3],-1,1)))\n"
                         "          for i in 1:t.npix for k in 1:5]\n"
                         "  writedlm(path, [r[j] for r in rows, j in 1:4])")
    F = np.loadtxt(path)
    n = int(F[:, 0].max())
    return F, n

def main():
    im_arr = ex.raster()
    (xc, yc, a, b), rms = ex.fit_ellipse(im_arr)
    print(f"ellipse aspect {a/b:.4f} (Mollweide 2.0), rms {rms:.2f} px")
    def project(lon, lat):
        xn, yn = ex.mollweide(np.asarray(lon, float), np.asarray(lat, float))
        return xc + xn*a, yc - yn*b

    orig = Image.open(os.path.join(ZET, "zetAnd2019A_orig.png")).convert("RGB")
    big = orig.resize((orig.width*SCALE, orig.height*SCALE), Image.LANCZOS)
    draw = ImageDraw.Draw(big)
    F, n = facets()
    for i in range(1, n+1):
        rows = F[F[:, 0] == i]
        corners = rows[rows[:, 1] <= 4]
        lon, lat = corners[:, 2].copy(), corners[:, 3]
        branches = ([np.where(lon < 0, lon+360, lon), np.where(lon > 0, lon-360, lon)]
                    if lon.max() - lon.min() > 180 else [lon])
        for L in branches:
            X, Y = project(np.clip(L, -180, 180), lat)
            draw.polygon([(x*SCALE, y*SCALE) for x, y in zip(X, Y)],
                         outline=(0, 160, 255), width=2)
        c = rows[rows[:, 1] == 5][0]
        cx, cy = project([c[2]], [c[3]])
        draw.ellipse([cx[0]*SCALE-4, cy[0]*SCALE-4, cx[0]*SCALE+4, cy[0]*SCALE+4],
                     fill=(0, 255, 0))
    out = os.path.join(ZET, "zetAnd2019A_facets.png")
    big.resize((big.width//2, big.height//2), Image.LANCZOS).save(out)
    print("wrote", out)
    tip = big.crop((int((xc-a)*SCALE)-20, int((yc-55)*SCALE), int((xc-a+75)*SCALE), int((yc+55)*SCALE)))
    out = os.path.join(ZET, "zetAnd2019A_facets_tip.png")
    tip.save(out); print("wrote", out)

    # the same tip, original beside our re-render
    def crop_tip(path):
        im = Image.open(path).convert("RGB"); A = np.asarray(im).astype(int)
        ink = A.sum(2) < 700
        ys = np.where(ink.any(1))[0]
        rows = [(y, np.where(ink[y])[0]) for y in ys if ink[y].sum() > 50]
        y_eq, xr = max(rows, key=lambda t: t[1].max()-t[1].min())
        box = (max(0, xr.min()-6), max(0, y_eq-70),
               min(im.width, xr.min()+110), min(im.height, y_eq+70))
        return im.crop(box).resize(((box[2]-box[0])*5, (box[3]-box[1])*5), Image.NEAREST)
    L = crop_tip(os.path.join(ZET, "zetAnd2019A_orig.png"))
    R = crop_tip(os.path.join(ZET, "zetAnd2019A_map_rotir.png"))
    cv = Image.new("RGB", (L.width + R.width + 30, max(L.height, R.height)), "white")
    cv.paste(L, (0, 0)); cv.paste(R, (L.width+30, 0))
    out = os.path.join(ZET, "zetAnd2019A_tip_compare.png")
    cv.save(out); print("wrote", out)

if __name__ == "__main__":
    main()
