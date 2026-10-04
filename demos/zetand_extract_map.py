#!/usr/bin/env python3
"""Read the ζ And surface map off the figure in `demos/zetAnd/zetAnd2019A.pdf`.

The PDF holds a single 836x514 raster: a Mollweide projection in matplotlib's `afmhot`
colormap with a colourbar labelled 0.8025 to 0.9977. Nothing vector, so the values have to come
back out through the colours.

    python3 demos/zetand_extract_map.py            # writes demos/zetAnd/zetAnd2019A_map.txt

Needs numpy, PIL, scipy, and `pdfimages` (poppler). The result is resampled onto HEALPix
nside_exp = 3 (768 tessels), which is close to the facet size visible in the figure.

HOW IT IS KNOWN TO WORK, since a colour inversion is easy to get subtly wrong:

  * the colourbar ramp is strictly monotone, so RGB -> value is single valued, and the sampled
    pixels sit within 1/255 of it (median colour distance 1.0 over 237 758 pixels);
  * the plot boundary fits an ellipse of axis ratio 1.9926 with 0.29 px rms, against the 2.0
    Mollweide requires — so the projection and the plot box are both confirmed;
  * the resampled HEALPix histogram matches the image's own pixel histogram (mean 0.9493 vs
    0.9459, sd 0.0385 vs 0.0419, median 0.9598 vs 0.9583), which is the check that the
    resampling did not distort the distribution;
  * graticule contamination is bounded rather than assumed: see the comment in the sampling
    loop. Purging grey pixels moves the map by rms 0.0035 of a 0.195 range, and by 0.085 for
    the one tessel that lies along the solid equator;
  * and the map FITS ζ And's visibilities, which no artefact would. See
    `demos/zetand_published_map.jl`: of the four latitude/longitude parities, exactly one
    reaches chi2/n = 2.29 against 3.63 for a uniform surface on the same figure, the other
    three being WORSE than uniform (3.89, 4.20, 4.37).

LATITUDE IS FLIPPED relative to the plot. That is the parity the data select, and it is what
this script writes; whether the flip lives in the figure's convention or in ROTIR's is not
something the data can say.
"""
import subprocess, sys, os
import numpy as np
from PIL import Image
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
PDF  = os.path.join(HERE, "zetAnd", "zetAnd2019A.pdf")
OUT  = os.path.join(HERE, "zetAnd", "zetAnd2019A_map.txt")
VMIN, VMAX = 0.8025, 0.9977          # printed under the colourbar
NSIDE_EXP = 3

def raster(tmp="/tmp/zetand_extract"):
    subprocess.run(["pdfimages", "-png", PDF, tmp], check=True)
    return np.asarray(Image.open(tmp + "-000.png").convert("RGB")).astype(float)

def fit_ellipse(im):
    """The Mollweide boundary, from the rows wide enough to be the plot and not the title."""
    ink = (im[:455].sum(2) < 700)
    rows = []
    for y in range(ink.shape[0]):
        xs = np.where(ink[y])[0]
        if len(xs) and xs.max() - xs.min() > 300:
            rows.append((y, xs.min(), xs.max()))
    rows = np.array(rows)
    def resid(p):
        xc, yc, a, b = p
        t = np.clip((rows[:, 0] - yc)/b, -1, 1)
        w = a*np.sqrt(np.maximum(0, 1 - t**2))
        return np.concatenate([(xc - w) - rows[:, 1], (xc + w) - rows[:, 2]])
    p0 = [(rows[:, 1].min() + rows[:, 2].max())/2, rows[:, 0].mean(),
          (rows[:, 2].max() - rows[:, 1].min())/2, (rows[-1, 0] - rows[0, 0])/2]
    sol = least_squares(resid, p0)
    return sol.x, np.sqrt(np.mean(sol.fun**2))

def mollweide(lon_deg, lat_deg):
    """Normalised to the unit box: x in [-1,1] across the ellipse, y in [-1,1] pole to pole."""
    lam, phi = np.radians(lon_deg), np.radians(lat_deg)
    th = phi.copy()
    for _ in range(60):                       # Newton on 2θ + sin 2θ = π sin φ
        th -= (2*th + np.sin(2*th) - np.pi*np.sin(phi)) / (2 + 2*np.cos(2*th) + 1e-12)
    return lam*np.cos(th)/np.pi, np.sin(th)

def main():
    lonlat_file = os.path.join(HERE, "zetAnd", "tessel_lonlat_3.txt")
    if not os.path.exists(lonlat_file):
        sys.exit(f"missing {lonlat_file}; write it from Julia with\n"
                 "  using ROTIR, DelimitedFiles\n"
                 "  t = tessellation_healpix(3); c = t.unit_xyz[:, 1, :]\n"
                 "  writedlm(path, hcat([atand(c[i,2], c[i,1]) for i in 1:t.npix],\n"
                 "                      [asind(clamp(c[i,3], -1, 1)) for i in 1:t.npix]))")
    im = raster()
    (xc, yc, a, b), rms = fit_ellipse(im)
    print(f"ellipse xc={xc:.2f} yc={yc:.2f} a={a:.2f} b={b:.2f} "
          f"aspect {a/b:.4f} (Mollweide 2.0) rms {rms:.2f} px")
    bar  = im[470, 215:622]                   # the colourbar ramp, left = VMIN
    vals = np.linspace(VMIN, VMAX, len(bar))
    lonlat = np.loadtxt(lonlat_file)
    xn, yn = mollweide(lonlat[:, 0], -lonlat[:, 1])      # LATITUDE FLIPPED — see the docstring
    px, py = xc + xn*a, yc - yn*b
    out = np.full(len(px), np.nan)
    for i, (X, Y) in enumerate(zip(px, py)):
        xs = np.clip(np.arange(int(X)-3, int(X)+4), 0, im.shape[1]-1)
        ys = np.clip(np.arange(int(Y)-3, int(Y)+4), 0, im.shape[0]-1)
        patch = im[np.ix_(ys, xs)].reshape(-1, 3)
        d = ((patch[:, None, :] - bar[None, :, :])**2).sum(2)
        j = d.argmin(1)
        good = np.sqrt(d[np.arange(len(patch)), j]) < 25
        # THE COLOUR FILTER DOES NOT CATCH THE GRATICULE. `afmhot` starts at pure black, so a
        # black gridline pixel sits exactly ON the ramp and inverts to the colourbar floor.
        # Measured: 2.1 % of accepted pixels per tessel are grey on average and 57 % for the
        # worst, which is a tessel lying along the solid equator — its median came out at the
        # floor, 0.8025, where purging gives 0.887. Grey is the signature, because afmhot is
        # grey only at its two endpoints.
        grey = (patch.max(1) - patch.min(1) < 12) & (patch.sum(1) < 690)
        keep = good & ~grey
        keep = keep if keep.sum() >= 5 else good          # all-graticule: take what there is
        if keep.sum() >= 5:
            out[i] = np.median(vals[j][keep])             # median too: dotted lines are 1-2 px
    print(f"sampled {np.isfinite(out).sum()}/{len(out)} tessels, "
          f"range {np.nanmin(out):.4f}..{np.nanmax(out):.4f}, mean {np.nanmean(out):.4f}")
    np.savetxt(OUT, out)
    print("wrote", OUT)

if __name__ == "__main__":
    main()
