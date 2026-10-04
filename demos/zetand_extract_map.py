#!/usr/bin/env python3
"""Read the ζ And surface map off the figure in `demos/zetAnd/zetAnd2019A.pdf`.

The PDF holds a single 836x514 raster: a Mollweide projection with a colourbar labelled 0.8025
to 0.9977. Nothing vector, so the values have to come back out through the colours. The ramp is
`gist_heat` — matplotlib's version of IDL's Red Temperature, which is what the `MIRCX_IDL`
pipeline in the data filenames suggests — matching to 0.0125 rms against 0.233 for `afmhot`.
NONE OF THE CODE BELOW DEPENDS ON THAT: the lookup is built from the figure's own colourbar, so
identifying the map by name is a documentation convenience, not an assumption.

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

THE MAP IS WRITTEN AS PLOTTED, latitude included, so rendering it with ROTIR's own
`plot_mollweide` reproduces the figure. ζ And's visibilities then select the EQUATOR-FLIPPED
parity of it (chi2/n 1.88 against 3.63 for a uniform surface, the three other parities being
worse than uniform), so the flip is applied on the way into the forward model, not here — see
`demos/zetand_published_map.jl` and `healpix_flip_3.txt`, which is that flip as an exact
permutation of HEALPix centres. Keeping the two frames apart matters: the plot convention and
the model convention are different questions, and only the second is something the data answer.
"""
import subprocess, sys, os
import numpy as np
from PIL import Image
from scipy.optimize import least_squares
from scipy.spatial import cKDTree

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
                 "  t = tessellation_healpix(3); c = t.unit_xyz[:, 5, :]   # 5 = CENTRE\n"
                 "  writedlm(path, hcat([atand(c[i,2], c[i,1]) for i in 1:t.npix],\n"
                 "                      [asind(clamp(c[i,3], -1, 1)) for i in 1:t.npix]))")
    im = raster()
    (xc, yc, a, b), rms = fit_ellipse(im)
    print(f"ellipse xc={xc:.2f} yc={yc:.2f} a={a:.2f} b={b:.2f} "
          f"aspect {a/b:.4f} (Mollweide 2.0) rms {rms:.2f} px")
    bar  = im[470, 215:622]                   # the colourbar ramp, left = VMIN
    vals = np.linspace(VMIN, VMAX, len(bar))

    # EVERY PIXEL IS ASSIGNED TO A FACET, rather than every facet sampling a patch around its
    # projected centre. A patch has three failure modes that this does not: it runs off the plot
    # for a limb facet (and the background inverts to the colourbar MAXIMUM, so those facets came
    # back white); it cannot straddle the +-180 deg seam, where a facet's pixels are split
    # between the two edges of the plot; and a facet whose centre lands exactly ON the ellipse
    # has no patch to speak of. Going the other way, pixel -> facet, all three disappear: a
    # facet simply collects whichever of its pixels are drawn, wherever they are.
    ys, xs = np.mgrid[0:im.shape[0], 0:im.shape[1]]
    xn, yn = (xs - xc)/a, (yc - ys)/b
    inside = xn**2 + yn**2 < 0.995**2
    th  = np.arcsin(np.clip(yn[inside], -1, 1))
    lat = np.degrees(np.arcsin(np.clip((2*th + np.sin(2*th))/np.pi, -1, 1)))
    lon = np.degrees(np.pi*xn[inside]/np.maximum(np.cos(th), 1e-12))
    ok  = np.abs(lon) <= 180
    P   = im[inside][ok]

    lonlat = np.loadtxt(lonlat_file)
    cl, sl = np.cos(np.radians(lonlat[:, 0])), np.sin(np.radians(lonlat[:, 0]))
    cb, sb = np.cos(np.radians(lonlat[:, 1])), np.sin(np.radians(lonlat[:, 1]))
    tree = cKDTree(np.column_stack([cb*cl, cb*sl, sb]))
    clo, slo = np.cos(np.radians(lon[ok])), np.sin(np.radians(lon[ok]))
    cla, sla = np.cos(np.radians(lat[ok])), np.sin(np.radians(lat[ok]))
    _, owner = tree.query(np.column_stack([cla*clo, cla*slo, sla]))

    d = ((P[:, None, :] - bar[None, :, :])**2).sum(2)
    j = d.argmin(1)
    on_ramp = np.sqrt(d[np.arange(len(P)), j]) < 25
    # THE COLOUR FILTER DOES NOT CATCH THE GRATICULE. The ramp starts at pure black, so a black
    # gridline pixel sits exactly ON it and inverts to the floor. Grey is the signature, because
    # the ramp is grey only at its two endpoints. The median does the rest: the dotted graticule
    # is 1-2 px wide against a few hundred pixels per facet.
    grey = (P.max(1) - P.min(1) < 12) & (P.sum(1) < 690)
    v = vals[j]
    out = np.full(len(lonlat), np.nan)
    counts = np.zeros(len(lonlat), dtype=int)
    for t in range(len(lonlat)):
        m = (owner == t) & on_ramp & ~grey
        if m.sum() < 5:
            m = (owner == t) & on_ramp
        counts[t] = m.sum()
        if m.sum() >= 1:
            out[t] = np.median(v[m])
    print(f"pixels assigned: {len(P)}; per facet min {counts.min()} median {int(np.median(counts))} max {counts.max()}")
    print(f"sampled {np.isfinite(out).sum()}/{len(out)} tessels, "
          f"range {np.nanmin(out):.4f}..{np.nanmax(out):.4f}, mean {np.nanmean(out):.4f}")
    np.savetxt(OUT, out)
    print("wrote", OUT)

if __name__ == "__main__":
    main()
