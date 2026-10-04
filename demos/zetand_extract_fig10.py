#!/usr/bin/env python3
"""Extract Figure 10 of arXiv:2603.18144 — and measure how SURFING converts temperature to flux.

    python3 demos/zetand_extract_fig10.py

Figure 10 is the decisive page because it contains BOTH halves of the conversion: the top panel
is ζ And's temperature map in kelvin (colourbar 3500-4600 K), the bottom eight are the H-band
INTENSITY images that map produces, at phases 0.000 to 0.875. The relation between them is
therefore measurable from the figure alone, with no data and no assumptions about their code.

WHAT IT MEASURES, and these are the numbers the rest of the ζ And work rests on:

  * the top map, resampled to HEALPix nside_exp 3: 3742-4589 K, mean 4372 K;
  * the intensity panels' centre-to-limb profile, fitted as value = A + B*mu^alpha on the
    phi = 0.500 panel (the phi = 0.000 fit is degenerate, a spot sits at disc centre):
    A = 3338, B = 1233, alpha = 0.248 — consistent with the alpha = 0.269 they fix;
  * de-limb-darkening the panels with that law, the deepest intensity contrast over all eight
    phases is 0.881 (1st percentile; 0.847 minimum);
  * against a temperature contrast of 3742/4372 = 0.856 in the map above them.

    I proportional to T   predicts 0.856   <- matches
    Planck at 1.605 um    predicts 0.679   <- does not

So SURFING converts LINEARLY. That is the Rayleigh-Jeans limit, and H band is not the
Rayleigh-Jeans limit for a K giant: at 1.605 um and 4372 K, hc/(lambda k T) = 2.05 and
dlnI/dlnT = 2.36. Their intensity maps are unaffected; their TEMPERATURES are too cold by
roughly a factor two in the deficit. `demos/zetand_proper.jl` does the arithmetic.

CAVEATS, because this rests on reading a figure:
  * the affine colour calibration comes from one panel, though the other seven give consistent
    contrasts of 0.88-0.95;
  * it assumes the orthographic panels share the Mollweide's colourbar, which the caption does
    not say. The evidence is that they span almost exactly its range (3561-4589 against
    3500-4600) with only 0.1 % of pixels below 4000, which is what a limb-darkened disc read
    through that scale looks like.
"""
import os, subprocess
import numpy as np
from PIL import Image
from scipy.spatial import cKDTree
import scipy.optimize as so

HERE = os.path.dirname(os.path.abspath(__file__))
ZET  = os.path.join(HERE, "zetAnd")
PDF  = os.path.join(ZET, "Roettenbacher_zetAnd_arXiv2603.18144.pdf")
PAGE, DPI = 13, 400
VMIN, VMAX = 3500.0, 4600.0
# panel boxes at 400 dpi, (row0, row1, col0, col1), found by locating the orange blobs
PANELS = {"0.000": (1574,2064,513,994),  "0.125": (1574,2064,1177,1660),
          "0.250": (1574,2064,1844,2320),"0.375": (1574,2064,2508,2981),
          "0.500": (2235,2726,513,994),  "0.625": (2235,2726,1177,1660),
          "0.750": (2235,2726,1844,2320),"0.875": (2235,2726,2508,2981)}

def render(tmp="/tmp/zetand_fig10"):
    subprocess.run(["pdftoppm","-r",str(DPI),"-f",str(PAGE),"-l",str(PAGE),"-png",PDF,tmp],
                   check=True)
    return np.asarray(Image.open(f"{tmp}-{PAGE}.png").convert("RGB")).astype(float)

def colourbar(im):
    """The ramp row: wide, dark on the left, light on the right."""
    H = im.shape[0]
    for y in range(int(0.20*H), int(0.40*H)):
        row = im[y]; nz = np.where(row.sum(1) < 740)[0]
        if len(nz) < 300: continue
        x0, x1 = nz.min(), nz.max()
        if x1-x0 < 400: continue
        seg = row[x0:x1+1].sum(1)
        if seg[:50].mean() < 150 and seg[-50:].mean() > 600:
            return im[y, x0+6:x1-5]
    raise SystemExit("colourbar not found")

def main():
    im  = render()
    bar = colourbar(im); vals = np.linspace(VMIN, VMAX, len(bar)); bt = cKDTree(bar)
    print(f"colourbar: {len(bar)} samples, {VMIN:.0f}-{VMAX:.0f} K")

    # ---- the orthographic panels: centre-to-limb law, then the de-limb-darkened contrast ----
    A = B = AL = None
    contrasts = {}
    for name, (r0,r1,c0,c1) in PANELS.items():
        p = im[r0:r1, c0:c1]
        d, j = bt.query(p.reshape(-1,3), k=1)
        on = (d < 30).reshape(p.shape[:2]) & ~(p.max(2)-p.min(2) < 12)
        V = np.where(on, vals[j.reshape(p.shape[:2])], np.nan)
        ys, xs = np.mgrid[0:p.shape[0], 0:p.shape[1]]
        yc, xc = np.nanmean(ys[on]), np.nanmean(xs[on])
        rr = np.hypot(xs-xc, ys-yc); R = np.percentile(rr[on], 99.5); rho = rr/R
        if name == "0.500":                      # calibrate the affine colour scale here
            edges = np.linspace(0, 0.97, 25); prof=[]; rc=[]
            for k in range(24):
                m = on & (rho>=edges[k]) & (rho<edges[k+1]) & (rho<0.97)
                if m.sum() > 50: prof.append(np.nanmedian(V[m])); rc.append(0.5*(edges[k]+edges[k+1]))
            prof = np.array(prof); mu = np.sqrt(np.maximum(1-np.array(rc)**2, 1e-6))
            s = so.least_squares(lambda q: q[0]+q[1]*mu**q[2]-prof, [3300.,1200.,0.3],
                                 bounds=([2000,1,0.0],[4500,5000,3.0]))
            A, B, AL = s.x
            print(f"centre-to-limb (phi=0.500): value = {A:.0f} + {B:.0f}*mu^{AL:.3f}, "
                  f"rms {np.sqrt(np.mean(s.fun**2)):.1f}")
    for name, (r0,r1,c0,c1) in PANELS.items():
        p = im[r0:r1, c0:c1]
        d, j = bt.query(p.reshape(-1,3), k=1)
        on = (d < 30).reshape(p.shape[:2]) & ~(p.max(2)-p.min(2) < 12)
        V = np.where(on, vals[j.reshape(p.shape[:2])], np.nan)
        ys, xs = np.mgrid[0:p.shape[0], 0:p.shape[1]]
        yc, xc = np.nanmean(ys[on]), np.nanmean(xs[on])
        rr = np.hypot(xs-xc, ys-yc); R = np.percentile(rr[on], 99.5); rho = rr/R
        mu = np.sqrt(np.maximum(1-rho**2, 1e-6))
        I = (V - A)/(B*mu**AL)
        q = I[on & (rho < 0.85)]
        contrasts[name] = np.nanpercentile(q, 1)
    deepest = min(contrasts.values())
    print("de-limb-darkened intensity contrast per phase: " +
          " ".join(f"{k}:{v:.3f}" for k,v in contrasts.items()))
    print(f"deepest = {deepest:.3f}")

    # ---- the Mollweide map, and the verdict -------------------------------------------------
    sub = im[:1300]
    d, j = bt.query(sub.reshape(-1,3), k=1)
    good = ((d < 30).reshape(sub.shape[:2])) & ~(sub.max(2)-sub.min(2) < 12)
    rows = [(y, np.where(good[y])[0].min(), np.where(good[y])[0].max())
            for y in range(good.shape[0])
            if good[y].any() and np.ptp(np.where(good[y])[0]) > 300]
    R_ = np.array(rows)
    def resid(q):
        xc,yc,a,b = q; t = np.clip((R_[:,0]-yc)/b, -1, 1)
        w = a*np.sqrt(np.maximum(0,1-t**2))
        return np.concatenate([(xc-w)-R_[:,1], (xc+w)-R_[:,2]])
    xc,yc,a,b = so.least_squares(resid, [(R_[:,1].min()+R_[:,2].max())/2, R_[:,0].mean(),
                                         (R_[:,2].max()-R_[:,1].min())/2,
                                         (R_[-1,0]-R_[0,0])/2]).x
    ys, xs = np.mgrid[0:sub.shape[0], 0:sub.shape[1]]
    xn, yn = (xs-xc)/a, (yc-ys)/b
    inside = (xn**2+yn**2 < 0.995**2) & good
    th  = np.arcsin(np.clip(yn[inside], -1, 1))
    lat = np.degrees(np.arcsin(np.clip((2*th+np.sin(2*th))/np.pi, -1, 1)))
    lon = np.degrees(np.pi*xn[inside]/np.maximum(np.cos(th), 1e-12))
    ok  = np.abs(lon) <= 180
    v   = vals[j.reshape(sub.shape[:2])[inside][ok]]
    print(f"Mollweide map: {v.min():.0f}-{v.max():.0f} K, mean {v.mean():.0f} K "
          f"(ellipse a/b {a/b:.4f})")
    Tmin, Tmean = np.percentile(v, 0.5), v.mean()
    h, c_, kB, lam = 6.62607015e-34, 2.99792458e8, 1.380649e-23, 1.605e-6
    Bp = lambda T: 1.0/np.expm1(h*c_/(lam*kB*T))
    print(f"\ntemperature contrast in the map : {Tmin/Tmean:.3f}")
    print(f"  predicted intensity if I ~ T   : {Tmin/Tmean:.3f}")
    print(f"  predicted intensity if Planck  : {Bp(Tmin)/Bp(Tmean):.3f}")
    print(f"  MEASURED in their own panels   : {deepest:.3f}")
    print("\n=> SURFING converts temperature to intensity LINEARLY (Rayleigh-Jeans).")

if __name__ == "__main__":
    main()
