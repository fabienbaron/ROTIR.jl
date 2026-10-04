#!/usr/bin/env python3
"""Step 1 of 2 — measure Figure 10 of arXiv:2603.18144 and export its orthographic panels.

    python3 demos/zetand_fig10_panels.py          # measure + export
    julia --project=. demos/zetand_fig10_match.jl # step 2: match, then come back here
    python3 demos/zetand_fig10_panels.py --figure # build the comparison figure

WHAT THIS ESTABLISHES, from the published figure alone, with no data and no assumption about
SURFING beyond what Table 6 prints:

  * the silhouette is PHASE-DEPENDENT -- semi-major axis 251.8, 245.4, 238.8, 245.2 px at
    phi = 0.000, 0.125, 0.250, 0.375 while the semi-minor stays at 236.9. A rotationally
    oblate spheroid projects to a FIXED ellipse, so the figure is not one. The period-0.5
    modulation with a constant minor axis is a tidally elongated (prolate) figure, which is
    what the paper's own text describes: the light curve is "dominated by the star's
    ellipsoidal shape". zeta And is a synchronised RS CVn binary, so this is expected.

  * fitting the orthographic shadow of an ellipsoid to the eight measured second-moment
    matrices (24 numbers, 6 parameters) reproduces them to 0.1 px out of 240:

        radius_x : radius_y : radius_z  =  1.0629 : 1.0008 : 1
        inclination                     =  69.44 deg       (Table 6: 69.3 +/- 6.9)
        major-to-minor ratio            =  1.0629          (Table 6: 1.063 +/- 0.013)
        rotational flattening           =  1.0008          i.e. NONE -- a prolate spheroid

  * the plate scale comes from the figure's own x ticks, 7 of them at 0.5 mas spacing and
    94.917 px apart, NOT from the axis limits (the +/-1.5 ticks sit 22 px inside the spines):

        2 * radius_z = 2.4960 mas       (Table 6 theta_LD: 2.495 +/- 0.014)

    So their quoted theta_LD is the MINOR axis, and the major axis is 1.063 times larger.

CAVEATS. Everything here is read off a rasterised PDF page at 400 dpi. The silhouette fit is
good to ~0.1 px but the colour values are read through a lookup of the figure's colourbar, and
the orthographic panels SATURATE against its 4600 K ceiling over most of the disc, so per-pixel
values there are quantisation-limited. That is why step 2 registers on facet EDGES, which need
no colour calibration at all.
"""
import os, sys, subprocess
import numpy as np
from PIL import Image
from scipy.spatial import cKDTree
from scipy.ndimage import binary_fill_holes, binary_erosion, label, gaussian_filter
from scipy.ndimage import shift as ndshift
import scipy.optimize as so

HERE = os.path.dirname(os.path.abspath(__file__))
ZET  = os.path.join(HERE, "zetAnd")
PDF  = os.path.join(ZET, "Roettenbacher_zetAnd_arXiv2603.18144.pdf")
PAGE, DPI, PAD, NX = 13, 400, 60, 641
VMIN, VMAX = 3500.0, 4600.0
# panel boxes at 400 dpi, (row0, row1, col0, col1): tight bounding boxes of the orange blobs.
# Everything below pads them by PAD and re-segments, so a tight box cannot truncate a fit.
PANELS = {"0.000": (1574,2064,513,994),  "0.125": (1574,2064,1177,1660),
          "0.250": (1574,2064,1844,2320),"0.375": (1574,2064,2508,2981),
          "0.500": (2235,2726,513,994),  "0.625": (2235,2726,1177,1660),
          "0.750": (2235,2726,1844,2320),"0.875": (2235,2726,2508,2981)}
TICK_PX  = 94.9167        # pixels per 0.5 mas, measured from the bottom-left panel's x ticks
SCALE    = 0.5/TICK_PX    # mas per pixel
THETA_LD = 2.495          # Table 6, mas


def page(tmp="/tmp/zetand_fig10"):
    if not os.path.exists(f"{tmp}-{PAGE}.png"):
        subprocess.run(["pdftoppm", "-r", str(DPI), "-f", str(PAGE), "-l", str(PAGE),
                        "-png", PDF, tmp], check=True)
    return np.asarray(Image.open(f"{tmp}-{PAGE}.png").convert("RGB")).astype(float)


def colourbar(im):
    """The ramp row between the Mollweide and the panels: wide, dark left, light right."""
    H = im.shape[0]
    for y in range(int(0.20*H), int(0.40*H)):
        row = im[y]; nz = np.where(row.sum(1) < 740)[0]
        if len(nz) < 300: continue
        x0, x1 = nz.min(), nz.max()
        if x1 - x0 < 400: continue
        seg = row[x0:x1+1].sum(1)
        if seg[:50].mean() < 150 and seg[-50:].mean() > 600:
            return im[y, x0+6:x1-5]
    raise SystemExit("colourbar not found")


def disc(p):
    """The star's silhouette: the largest connected coloured component, holes filled."""
    m = binary_fill_holes((p.max(2) - p.min(2)) > 25)
    lab, n = label(m)
    s = np.bincount(lab.ravel()); s[0] = 0
    return lab == s.argmax()


def moments(mask):
    ys, xs = np.nonzero(mask)
    c = np.cov(np.vstack([xs - xs.mean(), ys - ys.mean()]))
    return np.array([c[0,0], c[1,1], c[0,1]]), (ys.mean(), xs.mean())


def axes_of(m):
    C = np.array([[m[0], m[2]], [m[2], m[1]]])
    w, v = np.linalg.eigh(C)
    return 2*np.sqrt(w[1]), 2*np.sqrt(w[0]), np.degrees(np.arctan2(v[1,1], v[0,1])) % 180


def rot_vertex(a1, a2, a3):
    """ROTIR's `rot_vertex`, so a fitted angle means what `rotate_star` means by it."""
    c1,s1,c2,s2,c3,s3 = np.cos(a1),np.sin(a1),np.cos(a2),np.sin(a2),np.cos(a3),np.sin(a3)
    return np.array([[-s1*c2*s3+c1*c3,  s1*c3*c2+c1*s3, -s1*s2],
                     [-c1*c2*s3-s1*c3,  c1*c3*c2-s1*s3, -c1*s2],
                     [         -s2*s3,           s2*c3,     c2]])


def predict(th, phases, sgn, sig):
    """Second moments of an ellipsoid's orthographic shadow.

    ROTIR right-multiplies, s = R' x, so the sky shape matrix is W = R' M R with
    M = diag(rx^2, ry^2, rz^2); the outline is the ellipse of W's (west, north) block, and a
    uniformly filled ellipse has pixel second moments exactly W_2x2 / 4.
    """
    rx, ry, rz, inc, pa, dphi = th
    M = np.diag([rx**2, ry**2, rz**2]); out = []
    for ph in phases:
        R = rot_vertex(2*np.pi*(sgn*ph + dphi), np.radians(inc), np.radians(pa))
        W = R.T @ M @ R
        out.append([W[0,0]/4, W[1,1]/4, sig*W[0,1]/4])
    return np.array(out)


def measure(im):
    """Silhouettes of the eight panels, and the ellipsoid that reproduces all of them."""
    phases = np.array([float(k) for k in PANELS])
    obs, cen = [], []
    print("published silhouettes (400 dpi pixels)")
    print("   phi        a        b      b/a      PA")
    for name, (r0,r1,c0,c1) in PANELS.items():
        m = disc(im[r0-PAD:r1+PAD, c0-PAD:c1+PAD])
        mo, ct = moments(m); obs.append(mo); cen.append(ct)
        a, b, pa = axes_of(mo)
        print(f"  {name}  {a:7.2f}  {b:7.2f}  {b/a:.4f}  {pa:6.2f}")
    obs = np.array(obs)
    best = None
    for sgn in (1, -1):
        for sig in (1, -1):
            for d0 in np.arange(0, 1, 0.05):
                for i0 in (40.0, 69.3, 85.0):
                    for pa0 in (34.0, 124.0):
                        f = lambda t: (predict(t, phases, sgn, sig) - obs).ravel()
                        s = so.least_squares(f, [250.,240.,237.,i0,pa0,d0],
                                             bounds=([100,100,100,0.1,-360,-2],
                                                     [400,400,400,89.9,360,2]))
                        if best is None or s.cost < best[0].cost:
                            best = (s, sgn, sig)
    s, sgn, sig = best
    rx, ry, rz, inc, pa, dphi = s.x
    print(f"\nellipsoid fitted to all eight (rms {np.sqrt(np.mean(s.fun**2)):.2f} px^2):")
    print(f"  rx:ry:rz = {rx/rz:.4f} : {ry/rz:.4f} : 1     "
          f"major/minor = {rx/rz:.4f}   (Table 6: 1.063 +/- 0.013)")
    print(f"  inclination {inc:.2f} deg                     (Table 6: 69.3 +/- 6.9)")
    print(f"  |PA| {pa % 180:.2f} deg or {180 - pa % 180:.2f} deg -- the silhouette cannot")
    print(f"       distinguish a reflection; step 2 does.    (Table 6: 124.0 +/- 2.9)")
    print(f"\nplate scale {SCALE*1e3:.4f} muas/px, from the figure's own 0.5 mas x ticks:")
    for n, v in (("2*rx (major)", rx), ("2*ry", ry), ("2*rz (minor)", rz)):
        print(f"  {n:14s} = {2*v*SCALE:.4f} mas")
    print(f"  Table 6 theta_LD = {THETA_LD} +/- 0.014 mas  -- i.e. their theta_LD is the MINOR axis")
    return np.array(cen), (rx, ry, rz, inc, pa)


def export(im):
    """Panel values and facet-edge maps, each centred on its own disc."""
    B = colourbar(im); vals = np.linspace(VMIN, VMAX, len(B)); bt = cKDTree(B)
    V = np.full((8, NX, NX), np.nan, np.float32)
    M = np.zeros((8, NX, NX), bool)
    E = np.zeros((8, NX, NX), np.float32)
    for k, (name, (r0,r1,c0,c1)) in enumerate(PANELS.items()):
        p = im[r0-PAD:r1+PAD, c0-PAD:c1+PAD]
        m = disc(p); ys, xs = np.nonzero(m); yc, xc = ys.mean(), xs.mean()
        d, j = bt.query(p.reshape(-1,3), k=1)
        v = np.where(m, vals[j.reshape(p.shape[:2])], np.nan)
        gy = np.abs(np.diff(p, axis=0, prepend=p[:1])).sum(2)
        gx = np.abs(np.diff(p, axis=1, prepend=p[:,:1])).sum(2)
        e = np.hypot(gx, gy); e[~binary_erosion(m, iterations=4)] = 0.0

        iy, ix, h = int(round(yc)), int(round(xc)), NX//2
        pad = ((max(0,h-iy), max(0,iy+h+1-p.shape[0])), (max(0,h-ix), max(0,ix+h+1-p.shape[1])))
        iyp, ixp = iy + pad[0][0], ix + pad[1][0]
        cut = lambda a, fill: np.pad(a, pad, constant_values=fill)[iyp-h:iyp+h+1, ixp-h:ixp+h+1]
        dy, dx = round(yc)-yc, round(xc)-xc
        vc = ndshift(np.nan_to_num(cut(v, np.nan)), (dy,dx), order=1, cval=0.0)
        mc = ndshift(cut(m.astype(float), 0.0),     (dy,dx), order=1, cval=0.0)
        ec = ndshift(cut(e, 0.0),                   (dy,dx), order=1, cval=0.0)
        good = mc > 0.999
        V[k] = np.where(good, vc, np.nan); M[k] = good; E[k] = gaussian_filter(ec, 2.0)
    for nm, a in (("V", V), ("M", M.astype(np.uint8)), ("E", E)):
        a.astype(a.dtype).tofile(os.path.join(ZET, f"fig10_panels_{nm}.bin"))
    np.save(os.path.join(ZET, "fig10_panels_V.npy"), V)
    np.save(os.path.join(ZET, "fig10_panels_M.npy"), M)
    print(f"\nexported eight {NX}x{NX} panels to {ZET}/fig10_panels_*.bin")
    print(f"  NOTE the panels saturate: interiors sit at 4530-4599 K against the 4600 K ceiling,")
    print(f"  so step 2 registers on facet edges rather than on these values.")


def figure():
    """Step 3: published vs ROTIR vs residual, once step 2 has dumped the model panels."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    need = os.path.join(ZET, "fig10_model_panels.bin")
    if not os.path.exists(need):
        raise SystemExit(f"run `julia --project=. demos/zetand_fig10_match.jl` first ({need})")
    V = np.load(os.path.join(ZET, "fig10_panels_V.npy"))
    M = np.load(os.path.join(ZET, "fig10_panels_M.npy"))
    rd = lambda f: np.fromfile(os.path.join(ZET, f), dtype=np.float32
                               ).reshape(NX, NX, 8, order="F").transpose(2,0,1)
    # the render has row ~ North; Figure 10 has North DOWN (its y axis runs -1.5 at the top to
    # +1.5 at the bottom), so the model's rows reverse into the published frame.
    m  = rd("fig10_model_panels.bin")[:, ::-1, :]
    mu = rd("fig10_model_mu.bin")[:, ::-1, :]
    phi = list(PANELS)
    good = M & (m > 0)
    fig, ax = plt.subplots(3, 8, figsize=(21, 8.4))
    for k in range(8):
        th = np.where(M[k], V[k], np.nan); ou = np.where(m[k] > 0, m[k], np.nan)
        d  = np.where(good[k], V[k]-m[k], np.nan)
        for r, (img, lo, hi, cm) in enumerate(((th,3500,4600,"gist_heat"),
                                               (ou,3500,4600,"gist_heat"),
                                               (d,-60,60,"RdBu_r"))):
            h = ax[r,k].imshow(img, vmin=lo, vmax=hi, cmap=cm, origin="upper")
            ax[r,k].set_xticks([]); ax[r,k].set_yticks([])
            if k == 7 and r == 2:
                fig.colorbar(h, ax=ax[2,:], fraction=0.012, pad=0.01, label="residual (K)")
        ax[0,k].set_title(f"$\\phi$ = {phi[k]}", fontsize=12)
    for r, lab in enumerate(("Roettenbacher et al.\nFigure 10 (published)",
                             "ROTIR, rendered from\ntheir Figure 10 map", "published $-$ ROTIR")):
        ax[r,0].set_ylabel(lab, fontsize=10)
    fig.suptitle("$\\zeta$ And Figure 10 reproduced from its own map — prolate 1.0629:1.0008:1, "
                 "$2r_{\\rm minor}$ = 2.496 mas, $i$ = 69.44$^\\circ$, PA = 124.3$^\\circ$, "
                 "$I \\propto T\\mu^{0.269}$", fontsize=13)
    fig.tight_layout(rect=[0,0,1,0.95])
    out = os.path.join(ZET, "zetAnd_fig10_ortho_compare.png")
    fig.savefig(out, dpi=110)
    res = (V-m)[good]; mus = mu[good]
    print(f"wrote {out}")
    print("\nresidual against limb angle -- the test that the projection is right:")
    print("   mu range        n      mean(K)   rms(K)")
    for lo, hi in ((0.95,1.0),(0.9,0.95),(0.8,0.9),(0.7,0.8),(0.6,0.7),(0.5,0.6),
                   (0.4,0.5),(0.3,0.4),(0.2,0.3),(0.1,0.2),(0.0,0.1)):
        s = (mus >= lo) & (mus < hi)
        if s.sum() > 100:
            print(f"   {lo:.2f}-{hi:.2f}  {s.sum():8d}  {res[s].mean():+8.2f}  {res[s].std():7.2f}")
    print(f"\n   all pixels: rms {res.std():.2f} K, correlation "
          f"{np.corrcoef(V[good], m[good])[0,1]:.5f}")


if __name__ == "__main__":
    if "--figure" in sys.argv:
        figure()
    else:
        im = page()
        measure(im)
        export(im)
        print("\nnext: julia --project=. demos/zetand_fig10_match.jl")
