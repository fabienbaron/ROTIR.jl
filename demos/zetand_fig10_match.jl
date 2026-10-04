#!/usr/bin/env julia
# Step 2 of 2 — match ROTIR's renders to Figure 10's orthographic panels, and read off what
# SURFING actually plotted.
#
#     python3 demos/zetand_fig10_panels.py           # first: measure + export the panels
#     julia --project=. demos/zetand_fig10_match.jl   # this file
#     python3 demos/zetand_fig10_panels.py --figure   # then: the comparison figure
#
# Step 1 pinned the shape and the inclination from the silhouettes, which it can do without any
# colour calibration. What it CANNOT do is tell a configuration from its mirror image, because a
# reflected ellipsoid has the same outline. This file settles that, and everything else that is
# discrete, in two passes:
#
#   PASS 1, facet EDGES. The published panels are visibly facetted, so the projected HEALPix
#   mesh is itself a registration signal -- one that needs no colourbar, no limb-darkening law
#   and no map. Correlating the rendered facet outlines against the panels' colour-gradient map
#   gives render_phase = -panel_phase, i.e. their phase runs the other way, with a residual
#   0.5 ambiguity that the mesh's own symmetry cannot break.
#
#   PASS 2, facet VALUES. Assign every published pixel to the facet covering it, take the median
#   published value per facet, and regress it on the map value. This breaks the remaining
#   ambiguity decisively: the winner scores r = 0.991 and the next distinct configuration 0.916.
#
# WHAT COMES OUT. The top four configurations tie exactly -- they are one solution written four
# ways -- and the member that uses the paper's own numbers is:
#
#     position angle  124.30 deg E of N   (Table 6: 124.0 +/- 2.9)
#     inclination      69.44 deg          (Table 6: 69.3 +/- 6.9)
#     figure           prolate 1.0629 : 1.0008 : 1, 2*r_minor = 2.496 mas (theta_LD 2.495)
#     map              used exactly as the Mollweide plots it -- no parity surgery
#     phase            sub-observer longitude = 90 deg - 360 deg * phi
#     orientation      East LEFT and North DOWN: the panels' y axis runs -1.5 at the top to
#                      +1.5 at the bottom, so the vertical axis is inverted from the usual
#                      convention. This is why a naive overlay mirrors.
#
# AND WHAT THEY PLOTTED. Regressing per-facet published value on the map gives
#
#     published = 3496 + 0.2478 * (T * mu^0.269),  rms 13 K
#
# a single affine relation holding across all eight panels -- fitting it panel by panel gives
# A = 3460-3521 and B = 0.2427-0.2540, a 2 % spread, so they did NOT auto-scale each panel and
# one calibration can be quoted for the whole figure. Two consequences:
#
#   * the panels are INTENSITY on their own linear scale. The "Temperature (K)" colourbar sits
#     between the Mollweide and the panels and belongs to the Mollweide; read against it the
#     panel interiors come out at 4530-4599 K, which is saturation, not physics.
#   * the conversion is LINEAR in temperature. Against the same map and geometry, I = T gives
#     r = 0.9909 and rms 13.0 K, Planck at 1.605 um gives 0.9152 and 39.1 K, T^4 gives 0.7795
#     and 60.7 K; the best power law is I ~ T^0.96.
#
#     DO NOT READ THE PLANCK ROW AS "THE CORRECT ANSWER". A blackbody has dlnI/dlnT = 2.35 at
#     4372 K, but a real H-band atmosphere is much shallower — Kurucz ATLAS9 at logg 2.5 gives
#     **1.51**, since H- opacity falls with temperature and partly cancels the Planck
#     steepening. Scored the same way, Kurucz gives r = 0.97597 and rms 21.1 K: still behind
#     linear, but far ahead of a blackbody. Inverting their map through each, their spot
#     deficits are overstated by **1.51x** against a Kurucz atmosphere (417 K rather than the
#     published 630 K) and by 2.31x against a blackbody. The honest figure is ~1.5x.
#
#     Their own 2016 methods section says they convert WITH Kurucz models ("for [Fe/H]=-0.25
#     and appropriate log g"), which the measured exponent of 0.96 does not match. That tension
#     is unresolved; this figure is the only evidence either way, and it says 0.96.
#
# The limb-darkening exponent their figure is drawn with is alpha = 0.285 best-fitting, with
# their fixed 0.269 statistically indistinguishable (r 0.99103 vs 0.99091) and 0.145 excluded
# (0.9555). So the figure is internally consistent with Table 6; the disagreement with the
# alpha the DATA prefer is real and is examined in `demos/zetand_fig10_family.jl`.

using ROTIR, Printf, Statistics, DelimitedFiles, Base.Threads

const PROT = 17.769426
const ZET  = joinpath(pkgdir(ROTIR), "demos", "zetAnd")
const NX   = 641
const PIX  = 0.5/94.9167                       # mas/px, Figure 10's own x ticks
const RX, RY, RZ = 2.6529/2, 2.4979/2, 2.4959/2
const INC  = 69.44
const LAM  = 1.605e-6
const PHI  = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875]

# ------------------------------------------------------------------------------------------
# rendering helpers
# ------------------------------------------------------------------------------------------
params(pa, ld1 = 0.269) =
    (surface_type = 1, radius_x = RX, radius_y = RY, radius_z = RZ, tpole = 4372.0,
     ldtype = 3, ld1 = ld1, ld2 = 0.0, beta = 0.25, inclination = INC,
     position_angle = pa, rotation_period = PROT)

"""Per-pixel facet index: the visible facet whose projected quad contains the pixel centre,
frontmost first. This is the answer to "which published pixels show which facet"."""
function facet_index(star, idx, nx, pix)
    fi = zeros(Int, nx, nx); best = fill(-Inf, nx, nx); c = nx ÷ 2 + 1
    @inbounds for q in idx
        xs = (star.proj_west[q,:] ./ pix) .+ c
        ys = (star.proj_north[q,:] ./ pix) .+ c
        i0 = max(1, floor(Int, minimum(ys))); i1 = min(nx, ceil(Int, maximum(ys)))
        j0 = max(1, floor(Int, minimum(xs))); j1 = min(nx, ceil(Int, maximum(xs)))
        nz = star.normals[q,3]
        for i in i0:i1, j in j0:j1
            nz <= best[i,j] && continue
            inside = true; px = Float64(j); py = Float64(i); s0 = 0.0
            for v in 1:4
                w = v == 4 ? 1 : v+1
                s = (xs[w]-xs[v])*(py-ys[v]) - (ys[w]-ys[v])*(px-xs[v])
                if v == 1; s0 = s; elseif s*s0 < 0; inside = false; break; end
            end
            inside && (fi[i,j] = q; best[i,j] = nz)
        end
    end
    fi
end

"""Rasterise the projected facet outlines: the colour-free registration signal."""
function edge_image(star, idx, nx, pix)
    img = zeros(nx, nx); c = nx ÷ 2 + 1
    @inbounds for q in idx, v in 1:4
        w = v == 4 ? 1 : v+1
        x1 = star.proj_west[q,v]/pix + c; y1 = star.proj_north[q,v]/pix + c
        x2 = star.proj_west[q,w]/pix + c; y2 = star.proj_north[q,w]/pix + c
        n = max(2, ceil(Int, 2*hypot(x2-x1, y2-y1)))
        for s in 0:n
            t = s/n
            ix = round(Int, x1 + t*(x2-x1)); iy = round(Int, y1 + t*(y2-y1))
            (1 <= ix <= nx && 1 <= iy <= nx) && (img[iy,ix] += 1.0)
        end
    end
    img
end

function gaussblur(A::Matrix{Float64}, s::Float64)
    r = ceil(Int, 3s); k = [exp(-0.5*(i/s)^2) for i in -r:r]; k ./= sum(k)
    ny, nx = size(A); B = similar(A); C = similar(A)
    @inbounds for j in 1:nx, i in 1:ny
        a = 0.0; for t in -r:r; a += k[t+r+1]*A[clamp(i+t,1,ny), j]; end; B[i,j] = a
    end
    @inbounds for j in 1:nx, i in 1:ny
        a = 0.0; for t in -r:r; a += k[t+r+1]*B[i, clamp(j+t,1,nx)]; end; C[i,j] = a
    end
    C
end

function erode(m::BitMatrix, n::Int)
    out = copy(m)
    for _ in 1:n
        p = copy(out)
        @inbounds for j in 2:size(m,2)-1, i in 2:size(m,1)-1
            out[i,j] = p[i,j] & p[i-1,j] & p[i+1,j] & p[i,j-1] & p[i,j+1]
        end
        out[1,:] .= false; out[end,:] .= false; out[:,1] .= false; out[:,end] .= false
    end
    out
end

"""The four parity variants of a HEALPix map, built by reflecting facet centres -- so this is
independent of any stored permutation file."""
function parity_perms(tess)
    c = tess.unit_xyz[:, 5, :]; n = size(c, 1)
    p(f) = [argmin(sum(abs2, c .- permutedims(f(c[i,:])), dims=2))[1] for i in 1:n]
    (:identity => collect(1:n), :latflip => p(v -> [v[1],  v[2], -v[3]]),
     :lonflip  => p(v -> [v[1], -v[2], v[3]]), :both => p(v -> [v[1], -v[2], -v[3]]))
end

flipimg(A, fw, fn) = (B = fw > 0 ? A : reverse(A, dims=2); fn > 0 ? B : reverse(B, dims=1))

function corr(a, b, m)
    i = findall(m); length(i) < 500 && return -2.0
    x = [a[k] for k in i]; y = [b[k] for k in i]
    (std(x) == 0 || std(y) == 0) ? -2.0 : cor(x, y)
end

readpanels() = begin
    rd(f, T) = permutedims(reshape(reinterpret(T, read(joinpath(ZET, f))), NX, NX, 8), (3,2,1))
    (rd("fig10_panels_V.bin", Float32), rd("fig10_panels_M.bin", UInt8),
     rd("fig10_panels_E.bin", Float32))
end

# ------------------------------------------------------------------------------------------
function pass1_edges(tess, Ep)
    println("PASS 1 — facet edges (no colour calibration, no map)")
    println("  panel     corr   render-phase    PA")
    phases = collect(range(0, 1 - 1/240, length = 240))
    best = [(-2.0, 0.0, 0.0, 0, 0) for _ in 1:8]; lk = ReentrantLock()
    for pa in (55.70, 124.30)
        @threads for ip in eachindex(phases)
            st  = create_star(tess, params(pa), phases[ip]*PROT)
            idx = findall(st.normals[:,3] .> 0)
            E   = gaussblur(edge_image(st, idx, NX, PIX), 2.0)
            fl  = rasterize_polygon_image(st.proj_west[idx,:], st.proj_north[idx,:],
                                          ones(length(idx)), PIX, NX)
            msk = erode(BitMatrix(fl .> 0.5*maximum(fl)), 6)
            for k in 1:8, fw in (1,-1), fn in (1,-1)
                c = corr(flipimg(E, fw, fn), Float64.(Ep[k,:,:]), flipimg(msk, fw, fn))
                lock(lk) do
                    c > best[k][1] && (best[k] = (c, phases[ip], pa, fw, fn))
                end
            end
        end
    end
    for k in 1:8
        @printf("   %.3f  %+.4f    %7.4f    %7.2f\n", PHI[k], best[k][1], best[k][2], best[k][3])
    end
    println("  => render_phase = -panel_phase (their phase runs the other way), " *
            "0.5 ambiguity left\n")
end

function pass2_values(tess, Vp, Mp, T10)
    println("PASS 2 — facet values (breaks the reflection and the 0.5 ambiguity)")
    P = parity_perms(tess); out = Vector{Any}(); lk = ReentrantLock()
    pas = [55.70, 124.30, 235.70, 304.30]
    @threads for ipa in eachindex(pas)
        pa = pas[ipa]
        for delta in (0.0, 0.5), sense in (-1, 1)
            FI = Vector{Matrix{Int}}(undef, 8); LD = Vector{Vector{Float64}}(undef, 8)
            NZ = Vector{Vector{Float64}}(undef, 8); IX = Vector{Vector{Int}}(undef, 8)
            for k in 1:8
                st = create_star(tess, params(pa), mod(sense*PHI[k] + delta, 1.0)*PROT)
                IX[k] = findall(st.normals[:,3] .> 0)
                FI[k] = facet_index(st, IX[k], NX, PIX)
                LD[k] = st.ldmap; NZ[k] = st.normals[:,3]
            end
            for (pn, perm) in P, fw in (1,-1), fn in (1,-1)
                X = Float64[]; Y = Float64[]
                for k in 1:8
                    fi = flipimg(FI[k], fw, fn); x = T10[perm]
                    for q in IX[k]
                        NZ[k][q] < 0.25 && continue
                        s = findall(i -> fi[i] == q && Mp[k,i] != 0x00, CartesianIndices(fi))
                        length(s) < 40 && continue
                        push!(Y, median([Float64(Vp[k,i]) for i in s]))
                        push!(X, x[q]*LD[k][q])
                    end
                end
                length(X) < 50 && continue
                b = cov(X,Y)/var(X); a = mean(Y) - b*mean(X)
                lock(lk) do
                    push!(out, (cor(X,Y), sqrt(mean((Y .- (a .+ b.*X)).^2)),
                                pa, delta, sense, pn, fw, fn, a, b))
                end
            end
        end
    end
    sort!(out, by = r -> -r[1])
    println("      r      rms     PA    delta sense  parity    col~W row~N")
    for r in out[1:6]
        @printf("   %+.5f %6.2f  %6.2f   %.1f   %+d   %-9s %+d   %+d\n",
                r[1], r[2], r[3], r[4], r[5], string(r[6]), r[7], r[8])
    end
    @printf("\n  winner: published = %.1f + %.4f * (T * mu^0.269)\n", out[1][9], out[1][10])
    out[1]
end

function pass3_laws(tess, Vp, Mp, T10)
    println("\nWHAT THEY PLOTTED — per-facet samples at the winning configuration")
    PV = Float64[]; MU = Float64[]; TT = Float64[]; PK = Int[]
    for k in 1:8
        st = create_star(tess, params(124.30), mod(PHI[k] + 0.5, 1.0)*PROT)
        idx = findall(st.normals[:,3] .> 0)
        # the winning configuration is (PA 124.30, delta 0.5, sense +1, identity map) with
        # col ~ +West and row ~ -North, so the render's rows reverse into the published frame.
        fi = flipimg(facet_index(st, idx, NX, PIX), 1, -1)
        for q in idx
            st.normals[q,3] < 0.25 && continue
            s = findall(i -> fi[i] == q && Mp[k,i] != 0x00, CartesianIndices(fi))
            length(s) < 40 && continue
            push!(PV, median([Float64(Vp[k,i]) for i in s]))
            push!(MU, st.normals[q,3]); push!(TT, T10[q]); push!(PK, k)
        end
    end
    score(X) = (b = cov(X,PV)/var(X); a = mean(PV) - b*mean(X);
                (cor(X,PV), sqrt(mean((PV .- (a .+ b.*X)).^2))))
    # Did they auto-scale each panel, or share one normalisation? Fit the affine relation panel
    # by panel: a spread in A and B would mean each imshow set its own vmin/vmax, and then no
    # single calibration could be quoted for the figure.
    println("  colour normalisation, fitted separately per panel:")
    for k in 1:8
        sel = PK .== k
        X = TT[sel] .* MU[sel].^0.269; Y = PV[sel]
        b = cov(X,Y)/var(X); a = mean(Y) - b*mean(X)
        @printf("    phi %.3f   A = %7.1f   B = %.4f\n", PHI[k], a, b)
    end
    println("  intensity law at their fixed alpha = 0.269:")
    for (nm, I) in (("I = T (linear)", TT),
                    ("I = Planck(T, 1.605um)", [ROTIR.planck_and_dT(t, LAM)[1] for t in TT]),
                    ("I = T^4", TT.^4))
        r, rms = score(I .* MU.^0.269)
        @printf("    %-24s r = %+.5f   rms = %6.2f K\n", nm, r, rms)
    end
    best = (-2.0, 0.0)
    for al in 0.0:0.005:0.60
        r, _ = score(TT .* MU.^al); r > best[1] && (best = (r, al))
    end
    println("  limb-darkening exponent, I = T * mu^alpha:")
    for al in (0.145, 0.269, best[2], 0.400)
        r, rms = score(TT .* MU.^al)
        @printf("    alpha = %.3f%s   r = %+.5f   rms = %6.2f K\n",
                al, al == best[2] ? " (best)" : "      ", r, rms)
    end
end

function dump_model(tess, T10)
    A, B = 3495.8, 0.2478              # the global per-facet regression of pass 2
    mdl = zeros(Float32, NX, NX, 8); muim = zeros(Float32, NX, NX, 8)
    println("\nsub-observer geometry at the winning configuration")
    println("  panel   sub-observer longitude   visible facets")
    for k in 1:8
        st = create_star(tess, params(124.30), mod(PHI[k] + 0.5, 1.0)*PROT)
        idx = findall(st.normals[:,3] .> 0); fi = facet_index(st, idx, NX, PIX)
        for i in CartesianIndices(fi)
            q = fi[i]; q == 0 && continue
            mdl[i,k] = A + B*T10[q]*st.ldmap[q]; muim[i,k] = st.normals[q,3]
        end
        qb = idx[argmax(st.normals[idx,3])]; c = tess.unit_xyz[qb,5,:]
        @printf("  %.3f          %6.1f deg               %4d\n",
                PHI[k], mod(atand(c[2], c[1]), 360), length(idx))
    end
    write(joinpath(ZET, "fig10_model_panels.bin"), mdl)
    write(joinpath(ZET, "fig10_model_mu.bin"), muim)
    println("  => sub-observer longitude = 90 deg - 360 deg * phi")
    println("\nwrote fig10_model_panels.bin, fig10_model_mu.bin")
    println("next: python3 demos/zetand_fig10_panels.py --figure")
end

function main()
    isfile(joinpath(ZET, "fig10_panels_V.bin")) ||
        error("run `python3 demos/zetand_fig10_panels.py` first")
    Vp, Mp, Ep = readpanels()
    T10  = vec(readdlm(joinpath(ZET, "zetAnd_fig10_epochA_T.txt")))
    tess = tessellation_healpix(3; T = Float64)
    pass1_edges(tess, Ep)
    pass2_values(tess, Vp, Mp, T10)
    pass3_laws(tess, Vp, Mp, T10)
    dump_model(tess, T10)
end

main()
