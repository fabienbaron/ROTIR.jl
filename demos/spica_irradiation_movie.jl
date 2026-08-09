# spica_irradiation_movie.jl
#
# A movie of the two Spica components through one orbit: tidally distorted (Roche)
# surfaces whose shape breathes with D(t), von Zeipel gravity darkening, and mutual
# irradiation on top.
#
# Run:  julia --project=demos demos/spica_irradiation_movie.jl
#
# Knobs (environment variables):
#   NFRAMES=200   frames over one orbital period
#   NSIDE1=4      HEALPix level of the primary   (12·4^n tessels: 4 → 3072)
#   NSIDE2=3      HEALPix level of the secondary (3 → 768)
#   ALBEDO=0.6    bolometric albedo
#   OUTDIR=...    where the frames and mp4 go
#   VCONS=1       volume-conserving Roche shape (0 for the fixed-rpole convention)
#
# Output (in OUTDIR):
#   single/  — the irradiated system alone
#   compare/ — irradiated | intrinsic only | ΔT, which is the pair that actually shows
#              the effect: on its own a heated star just looks like a slightly different
#              star.
#   substellar_dT.png — peak ΔT on each component vs orbital phase.

ENV["MPLBACKEND"] = get(ENV, "MPLBACKEND", "Agg")   # headless by default
using ROTIR, PyPlot, Printf
include(joinpath(@__DIR__, "spica_params.jl"))

nframes = parse(Int,     get(ENV, "NFRAMES", "200"))
nside1  = parse(Int,     get(ENV, "NSIDE1",  "4"))
nside2  = parse(Int,     get(ENV, "NSIDE2",  "3"))
albedo  = parse(Float64, get(ENV, "ALBEDO",  string(ALBEDO_DEFAULT)))
vcons   = get(ENV, "VCONS", "1") == "1"
outdir  = get(ENV, "OUTDIR", joinpath(@__DIR__, "results", "spica_movie"))

spica_audit()

tes1 = tessellation_healpix(nside1)
tes2 = tessellation_healpix(nside2)
@info "tessellation: primary $(tes1.npix) tessels, secondary $(tes2.npix) tessels; " *
      "albedo = $albedo, volume-conserving = $vcons, $nframes frames"

t0, t1 = T0_ORB, T0_ORB + P_ORB      # one full orbit, starting at periastron

# ---------------------------------------------------------------------------------------
# 1. Movies
# ---------------------------------------------------------------------------------------
# `band` is the H-band effective wavelength of the CHARA/MIRC data used in the companion
# detectability demo; with intensity_model = :planck the two components' relative
# brightness is a true Planck ratio rather than the T-proportional approximation.
for panels in (:single, :compare)
    dir = joinpath(outdir, String(panels))
    @info "rendering $(panels) movie into $dir"
    _, mp4 = binary_movie(SPICA_BP, tes1, SPICA_P1, tes2, SPICA_P2;
                          tstart=t0, tstop=t1, nframes=nframes, panels=panels,
                          outdir=dir, fps=20, dpi=110,
                          reflection=true, albedo1=albedo, albedo2=albedo,
                          method=:horvat, volume_conserving=vcons,
                          intensity=true, intensity_model=:planck, band=1.609e-6,
                          verbose=true)
    mp4 === nothing || @info "wrote $mp4"
end

# ---------------------------------------------------------------------------------------
# 2. Peak ΔT vs orbital phase
# ---------------------------------------------------------------------------------------
# Quantifies what the ΔT panel of the movie shows. The curve should peak at periastron
# (phase 0) and be a factor ((1+e)/(1-e))² ≈ 1.6 lower at apastron.
nph = 121
phases = range(0, 1, length=nph)
dT1 = zeros(nph); dT2 = zeros(nph); seps = zeros(nph)
Ω1tab = vcons ? roche_omega_table(tes1, SPICA_P1, SPICA_BP; secondary=false) : nothing
Ω2tab = vcons ? roche_omega_table(tes2, SPICA_P2, SPICA_BP; secondary=true)  : nothing
for (k, ph) in enumerate(phases)
    t = T0_ORB + ph * P_ORB
    D = compute_separation(SPICA_BP, t)
    o1 = Ω1tab === nothing ? nothing : Ω1tab(D)
    o2 = Ω2tab === nothing ? nothing : Ω2tab(D)
    s1, s2, h1, h2, i1, i2 = binary_frame_maps(tes1, SPICA_P1, tes2, SPICA_P2, SPICA_BP, t;
                                               albedo1=albedo, albedo2=albedo,
                                               omega1=o1, omega2=o2)
    dT1[k] = maximum(h1 .- i1); dT2[k] = maximum(h2 .- i2)
    seps[k] = binary_frame(SPICA_BP, t)[2]
end

fig, ax = subplots(2, 1, figsize=(9, 7), sharex=true)
ax[1].plot(phases, dT2, "-", lw=2, label="secondary")
ax[1].plot(phases, dT1, "-", lw=2, label="primary")
ax[1].set_ylabel("peak ΔT (K)")
ax[1].legend(); ax[1].grid(alpha=0.3)
ax[1].set_title(@sprintf("Spica mutual irradiation, bolometric albedo = %.2f", albedo))
ax[2].plot(phases, seps, "k-", lw=2)
ax[2].set_ylabel("separation D (mas)")
ax[2].set_xlabel("orbital phase (0 = periastron)")
ax[2].grid(alpha=0.3)
tight_layout()
mkpath(outdir)
fig.savefig(joinpath(outdir, "substellar_dT.png"), dpi=130, bbox_inches="tight")
PyPlot.close(fig)

@printf("\npeak ΔT: secondary %.0f K (periastron) … %.0f K (apastron)\n",
        dT2[1], dT2[div(nph, 2) + 1])
@printf("peak ΔT: primary   %.0f K (periastron) … %.0f K (apastron)\n",
        dT1[1], dT1[div(nph, 2) + 1])
@info "done — output in $outdir"
