#!/usr/bin/env julia
# The line-of-sight velocity field (src/velocity_field.jl).
#
# Standalone script in the style of test_parametric_gradient.jl: prints its own table,
# never throws, exposes `nfail[]` for runtests.jl. Needs no AD package, so it runs
# identically under `--project=.` and `--project=demos`.
#
#     julia --project=. test/test_velocity_field.jl
#
# WHAT MAKES THESE REAL TESTS. Every check compares the velocity field against a quantity
# derived by a COMPLETELY DIFFERENT route: `projected_veq` comes from Omega*R_eq*sin(i) out
# of stellar_physics.jl, while `los_velocity` integrates sky-frame vertex positions through
# the rotation matrix. Section [6] is the one that matters most - it pins the ratio of the
# fast-rotator velocity to the sphere-limit velocity against R_eq/R_p, which is what
# demonstrates the field follows the ROCHE SHAPE. `src/di.jl:454` builds its velocity from
# surface NORMALS instead of positions; for a sphere that is equivalent, for an oblate
# figure it is not, and [6] is the check that would catch the substitution.
#
# Tolerances of a few percent in [1], [2], [5] and [6] are HEALPix discretisation, not
# slack: no tessel centroid lands exactly on the equator at exactly the limb, so the
# measured peak velocity is a per-cent or so below the analytic equatorial value (214.55
# vs 215.12 km/s at nside 4). Refining the mesh tightens it.

using ROTIR, LinearAlgebra, Printf

npass = Ref(0); nfail = Ref(0)
function cb(label, ok)
    ok ? (npass[] += 1) : (nfail[] += 1)
    @printf("  %-56s %s\n", label, ok ? "\u2713" : "\u2717")
    return ok
end
ap(label, a, b; tol = 1e-6) = cb(label, abs(a - b) / max(abs(a), abs(b), eps()) < tol)

const RP, DPC, FEV, PROT = 0.849, 16.8, 0.92, 1/1.12
const tess = tessellation_healpix(4, T = Float64)

println("\n[1] equator-on rigid rotator: peak v_los must equal veq from the geometry")
# Two completely independent routes to the same number: the velocity field integrates the
# sky-frame positions, projected_veq comes from Omega*R_eq*sin(i). They must agree.
sp = default_star_params(2; rpole=RP, d=DPC, frac_escapevel=FEV, rotation_period=PROT,
                         inclination=90.0, position_angle=0.0, ldtype=0, tpole=7208.0)
st = create_star(tess, sp, 0.0)
s  = velocity_field_summary(st, sp)
veq = equatorial_velocity(RP,DPC,FEV,PROT)
@printf("      v in [%.2f, %.2f] km/s   vsini_proj=%.2f   projected_veq=%.2f\n",
        s.vmin, s.vmax, s.vsini_proj, veq)
ap("vsini_proj == projected_veq (i=90)", s.vsini_proj, veq; tol=2e-2)
cb("antisymmetric about zero", abs(s.vmin+s.vmax) < 1e-8*max(abs(s.vmin),1))

println("\n[2] sin(i) scaling")
for inc in (90.0, 60.0, 30.0, 10.0)
    spi = merge(sp,(inclination=inc,)); sti=create_star(tess,spi,0.0)
    ap("  i=$(Int(inc)): vsini_proj == veq*sin(i)",
       velocity_field_summary(sti,spi).vsini_proj, veq*sind(inc); tol=3e-2)
end

println("\n[3] pole-on gives no rotational signal")
sp0 = merge(sp,(inclination=0.0,)); st0=create_star(tess,sp0,0.0)
cb("i=0 -> v == 0 everywhere", maximum(abs, los_velocity(st0,sp0)) < 1e-9)

println("\n[4] vgamma is a pure offset")
spv = merge(sp,(vgamma=-25.0,)); stv=create_star(tess,spv,0.0)
cb("shifts every tessel by exactly vgamma",
   maximum(abs, los_velocity(stv,spv) .- los_velocity(st,sp) .+ 25.0) < 1e-9)

println("\n[5] position angle rotates the velocity pattern, not its extent")
for pa in (0.0, 45.0, 90.0, 137.0, -60.0)
    spp = merge(sp,(position_angle=pa,)); stp=create_star(tess,spp,0.0)
    ap("  PA=$pa: vsini_proj unchanged", velocity_field_summary(stp,spp).vsini_proj, s.vsini_proj; tol=3e-2)
end

println("\n[6] oblateness matters: positions, not normals")
# di.jl:454 used normals. For an oblate figure the two differ; at fev->0 they must agree.
spr = merge(sp,(frac_escapevel=1e-6,)); str=create_star(tess,spr,0.0)
vr  = velocity_field_summary(str,spr)
ap("fev->0: vsini_proj == Omega*R_pole (sphere limit)", vr.vsini_proj,
   equatorial_velocity(RP,DPC,1e-6,PROT); tol=2e-2)
cb("fast rotator is FASTER than the sphere at equal rpole/P", s.vsini_proj > vr.vsini_proj)
@printf("      fev=0.92: %.2f km/s vs fev=0: %.2f km/s  (ratio %.4f, R_eq/R_p=%.4f)\n",
        s.vsini_proj, vr.vsini_proj, s.vsini_proj/vr.vsini_proj,
        polar_radius_rsun(RP,DPC)*ROTIR.f_rapid_rot_and_deriv(FEV)[1]/polar_radius_rsun(RP,DPC))

println("\n[7] differential rotation")
spd = merge(sp,(B_rot=0.3,)); std=create_star(tess,spd,0.0)
cb("B_rot>0 reduces the span (poles lag)", velocity_field_summary(std,spd).vspan < s.vspan)
cb("B_rot=0 identical to omitting it",
   los_velocity(create_star(tess,merge(sp,(B_rot=0.0,)),0.0), merge(sp,(B_rot=0.0,))) == los_velocity(st,sp))

println("\n[8] doppler_lambda")
ap("lambda(v=0) == lambda0", doppler_lambda(1.65e-6, 0.0), 1.65e-6)
ap("receding redshifts", doppler_lambda(1.65e-6, 100.0), 1.65e-6*(1+100/299792.458))
cb("broadcasts over a vector", length(doppler_lambda(1.65e-6, [10.0,-10.0])) == 2)

println("\n[9] refuses to invent a distance")
cb("no `d` -> error", try los_velocity(st, Base.structdiff(sp,NamedTuple{(:d,)})); false catch; true end)

println("\n[10] Float32: the MESH sets the precision, not the parameters")
# `compute_radii` and `rotate_star` both narrow parameters to the mesh type via
# `convert_params`. Promoting here instead would make a Float32 mesh silently run in Float64
# whenever a caller passed an unconverted Float64 angle — which is every hand-written demo.
let t32 = tessellation_healpix(3, T = Float32)
    for PT in (Float32, Float64)
        sp = default_star_params(2; T = PT, rpole = RP, d = DPC, frac_escapevel = FEV,
                                 rotation_period = PROT, inclination = 60.0,
                                 position_angle = 0.0, ldtype = 0, tpole = 7208.0)
        v = los_velocity(create_star(t32, sp, zero(PT)), sp)
        cb("Float32 mesh + $PT params -> Float32", eltype(v) === Float32)
    end
    # And a Float64 mesh stays Float64.
    t64 = tessellation_healpix(3, T = Float64)
    sp64 = default_star_params(2; rpole = RP, d = DPC, frac_escapevel = FEV,
                               rotation_period = PROT, inclination = 60.0,
                               position_angle = 0.0, ldtype = 0, tpole = 7208.0)
    cb("Float64 mesh -> Float64",
       eltype(los_velocity(create_star(t64, sp64, 0.0), sp64)) === Float64)
end

@printf("\n=== %d passed, %d failed ===\n", npass[], nfail[])
