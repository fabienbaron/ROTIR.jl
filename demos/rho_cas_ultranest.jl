#!/usr/bin/env julia
# ρ Cas posterior by nested sampling (UltraNest, via PythonCall).
#
#   julia --project=demos -t 1 demos/rho_cas_ultranest.jl
#   FREE=rpole,omega,inc,PA,beta,ld1 STEPSAMPLER=1 NLIVE=400 julia --project=demos -t 1 demos/rho_cas_ultranest.jl
#
# NOTE the `-t 1`. UltraNest is reached through PythonCall: `pydecref_`
# calls Py_DecRef from a GC finalizer with no GIL check, so a PyObject finalized on a worker
# thread segfaults the process mid-run. `fit_parametric_ultranest` now refuses to start with
# more than one thread rather than crashing hours in. Note JULIA_NUM_THREADS in the
# environment counts — `-t 1` on the command line is not enough if it is set.
#
# Needs the Python package:  ~/.julia/conda/3/x86_64/bin/pip install ultranest
#
# Nested sampling needs no starting point and no gradients, explores every mode of a
# multimodal posterior, and returns log(Z) with an uncertainty — so it is the method that
# can actually answer "is the oblate solution preferred over the sphere?" rather than just
# fitting one of them. The price is evaluations: it ignores the gradients ROTIR can supply
# and its cost grows steeply with dimension.
#
# Rules of thumb at n=3/Float64 on ρ Cas (~0.13 s per likelihood):
#   2 parameters, region sampler   ~10-20k evaluations   → under an hour
#   7 parameters, step sampler     ~1e5-1e6 evaluations  → hours to a day
# Above ~5 parameters set STEPSAMPLER=1; the region sampler degrades badly.

using Printf, Statistics
# UltraNest lives in ROTIRUltraNestExt, which PythonCall triggers: ROTIR does not load
# Python on its own any more.
using PythonCall
include(joinpath(@__DIR__, "rho_cas_model.jl"))
include(joinpath(@__DIR__, "posterior_utils.jl"))

const NLIVE       = parse(Int, get(ENV, "NLIVE", "400"))
const STEPSAMPLER = get(ENV, "STEPSAMPLER", NDIM > 4 ? "1" : "0") == "1"
const NSTEPS      = parse(Int, get(ENV, "NSTEPS", string(4 * NDIM)))
const FRAC_REMAIN = parse(Float64, get(ENV, "FRAC_REMAIN", "0.01"))

describe_model()
@printf("UltraNest: %d live points, stepsampler=%s%s\n", NLIVE, STEPSAMPLER,
        STEPSAMPLER ? " (nsteps=$NSTEPS)" : "")

# fit_parametric_ultranest lives in ROTIR core — nested sampling needs no AD backend, and
# PythonCall is already a ROTIR dependency. The box below is the same uniform prior the
# gradient-based samplers get through their logit transform, so the posteriors match.
res, wall = timed(() -> fit_parametric_ultranest(DATA, TESSELS, TEPOCHS, BASE;
    θ0 = THETA0, free = FREE_NAMES, lb = BOX_LO_FULL, ub = BOX_HI_FULL,
    min_num_live_points = NLIVE, frac_remain = FRAC_REMAIN,
    use_stepsampler = STEPSAMPLER, nsteps = NSTEPS,
    log_dir = joinpath(RESULTS_DIR, "ultranest_logs"), resume = "overwrite");
    label = "ultranest")

S = res.samples
@printf("\nwall time %.1f s   log(Z) = %.3f ± %.3f   %d samples\n",
        wall, res.logz, res.logzerr, size(S, 1))
summarise(S, LABELS; title = "UltraNest posterior:")
save_posterior("ultranest", S, LABELS; wall_seconds = wall, logz = res.logz,
               logzerr = res.logzerr, nlive = NLIVE, stepsampler = STEPSAMPLER)
corner_plot(S, LABELS, "ultranest_corner.png"; title = "ρ Cas — UltraNest")

println("\nlog(Z) is the quantity to compare across models (sphere vs oblate vs oblate +")
println("gravity darkening): rerun with FREE=... and difference the log evidences. Under a")
println("misspecified model (χ²ᵣ ≈ 6 here) a Bayes factor mostly ranks which model absorbs")
println("the systematics best, so read it as a diagnostic, not a verdict.")
