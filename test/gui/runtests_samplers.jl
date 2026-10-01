#!/usr/bin/env julia
# The GUI's SAMPLER paths, which test/gui/runtests.jl structurally cannot reach.
#
#     xvfb-run -a julia --project=bin test/gui/runtests_samplers.jl
#
# WHY A SECOND FILE, AND NOT A TESTSET IN THE FIRST.
#
# `runtests.jl`'s opening testset ASSERTS that Pigeons and PythonCall are not loaded, and that
# is not pedantry: with Pigeons loaded, one canvas build goes from 370 ms to 3166 ms, and with
# PythonCall from 338 ms to 2477 ms, on every GUI start. The launcher therefore does not load
# them, `hmc_available()` / `nautilus_available()` / `pigeons_available()` are all FALSE in that
# process, and `shell_fit("hmc", …)` there exercises a refusal rather than NUTS. Four of the six
# engines in `FIT_METHODS` were untested for exactly this reason — the suite could not load what
# it needed without breaking its own startup-cost guarantee.
#
# So the samplers get their own process, which loads them deliberately and pays the cost once.
#
# WHAT THIS IS AND IS NOT. A SMOKE test: that the GUI's job machinery drives each sampler to
# completion, records a fit with a posterior, and reports it. The budgets are deliberately tiny
# and NONE of these results is a measurement — a 20-draw NUTS chain on one night of λ And is not
# converged and is not meant to be. What would break without this file is a sampler wired up
# wrongly in the GUI (dropped kwargs, a θ layout mismatch, a posterior that never reaches
# `FitEntry.samples`), and those are failures the budget does not have to be large to catch.
#
# `bin/Project.toml` already has AdvancedHMC, Nautilus, Pigeons and Zygote, so nothing new is
# needed to run it.

using Test
using ROTIR
using GLMakie, QMLMakie, QML
using Makie
# THE POINT OF THIS FILE: loading these is what populates `_fit_hmc`, `_fit_nautilus` and
# `_fit_pigeons`, since availability is `!isempty(methods(...))` on a weakdep extension.
using Zygote, AdvancedHMC, LogDensityProblems, ADTypes, Distributions
using Nautilus
using Pigeons

const G = Base.get_extension(ROTIR, :ROTIRGUIExt)
G === nothing && error("ROTIRGUIExt did not load; GLMakie, QMLMakie and QML are all needed")

const DATA = joinpath(pkgdir(ROTIR), "demos", "data")
const LAM  = sort([joinpath(DATA, f) for f in readdir(DATA) if occursin("lam_And", f)])

"A ShellState with every canvas built, as `gui()` builds them. Kept in step with runtests.jl."
function fresh_shell()
    s = G.Session()
    figs = [Makie.Figure() for _ in 1:10]
    sh = G.ShellState(; session = s,
        sky    = G.build_sky_canvas(figs[1]),  star  = G.build_star_canvas(figs[2]),
        moll   = G.build_moll_canvas(figs[3]), chi2  = G.build_chi2_canvas(figs[4]),
        imsky  = G.build_sky_canvas(figs[5]),  immoll = G.build_moll_canvas(figs[6]),
        msky   = G.build_sky_canvas(figs[7]),
        obsmodel = Makie.Observable(Makie.Point2f[]),
        orbitcanvas = G.build_orbit_canvas(figs[8]),
        post   = G.build_post_canvas(figs[9]),
        imstar = G.build_star_canvas(figs[10]))
    G.SHELL[] = sh
    return sh
end

"Run a job callback to completion, the way the QML poll timer does."
function drain!(; limit = 3000)
    for _ in 1:limit
        G.shell_job_running() == "1" || break
        sleep(0.2); G.shell_job_poll()
    end
    G.shell_job_poll()
    return nothing
end

rows(s) = isempty(s) ? String[] : split(s, '\n')
cols(r) = split(r, '\t')

"The smallest model a sampler can be pointed at: a sphere with one free radius."
function one_free_sphere!()
    sh = fresh_shell()
    G.shell_open(LAM[1], "0")
    G.shell_set_tessellation("healpix", 2, "Float32")   # 192 tessels: the cheapest real mesh
    G.shell_add_model(0)
    G.shell_set_param("radius", "1.2")
    G.shell_set_param_state("radius", "free")
    G.shell_set_bound("radius", "0.8", "2.0")
    return sh
end

@testset "ROTIR GUI samplers" begin

@testset "loading the samplers is what makes them available" begin
    fresh_shell()          # `shell_fit` goes through `_sh()`, which asserts a shell exists
    # The inverse of runtests.jl's opening assertion, and the reason this file exists. If any of
    # these is false the rest of the file is testing refusals, which is what it is here to stop.
    @test hmc_available()
    @test nautilus_available()
    @test pigeons_available()
    # UltraNest stays absent even here: it needs PythonCall, and the GUI is Python-free BY
    # CONSTRUCTION — `:ultranest` is not in `FIT_METHODS` at all, so this is a design assertion,
    # not an availability one.
    @test !occursin("ultranest", G.shell_fit_methods())
    @test occursin("does not offer :ultranest", G.shell_fit("ultranest", 100))
end

@testset "every engine the combo offers is now offered for real" begin
    sh = one_free_sphere!()
    # With the samplers loaded, the method list must GROW — `shell_fit_methods` filters on
    # availability, and a sampler present but unlisted is unreachable from the window.
    listed = [cols(r)[1] for r in rows(G.shell_fit_methods())]
    for k in ("gradient", "neldermead", "bobyqa", "hmc", "nautilus", "pigeons")
        @test k in listed
    end
    # Each one carries a budget with a unit, which is what the Fit panel labels its spinbox
    # with — a sampler measured in "iterations" would mislead about what the number buys.
    for (k, unit) in (("hmc", "draws"), ("nautilus", "eff. samples"), ("pigeons", "rounds"))
        b = cols(G.shell_fit_budget(k))
        @test b[1] == unit && parse(Int, b[2]) > 0
    end
end

@testset "NUTS through the GUI leaves a posterior behind" begin
    sh = one_free_sphere!()
    @test G.shell_fit("hmc", 20) == ""
    drain!()
    @test G.shell_job_running() == "0"
    f = G.current_fit(sh.session)
    @test f !== nothing
    @test f.method === :hmc
    @test f.names == [:radius]
    # THE ASSERTION THAT MATTERS: a sampler must deposit SAMPLES, not just a point. This is what
    # distinguishes a working sampler path from one whose draws are thrown away on the way back
    # into the session.
    @test size(f.samples, 1) > 0
    @test size(f.samples, 2) == 1
    @test G.shell_fit_has_posterior() == "1"
    # A point estimate AND an interval, with the interval ordered — `best` is the posterior
    # median for a sampler, and `err` half the 16–84 range rather than NaN as for an optimiser.
    @test 0.8 <= f.best[1] <= 2.0
    @test f.q16[1] <= f.best[1] <= f.q84[1]
    @test isfinite(f.err[1]) && f.err[1] >= 0
    @test all(0.8 .<= f.samples[:, 1] .<= 2.0)       # the box reparameterisation held
    # The DIAGNOSTIC string, which is what the fits table shows and the only place a user sees
    # how many draws diverged. A sampler whose draws arrive but whose diagnostics do not is
    # indistinguishable in the panel from a converged one.
    @test occursin("draws", f.diagnostics) && occursin("divergences", f.diagnostics)
    # And the posterior panel has a pair to plot, which is a separate path from storing it.
    @test G.shell_posterior_pair() != ""
    @test occursin("radius", G.shell_fit_params())
end

@testset "Nautilus reports an evidence, which is why it is offered" begin
    sh = one_free_sphere!()
    @test G.shell_fit("nautilus", 300) == ""
    drain!()
    f = G.current_fit(sh.session)
    @test f !== nothing && f.method === :nautilus
    @test size(f.samples, 1) > 0
    @test G.shell_fit_has_posterior() == "1"
    # THE DISTINGUISHING OUTPUT. Nelder–Mead and VMLMB give a χ²; a nested sampler gives
    # log Z with an error, and that is the only thing in this GUI that can compare two models
    # with different numbers of parameters. A NaN here means the field is plumbed but unfilled.
    @test isfinite(f.logz)
    @test isfinite(f.logzerr) && f.logzerr >= 0
    # Columns 6 and 7 of the fits table are log Z and its error, and they print "—" when the
    # field is not finite — so a dash here is the visible symptom of an unfilled evidence.
    nrow = cols(rows(G.shell_fits())[end])
    @test nrow[6] != "—" && nrow[7] != "—"
end

@testset "Pigeons runs, and its round trips are reported" begin
    sh = one_free_sphere!()
    # `rounds`, not draws: Pigeons' budget is 2^n scans, so a small number is still a real run.
    @test G.shell_fit("pigeons", 3) == ""
    drain!()
    f = G.current_fit(sh.session)
    @test f !== nothing && f.method === :pigeons
    @test size(f.samples, 1) > 0
    @test G.shell_fit_has_posterior() == "1"
    @test all(0.8 .<= f.samples[:, 1] .<= 2.0)
    # THE DIAGNOSTIC THAT DECIDES WHETHER A PIGEONS RESULT MAY BE USED AT ALL. A real β Cas run
    # in this repo finished with `round_trips = 0` and a global barrier Λ ≈ 7.95 against 10
    # chains, and its posterior was therefore discarded — the χ² at the median was WORSE than a
    # local optimiser's in the same basin. If the GUI does not surface the barrier and the round
    # trips, a user cannot tell that case from a converged one, so their presence on the console
    # is part of the contract and not decoration.
    @test occursin("round trips", f.diagnostics)
    @test occursin("chains", f.diagnostics)
    @test occursin("pigeons", lowercase(G.shell_console()))
    # The table surfaces it too — the diagnostics are its last column.
    @test occursin("round trips", cols(rows(G.shell_fits())[end])[end])
end

@testset "a sampler and an optimiser can be compared, which is the point of keeping both" begin
    # Four fits on ONE model, then the comparison table. This is the integration the separate
    # testsets above cannot make: the fit library has to hold entries from different engines
    # side by side with their χ² on the same data and the same free set.
    sh = one_free_sphere!()
    for (meth, budget) in (("neldermead", 150), ("gradient", 30), ("hmc", 15), ("nautilus", 200))
        G.shell_fit(meth, budget); drain!()
    end
    @test length(sh.session.fits) == 4
    @test [f.method for f in sh.session.fits] == [:neldermead, :gradient, :hmc, :nautilus]
    # Every entry fitted the SAME parameter on the SAME data, so χ² and ndata are comparable —
    # which is exactly the mistake the β Cas audit made across models with different weight
    # vectors (1260 vs 1980 data points), and the field that makes it checkable.
    @test all(f.names == [:radius] for f in sh.session.fits)
    @test length(unique(f.ndata for f in sh.session.fits)) == 1
    @test all(isfinite(f.chi2) && f.chi2 > 0 for f in sh.session.fits)
    # The two samplers carry a posterior; the two optimisers do not, and `err` is NaN for them
    # rather than zero — a zero error would read as an infinitely precise measurement.
    bykind = Dict(f.method => f for f in sh.session.fits)
    @test size(bykind[:hmc].samples, 1) > 0 && size(bykind[:nautilus].samples, 1) > 0
    @test size(bykind[:neldermead].samples, 1) == 0 && size(bykind[:gradient].samples, 1) == 0
    @test all(isnan, bykind[:neldermead].err) && all(isnan, bykind[:gradient].err)
    # All four land in the same neighbourhood: on a one-parameter sphere there is no degeneracy
    # to disagree about, so a wide spread here means an engine is wired to a different model.
    best = [f.best[1] for f in sh.session.fits]
    @test maximum(best) - minimum(best) < 0.3
    # And the library lists all four with their method names, which is what the panel shows.
    tbl = G.shell_fits()
    for s in ("Nelder", "VMLMB", "NUTS", "Nautilus")
        @test occursin(s, tbl)
    end
end

@testset "a sampler is refused rather than attempted where it cannot work" begin
    # NUTS and Pigeons need an analytic gradient, so they are offered only for a model that HAS
    # one. A surface type without a gradient path must be filtered out of the list, not left to
    # fail inside the sampler after the user has waited for it.
    sh = fresh_shell()
    G.shell_open(LAM[1], "0")
    # No model at all: the gradient-requiring engines cannot be offered yet.
    listed = [cols(r)[1] for r in rows(G.shell_fit_methods())]
    @test !("gradient" in listed)
    @test !("hmc" in listed) && !("pigeons" in listed)
    # Nothing free is reported by every engine the same way, before any work is done.
    G.shell_add_model(0)
    for meth in ("hmc", "nautilus", "pigeons")
        @test occursin("nothing is free", G.shell_fit(meth, 10))
    end
    @test G.shell_job_running() == "0"
end

end
