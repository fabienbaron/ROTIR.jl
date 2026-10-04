#!/usr/bin/env julia
# A NEGATIVE RESULT, kept so it is not repeated.
#
#     julia --project=. -t 16 demos/zetand_ensemble_test.jl
#
# SURFING imposes no spatial prior — the paper describes only MCMC, and its published map is
# "an average of 50 instances (24 walkers over 500 iterations)". So its smoothness is emergent,
# from averaging, rather than imposed. The idea tested here was that an ensemble of
# positivity-only reconstructions from randomised starts, averaged, would reproduce the same
# smoothing and remove the regularisation weight as a free choice.
#
# IT DOES NOT. Measured over 50 members: chi2/n 0.989 +- 0.000, correlation with the published
# map 0.066 +- 0.002, contrast 0.00-8.21 — and the ensemble MEAN is identical to any single
# member (chi2/n 0.988, r 0.0660, contrast 0.000-8.212). Every randomised start converges to the
# same map, so there is nothing to average.
#
# The reason is the distinction the experiment was built to test: averaging over OPTIMISER
# STARTS is not averaging over POSTERIOR SAMPLES. With 14218 data against 741 visible tessels
# the minimum is unique — very ill-conditioned, but unique — so a deterministic optimiser finds
# it from anywhere. SURFING's smoothing comes from sampling the posterior's UNCERTAINTY in the
# ~690 data-unconstrained modes, which no amount of re-optimisation explores. Reproducing it
# needs an actual sampler, not an ensemble of fits.
using ROTIR, Printf, DelimitedFiles, Statistics, Random
const PROT=17.769426; const PHASE0=0.511
const ZET=joinpath(pkgdir(ROTIR),"demos","zetAnd")
const FIGURE=(surface_type=1, radius_x=1.3275, radius_y=1.2663, radius_z=1.2120,
              tpole=5000.0, ldtype=3, ld1=0.139, ld2=0.0, beta=0.25,
              inclination=67.8, position_angle=126.2, rotation_period=PROT)
function main()
    files=sort(filter(f->endswith(f,".oifits"),readdir(ZET)))
    data=[readoifits(joinpath(ZET,f);verbose=false)[1] for f in files]
    mjd=[minimum(d.v2_mjd) for d in data]; ts=(mjd .- minimum(mjd)) .- PHASE0*PROT
    tess=tessellation_healpix(3;T=Float64); stars=create_star_multiepochs(tess,FIGURE,ts)
    setup_oi!(data,stars); npix=stars[1].npix
    seen=falses(npix); for s in stars; seen[s.index_quads_visible].=true; end
    nd=sum(d.nv2+d.nt3amp+d.nt3phi for d in data)
    pub=vec(readdlm(joinpath(ZET,"zetAnd2019A_map.txt")))
    flip=vec(Int.(readdlm(joinpath(ZET,"healpix_flip_3.txt"))))
    xpub=pub[flip]; xpub ./= mean(xpub)
    chi2n(x)= begin
        c=[ROTIR.chi2_breakdown(x,stars[i],data[i]) for i in eachindex(data)]
        (sum(ci.v2 for ci in c)+sum(ci.t3amp for ci in c)+sum(ci.t3phi for ci in c))/nd
    end
    rng=MersenneTwister(11); N=50
    acc=zeros(npix); rs=Float64[]; cs=Float64[]
    for k in 1:N
        x0 = 5000.0 .* (0.5 .+ rand(rng, npix))          # randomised start, positivity only
        x = image_reconstruct_oi(x0, data, stars; regularizers=Any[], maxiter=1000, verbose=false)
        x ./= mean(x); acc .+= x
        push!(rs, cor(x[seen], xpub[seen])); push!(cs, chi2n(x .* 5000.0))
    end
    m = acc ./ N
    @printf("individual members : chi2/n %.3f +- %.3f   r_pub %.3f +- %.3f   contrast %.2f-%.2f\n",
            mean(cs), std(cs), mean(rs), std(rs), minimum(m), maximum(m))
    @printf("ENSEMBLE MEAN of %d: chi2/n %.3f   r_pub %.4f   contrast %.3f-%.3f\n",
            N, chi2n(m .* 5000.0), cor(m[seen], xpub[seen]), minimum(m), maximum(m))
    @printf("published map      : contrast %.3f-%.3f\n", minimum(xpub), maximum(xpub))
    writedlm(joinpath(ZET,"zetAnd2019A_image_only_ensemble.txt"), m)
    println("wrote ensemble mean")
end
main()
