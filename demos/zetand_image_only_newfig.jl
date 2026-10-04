# Image-only reconstruction on the figure that her map prefers (chi2/n = 2.064), rather than on
# the figure a UNIFORM surface prefers. Same data, same prior, better geometry.
using ROTIR, Printf, DelimitedFiles, Statistics
const PROT=17.769426; const ZET=joinpath(pkgdir(ROTIR),"demos","zetAnd")
# re-optimised against her fixed map, with T3phi in
const FIG=(surface_type=1, radius_x=1.3332, radius_y=1.2545, radius_z=1.2269, tpole=5000.0,
           ldtype=3, ld1=0.145, ld2=0.0, beta=0.25, inclination=71.3, position_angle=124.3,
           rotation_period=PROT)
const PHASE0=0.5061
function main()
    files=sort(filter(f->endswith(f,".oifits"),readdir(ZET)))
    data=[readoifits(joinpath(ZET,f);verbose=false)[1] for f in files]
    mjd=[minimum(d.v2_mjd) for d in data]; ts=(mjd .- minimum(mjd)) .- PHASE0*PROT
    tess=tessellation_healpix(3;T=Float64); stars=create_star_multiepochs(tess,FIG,ts)
    setup_oi!(data,stars); npix=stars[1].npix
    seen=falses(npix); for s in stars; seen[s.index_quads_visible].=true; end
    nd=sum(d.nv2+d.nt3amp+d.nt3phi for d in data)
    pub=vec(readdlm(joinpath(ZET,"zetAnd2019A_map.txt")))
    flip=vec(Int.(readdlm(joinpath(ZET,"healpix_flip_3.txt"))))
    xpub=pub[flip]; xpub ./= mean(xpub)
    chi2n(x)=begin
        c=[ROTIR.chi2_breakdown(x,stars[i],data[i]) for i in eachindex(data)]
        (sum(ci.v2 for ci in c)+sum(ci.t3amp for ci in c)+sum(ci.t3phi for ci in c))/nd
    end
    @printf("her map on this figure: chi2/n %.3f\n", chi2n(xpub .* 5000.0))
    sob=sobel_gradient_healpix(3;T=Float64)
    for w in (1e4,)
        x=image_reconstruct_oi(fill(5000.0,npix), data, stars;
                               regularizers=Any[["sobel",w,sob,1:npix]], maxiter=1000, verbose=false)
        x ./= mean(x)
        @printf("sobel %6.0f : chi2/n %.3f  r_pub %.4f  contrast %.3f-%.3f\n",
                w, chi2n(x .* 5000.0), cor(x[seen],xpub[seen]), minimum(x), maximum(x))
        w == 1e4 && writedlm(joinpath(ZET,"zetAnd2019A_image_only_newfig.txt"), x)
    end
    println("wrote zetAnd2019A_image_only_newfig.txt (sobel 3000)")
end
main()
