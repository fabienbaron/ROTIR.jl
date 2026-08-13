using LinearAlgebra, OptimPackNextGen
# The Fourier transform for polygonal surfaces
function poly_to_cvis(x, polyflux, polyft)
  flux = sum(polyflux.*x); # get the total flux
  cvis_model = polyft * x / flux;
end

"""
    poly_to_cvis(x, star; intensity_model=:linear, band=nothing)

Normalized complex visibilities of a single component.

`x` is the per-tessel map. With `intensity_model = :linear` (the historical default) it is
used directly as surface brightness — the Rayleigh–Jeans proxy, adequate in H and K.
With `:planck` it is treated as a *temperature* map and converted to band intensity by
[`intensity`](@ref); `band` is the wavelength in metres (see [`band_of`](@ref)). Use the
Planck form for R/V-band data, where the proxy misstates flux ratios by >10 %.
"""
function poly_to_cvis(x, star; intensity_model::Symbol = :linear, band = nothing)
   indx = star.index_quads_visible
   I = intensity_model === :linear ? x : intensity(x, intensity_model, band)
   xw = I[indx] .* star.vis_weights[indx] .* star.ldmap[indx]  # fold in limb darkening
   flux = dot(star.polyflux, xw); # get the total flux
   return star.polyft * xw / flux;
end

"""
    poly_to_flux(x, star; intensity_model=:linear, band=nothing)

Total (un-normalized) flux of a single component; see [`poly_to_cvis`](@ref) for the
intensity-model keywords.
"""
function poly_to_flux(x, star; intensity_model::Symbol = :linear, band = nothing)
   indx = star.index_quads_visible
   I = intensity_model === :linear ? x : intensity(x, intensity_model, band)
   xw = I[indx] .* star.vis_weights[indx] .* star.ldmap[indx]  # fold in limb darkening
   flux = dot(star.polyflux, xw); # get the total flux
   return flux;
end

# One epoch's flux/FT operators. Deliberately a separate function rather than the body of
# the threaded loop below: if `stars` is not concretely typed (e.g. a Vector{Any}), the
# loop body's locals become type-unstable, Julia boxes them into the closure, and the box
# is then SHARED by every thread — so one thread's `indx` can be read by another between
# the proj_west and proj_north lookups, producing a DimensionMismatch. A function call
# gives each invocation its own frame and cannot be boxed that way.
@views function _setup_oi_epoch!(star, dat, ::Type{T}; alt::Bool=false) where {T}
  indx = star.index_quads_visible
  pjx = star.proj_west[indx, :]
  pjy = star.proj_north[indx, :]
  star.polyflux = setup_polyflux_single(pjx, pjy)
  star.polyft = alt ? setup_polyft_single_alt(dat.uv, pjx, pjy; T=T) :
                      setup_polyft_single(dat.uv, pjx, pjy; T=T)
  return nothing
end

function setup_oi!(data, stars)
  nepochs = size(data,1);
  T = eltype(stars[1].vertices_xyz);
  if nepochs>1
    Threads.@threads for i=1:nepochs
      _setup_oi_epoch!(stars[i], data[i], T)
    end
  else # single epoch, thread over calculation
    _setup_oi_epoch!(stars[1], data[1], T; alt=true)
  end
end

# Same boxing hazard as _setup_oi_epoch! — keep the per-epoch work in its own function.
@views function _polygon_ft_epoch(star, dat, ::Type{T}) where {T}
  indx = star.index_quads_visible
  pjx = star.proj_west[indx, :]
  pjy = star.proj_north[indx, :]
  return setup_polyflux_single(pjx, pjy), setup_polyft_single(dat.uv, pjx, pjy; T=T)
end

function setup_polygon_ft(data, star_epoch_geom)
  nepochs = size(data,1);
  T = eltype(star_epoch_geom[1].vertices_xyz);
  polyflux = Array{Array{T,1}}(undef,nepochs);
  polyft = Array{Array{Complex{T},2}}(undef,nepochs);
  Threads.@threads for i=1:nepochs
    polyflux[i], polyft[i] = _polygon_ft_epoch(star_epoch_geom[i], data[i], T)
  end
  return polyflux, polyft
end

@views function setup_polyflux_single(pjx, pjy)
  # Polyflux is the projected area of each pixel (shoelace formula)
  T = eltype(pjx);
  polyflux = T(0.5)*(
    pjx[:,1].*pjy[:,2]
  - pjx[:,2].*pjy[:,1]
  + pjx[:,2].*pjy[:,3]
  - pjx[:,3].*pjy[:,2]
  + pjx[:,3].*pjy[:,4]
  - pjx[:,4].*pjy[:,3]
  + pjx[:,4].*pjy[:,1]
  - pjx[:,1].*pjy[:,4]);
  return polyflux;
end

function stcis(x1,x2,y1,y2,kx,ky)
  return sinc.(kx*(x2-x1) + ky*(y2-y1)).*cis.(-π*(kx*(x2+x1)+ ky*(y2+y1))).*(ky*(x2-x1)-kx*(y2-y1))
end

@views function setup_polyft_single_alt(uv, pjx, pjy; T = float(real(eltype(pjx))))
  kx =  uv[1,:] * T(-pi / (180*3600000))
  ky =  uv[2,:] * T( pi / (180*3600000))
  x = [Array(pjx[:,i])' for i in 1:4]
  y = [Array(pjy[:,i])' for i in 1:4]

  term1 = Threads.@spawn stcis(x[1], x[2], y[1], y[2], kx, ky)
  term2 = Threads.@spawn stcis(x[2], x[3], y[2], y[3], kx, ky)
  term3 = Threads.@spawn stcis(x[3], x[4], y[3], y[4], kx, ky)
  term4 = Threads.@spawn stcis(x[4], x[1], y[4], y[1], kx, ky)

  factor = -im * T(1/(2pi)) ./ (kx .* kx + ky .* ky)
  
  # Wait for the tasks to complete and get their results
  term1 = fetch(term1)
  term2 = fetch(term2)
  term3 = fetch(term3)
  term4 = fetch(term4)

  polyft = factor .* (term1 + term2 + term3 + term4)
  return polyft
end

@views function setup_polyft_single(uv, pjx, pjy; T = float(real(eltype(pjx))))
  kx =  uv[1,:] * T(-pi / (180*3600000));
  ky =  uv[2,:] * T( pi / (180*3600000));
  # note: check definition sinc(x) = sin(pi*x)/(pi*x)
  polyft = -im*T(1/(2pi))*(((sinc.( (kx*transpose(pjx[:,2]-pjx[:,1])) + (ky*transpose(pjy[:,2]-pjy[:,1])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,2]+pjx[:,1])) +  (ky*transpose(pjy[:,2]+pjy[:,1])) )).* ( (ky*transpose(pjx[:,2]-pjx[:,1]))  - (kx*transpose(pjy[:,2]-pjy[:,1])) ) + sinc.( (kx*transpose(pjx[:,3]-pjx[:,2])) + (ky*transpose(pjy[:,3]-pjy[:,2])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,3]+pjx[:,2])) +  (ky*transpose(pjy[:,3]+pjy[:,2])) )).* ( (ky*transpose(pjx[:,3]-pjx[:,2]))  - (kx*transpose(pjy[:,3]-pjy[:,2])) )+ sinc.( (kx*transpose(pjx[:,4]-pjx[:,3])) + (ky*transpose(pjy[:,4]-pjy[:,3])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,4]+pjx[:,3])) +  (ky*transpose(pjy[:,4]+pjy[:,3])) )).* ( (ky*transpose(pjx[:,4]-pjx[:,3]))  - (kx*transpose(pjy[:,4]-pjy[:,3])) ) + sinc.( (kx*transpose(pjx[:,1]-pjx[:,4])) + (ky*transpose(pjy[:,1]-pjy[:,4])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,1]+pjx[:,4])) +  (ky*transpose(pjy[:,1]+pjy[:,4])) )).* ( (ky*transpose(pjx[:,1]-pjx[:,4]))  - (kx*transpose(pjy[:,1]-pjy[:,4])) )))./((kx.*kx+ky.*ky)));
  return polyft;
end




@views function setup_polyft_single(kx, ky, pjx, pjy; T = float(real(eltype(pjx))))
  polyft = -im*T(1/(2pi))*(((sinc.( (kx*transpose(pjx[:,2]-pjx[:,1])) + (ky*transpose(pjy[:,2]-pjy[:,1])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,2]+pjx[:,1])) +  (ky*transpose(pjy[:,2]+pjy[:,1])) )).* ( (ky*transpose(pjx[:,2]-pjx[:,1]))  - (kx*transpose(pjy[:,2]-pjy[:,1])) ) + sinc.( (kx*transpose(pjx[:,3]-pjx[:,2])) + (ky*transpose(pjy[:,3]-pjy[:,2])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,3]+pjx[:,2])) +  (ky*transpose(pjy[:,3]+pjy[:,2])) )).* ( (ky*transpose(pjx[:,3]-pjx[:,2]))  - (kx*transpose(pjy[:,3]-pjy[:,2])) )+ sinc.( (kx*transpose(pjx[:,4]-pjx[:,3])) + (ky*transpose(pjy[:,4]-pjy[:,3])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,4]+pjx[:,3])) +  (ky*transpose(pjy[:,4]+pjy[:,3])) )).* ( (ky*transpose(pjx[:,4]-pjx[:,3]))  - (kx*transpose(pjy[:,4]-pjy[:,3])) ) + sinc.( (kx*transpose(pjx[:,1]-pjx[:,4])) + (ky*transpose(pjy[:,1]-pjy[:,4])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,1]+pjx[:,4])) +  (ky*transpose(pjy[:,1]+pjy[:,4])) )).* ( (ky*transpose(pjx[:,1]-pjx[:,4]))  - (kx*transpose(pjy[:,1]-pjy[:,4])) )))./((kx.*kx+ky.*ky)));
  return polyft;
end





function mod360(x)
  return mod.(mod.(x.+180,360).+360, 360) .- 180
end

function cvis_to_v2(cvis, indx)
  v2_model = abs2.(cvis[indx]);
end

# T follows the visibilities rather than defaulting to Float32: it sets the precision of the
# 180/π conversion, and Float32(180/π) = 57.2957802 is only 7 digits, which shifted χ²_t3phi
# by ~1e-8 relative against an otherwise-Float64 chain.
function cvis_to_t3(cvis, indx1, indx2, indx3; T = real(float(eltype(cvis))))
  t3 = cvis[indx1].*cvis[indx2].*cvis[indx3];
  t3amp = abs.(t3);
  t3phi = angle.(t3)*T(180/pi);
  return t3, t3amp, t3phi
end

# function cvis_to_t4(cvis, indx1, indx2, indx3, indx4)
#   t4 = cvis[indx1].*cvis[indx2]./(cvis[indx3]*conj(cvis[indx4]));
#   t4amp = abs.(t4);
#   t4phi = angle.(t4)*180.0/pi;
#   return t4, t4amp, t4phi
# end

"""
    cvis_to_obs(cvis, data) -> (v2, t3amp, t3phi)

Convert complex visibilities to interferometric observables (V2, T3amp, T3phi).
Common step for both single-star and binary forward models.
"""
function cvis_to_obs(cvis, data)
  v2_model = cvis_to_v2(cvis, data.indx_v2);
  _, t3amp_model, t3phi_model = cvis_to_t3(cvis, data.indx_t3_1, data.indx_t3_2, data.indx_t3_3);
  return v2_model, t3amp_model, t3phi_model
end

"V², T3amp, T3phi on; visamp, visphi, flux, differential phase off."
const OI_DEFAULT_WEIGHTS = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]

"""
    cvis_chi2(cvis, data, weights = OI_DEFAULT_WEIGHTS) -> chi2

χ² of complex visibilities against an `OIdata` block, **with a hand-written adjoint**.

Thin wrapper over OITOOLS' `cvis_to_chi2_fg`, which returns the value *and* the complex
adjoint source `g_cvis` from a single pass over the data. The `rrule` below simply hands
that back, so reverse-mode AD through this call costs one hand-written sweep rather than a
tape of `abs2`/`angle`/`mod360` broadcasts.

Use this instead of `cvis_to_obs` + three `sum(abs2, …)` whenever the χ² is going to be
differentiated. It is the analytic-component counterpart of
[`interferometric_chi2`](@ref), whose rrule is fused to the polygon-FT image path and so
cannot be reused here.

# Adjoint convention

OITOOLS documents `g_params = real(Jᵀ · g_cvis)` with `J = ∂V/∂p` — a **transpose, not an
adjoint**, i.e. no conjugation. ChainRules/Zygote instead want the cotangent
`∂L/∂Re(cvis) + i ∂L/∂Im(cvis)`, which is the *complex conjugate* of that. Hence the
`conj(g)` in the pullback.

Getting this backwards is not a uniformly wrong answer, which is what makes it dangerous:
it leaves flux-like parameters nearly right (they move `V` along the real direction) while
making phase-like parameters — separation, position angle, T0 — wrong by ~100%. Verified
against finite differences to 1e-11 in both directions; `test_reflection.jl` pins it.
"""
# Pass the visibilities through at whatever complex float type they already are — do NOT
# force ComplexF64. OITOOLS accepts `AbstractVector{<:Complex}` and (in recent versions)
# derives its working eltype from the inputs, so widening here would throw away the point
# of carrying Float32 through the forward model.
@inline _ascomplexfloat(v::AbstractVector{<:Complex{<:AbstractFloat}}) = v
@inline _ascomplexfloat(v::AbstractVector) = complex.(float.(v))

cvis_chi2(cvis, data) = cvis_chi2(cvis, data, OI_DEFAULT_WEIGHTS)
# Value-only path uses `cvis_to_chi2_f`, NOT `first(cvis_to_chi2_fg(...))`: the latter also
# builds and fills the adjoint source, which cost ~40 % on the primal for nothing whenever
# no gradient was wanted. The rrule below still uses `_fg`, which returns both in one sweep.
cvis_chi2(cvis, data, weights) =
    cvis_to_chi2_f(_ascomplexfloat(cvis), data; weights = weights)

function ChainRulesCore.rrule(::typeof(cvis_chi2), cvis, data, weights)
    chi2, g = cvis_to_chi2_fg(_ascomplexfloat(cvis), data; weights = weights)
    function cvis_chi2_pullback(c̄raw)
        return (NoTangent(), unthunk(c̄raw) .* conj(g), NoTangent(), NoTangent())
    end
    return chi2, cvis_chi2_pullback
end

"""
    observables(x, star, data; intensity_model=:linear, band=nothing) -> (v2, t3amp, t3phi)

Model observables for a single component. See [`poly_to_cvis`](@ref) for the
intensity-model keywords.
"""
function observables(x, star, data; intensity_model::Symbol = :linear, band = nothing)
  cvis_model = poly_to_cvis(x, star; intensity_model=intensity_model, band=band);
  return cvis_to_obs(cvis_model, data)
end

function chi2s(x, star, data; verbose::Bool = true,
               intensity_model::Symbol = :linear, band = nothing)
  v2_model, t3amp_model, t3phi_model = observables(x, star, data;
                                                   intensity_model=intensity_model, band=band);
  chi2_v2 = sum(abs2, (v2_model - data.v2)./data.v2_err);
  chi2_t3amp = sum(abs2, (t3amp_model - data.t3amp)./data.t3amp_err);
  chi2_t3phi = sum(abs2, mod360(t3phi_model - data.t3phi)./data.t3phi_err);
  if verbose == true
    printstyled(@sprintf("V2: %.4f ", chi2_v2/data.nv2), color=:red)
    printstyled(@sprintf("T3A: %.4f ", chi2_t3amp/data.nt3amp), color=:blue)
    printstyled(@sprintf("T3P: %.4f\n", chi2_t3phi/data.nt3phi), color=:green)
  end
  return chi2_v2, chi2_t3amp, chi2_t3phi
end

function chi2s2(x, star, data; verbose::Bool = true)
  v2_model, t3amp_model, t3phi_model = observables(x, star, data);
  chi2_v2 = n2((v2_model - data.v2)./data.v2_err);
  chi2_t3amp = n2((t3amp_model - data.t3amp)./data.t3amp_err);
  chi2_t3phi = n2( mod360(t3phi_model - data.t3phi)./data.t3phi_err);
  if verbose == true
    printstyled(@sprintf("V2: %.4f ", chi2_v2/data.nv2), color=:red)
    printstyled(@sprintf("T3A: %.4f ", chi2_t3amp/data.nt3amp), color=:blue)
    printstyled(@sprintf("T3P: %.4f\n", chi2_t3phi/data.nt3phi), color=:green)
  end
  return chi2_v2, chi2_t3amp, chi2_t3phi
end

function spheroid_chi2_f(x, star, data; verbose::Bool = false,
                         intensity_model::Symbol = :linear, band = nothing)
  v2_model, t3amp_model, t3phi_model = observables(x, star, data;
                                                   intensity_model=intensity_model, band=band);
  chi2_v2 = sum(abs2, (v2_model - data.v2)./data.v2_err);
  chi2_t3amp = sum(abs2, (t3amp_model - data.t3amp)./data.t3amp_err);
  chi2_t3phi = sum(abs2, mod360(t3phi_model - data.t3phi)./data.t3phi_err);
  if verbose == true
    printstyled(@sprintf("V2: %.4f ", chi2_v2/data.nv2), color=:red)
    printstyled(@sprintf("T3A: %.4f ", chi2_t3amp/data.nt3amp), color=:blue)
    printstyled(@sprintf("T3P: %.4f\n", chi2_t3phi/data.nt3phi), color=:green)
  end
  return chi2_v2 + chi2_t3amp + chi2_t3phi
end

@views function spheroid_chi2_fg(x, g, star, data; verbose::Bool = true)
  npix = star.npix;
  T = eltype(x);
  indx = star.index_quads_visible
  w = star.vis_weights[indx] .* star.ldmap[indx]  # soft visibility × limb darkening
  xw = x[indx] .* w           # weighted pixel values
  cvis_model = poly_to_cvis(x, star);
  v2_model = cvis_to_v2(cvis_model, data.indx_v2);
  t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
  chi2_v2 = sum(abs2, (v2_model - data.v2)./data.v2_err);
  chi2_t3amp = sum(abs2, (t3amp_model - data.t3amp)./data.t3amp_err);
  chi2_t3phi = sum(abs2, mod360(t3phi_model - data.t3phi)./data.t3phi_err);
  # Gradient w.r.t. weighted pixels (xw), computed via polyft^T @ adjoint_signal
  g_v2 = real(transpose(star.polyft[data.indx_v2,:])*(4*((v2_model-data.v2)./data.v2_err.^2).*conj(cvis_model[data.indx_v2])));
  g_t3amp = real(transpose(star.polyft[data.indx_t3_1,:]) *(2*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*conj(cvis_model[data.indx_t3_1])./abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_3]) ))+real(transpose(star.polyft[data.indx_t3_2,:])*(2*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*conj(cvis_model[data.indx_t3_2])./abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_3]) ))+real(transpose(star.polyft[data.indx_t3_3,:])*(2*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*conj(cvis_model[data.indx_t3_3])./abs.(cvis_model[data.indx_t3_3]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]) ));
  g_t3phi = T(360/pi)*imag(transpose(star.polyft[data.indx_t3_1,:])*(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*cvis_model[data.indx_t3_2].*cvis_model[data.indx_t3_3].*conj(t3_model))+transpose(star.polyft[data.indx_t3_2,:])*(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_3].*conj(t3_model))+transpose(star.polyft[data.indx_t3_3,:])*(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_2].*conj(t3_model)));
  gsum = g_v2 + g_t3amp + g_t3phi;
  flux = poly_to_flux(x, star);
  # Gradient w.r.t. x: chain rule through soft visibility weights
  # gsum is ∂χ²/∂(xw), so ∂χ²/∂x = w .* ∂χ²/∂(xw) (after flux normalization)
  g_normalized = (gsum .- dot(xw, gsum) * star.polyflux / flux) / flux;
  g[indx] = w .* g_normalized;
  if verbose == true
    printstyled(@sprintf("V2: %.4f ", chi2_v2/data.nv2), color=:red)
    printstyled(@sprintf("T3A: %.4f ", chi2_t3amp/data.nt3amp), color=:blue)
    printstyled(@sprintf("T3P: %.4f ", chi2_t3phi/data.nt3phi), color=:green)
    printstyled(@sprintf("Flux: %.4f\n", flux), color=:normal)
  end
  return chi2_v2 + chi2_t3amp + chi2_t3phi
end

function spheroid_chi2_allepochs_f(x, stars, data; epochs_weights=[], verbose=false)
nepochs = length(data)
chi2_t = zeros(eltype(x), nepochs);
Threads.@threads for i=1:nepochs # weighted sum -- should probably do the computation in parallel
  chi2_t[i] = spheroid_chi2_f(x, stars[i], data[i], verbose=verbose);
end
f = sum(chi2_t)
if epochs_weights!=[]
  f = f.*epochs_weights
end
if verbose == true 
  printstyled(@sprintf("All epochs, χ²r: %.4f\n\n", f), color=:white);
end
return f;
end

function spheroid_crit_allepochs_fg(x, g, stars, data; regularizers=[], epochs_weights=[],
                                    verbose=false, T=eltype(x))
  # (the keyword used to default to Float32 and then be overwritten by `T = eltype(x)` on the
  #  next line, so anything a caller passed was silently discarded)
  #g[:] .= T(0);
  nepochs = length(data)
  chi2_t = zeros(T, nepochs);
  #npix = stars[1].npix
  singleepoch_g = [zeros(T, length(x)) for i=1:nepochs];
  Threads.@threads for i=1:nepochs # weighted sum -- should probably do the computation in parallel
    chi2_t[i] = spheroid_chi2_fg(x, singleepoch_g[i], stars[i], data[i], verbose=verbose);
  end
  f = sum(chi2_t)
  if epochs_weights!=[]
  #  f = f.*epochs_weights
    @warn "Epoch weights not implemented"
   end 
  g[:] .= sum(singleepoch_g)
  if verbose == true
    printstyled(@sprintf("Total χ²r: %.4f\n", f), color=:white);
  end
  # Map regularization
  if regularizers!=[]
    reg_g = zeros(T, length(x));
    f += spheroid_regularization(x, reg_g, regularizers=regularizers, verbose = verbose);
    g[:] += reg_g
  end
  return f;
end

function parametric_temperature_map(parameters, star; secondary=false) # dispatches parametric
  if star.surface_type == 3
    return temperature_map_vonZeipel_roche_single(parameters, star, star.t, secondary=secondary);
  elseif star.surface_type == 2
    return temperature_map_vonZeipel_rapid_rotator(parameters,star);
  elseif star.surface_type == 1
    return temperature_map_vonZeipel_ellipsoid(parameters,star);
  elseif star.surface_type == 0 # sphere
    T = eltype(star.vertices_spherical)
    return T(parameters.tpole)*ones(T, star.npix);
  else
    println("Unimplemented parametric von Zeipel function")
  end
  return
end

function spheroid_parametric_f(parameters, tessels, data, tepochs) 
  stars = create_star_multiepochs(tessels, parameters, tepochs);
  setup_oi!(data, stars) 
  x = parametric_temperature_map(parameters, stars[1]);
  return spheroid_chi2_allepochs_f(x, stars, data)
end

function proj_positivity(ztilde)
z = copy(ztilde)
z[ztilde.>0]=0
return z
end


# function spheroid_total_variation_fg(x, tv_g, tvinfo; ϵ = 1e-13, T=Float32, verbose = true)
#   # Add total variation regularization
#   #tvinfo: neighbors,south_neighbors,west_neighbors,south_neighbors_reverse,west_neighbors_reverse
#   ϵ = T(ϵ)
#   npix = length(x)
#   xs = x[tvinfo[2]];
#   xw = x[tvinfo[3]];
#   tv_f = sum(sqrt.( (x-xs).^2 + (x-xw).^2) .+ ϵ )
#   tv_g[1:npix] = (2*x-xs-xw)./(sqrt.( (x-xs).^2 + (x-xw).^2) .+ ϵ)
#   for j=1:length(x)
#     k = tvinfo[4][j]
#     l = tvinfo[5][j]
#     if length(k)>0
#       kk=k[1]
#       tv_g[j] += (x[j]-x[kk])/sqrt((x[kk]-x[j])^2+(x[kk]-x[tvinfo[3][kk]])^2 +ϵ)
#     end
#     if length(l)>0
#       ll=l[1]
#       tv_g[j] += (x[j]-x[ll])/sqrt((x[ll]-x[tvinfo[2][ll]])^2+(x[ll]-x[j])^2 +ϵ)
#     end
#    end
# if verbose == true
#       println("TV: ", tv_f);
# end
#   return tv_f
# end


function spheroid_total_variation2_fg(x, tv_g, tvinfo; verbose = true)
  npix = length(x)
  tv_f = norm(tvinfo[6]*x)^2
  tv_g[:] = 2*(tvinfo[7]*x)
  if verbose == true
      printstyled(@sprintf("TV2: %.4f\n", tv_f), color=:yellow);
  end
  return tv_f
end

function spheroid_total_variation_fg(x, tv_g, tvinfo; verbose = true)
  npix = length(x)
  tv_f = norm(tvinfo[6]*x)
  if tv_f>0 
  tv_g[:] = (tvinfo[7]*x)/tv_f
  else
    tv_g[:] .= 0
  end
  if verbose == true
      println("TV: ", tv_f);
  end
  return tv_f
end

function spheroid_mean_fg(x, g; verbose = true)
f = sum(abs.(x.-sum(x)/length(x)))
g[:] = sign.(x.-sum(x)/length(x))
if verbose == true
println(" MeanReg: ", f);
end
return f;
end

function spheroid_harmon_bias_fg(x, g, B::Float64; verbose = true)
n = length(x);
avg_x = mean(x);
bcorr = (B-1.0)*Int.((x.-avg_x).>0).+1.0
#reg_f = sum(bcorr.*(x.-avg_x).^2);
reg_f = sum(bcorr.*(x.-avg_x).^2)/n;
#reg_g = 2*(n-1)/n*bcorr.*(x.-avg_x)
reg_g = 2/n*bcorr.*(x.-avg_x)
g[:] = reg_g .- mean(reg_g); # necessary ?

if verbose == true
println(" Bias: ", reg_f);
end
return reg_f;
end

function max_entropy_fg(x, g; verbose=false, ϵ=1e-9)
  # mmap = sum(x) / length(x)
  # nmap = x ./ mmap
  # reg_f = sum(nmap .* log.(nmap))
  # g[:] = sum(-nmap ./ sum(x) .* (log.(nmap) .+ 1.0)) .+ (log.(nmap) .+ 1)./mmap
  xm = x ./ (mean(x) + ϵ)
  reg_f = sum(xm .* log.(xm .+ ϵ))
  g[:] = ((mean(x) .- (x ./ length(x))) ./ (mean(x)^2 + ϵ)) .* (log.(xm .+ ϵ) .+ 1)
  if verbose == true
    println(" MEM: ", reg_f)
  end
  return reg_f
end

"""
    radflat_bins(star; nbins=6, subset=nothing) -> NamedTuple

Precompute the projected-radial binning that [`spheroid_radflat_fg`](@ref) needs. Pass the
result as the third element of a `"radflat"` entry in `regularizers`.

Patches are binned by the projected radius `ρ` of their centre (mean of the four projected
corners), normalised so the outermost visible patch sits at `ρ = 1`. Returns the per-patch
bin index, the per-bin weight `w_b = ρ_b² · Σ polyflux`, and the patch counts.

`ρ_b²` weights the limb, which is the whole point: that is where a spot is degenerate with
the limb-darkening coefficient. `flux_weight_b` (the projected area) stops a bin containing
a handful of slivers from carrying the same force as one containing half the disk.
"""
function radflat_bins(star; nbins::Int = 6, subset = nothing)
    idx = subset === nothing ? star.index_quads_visible : subset
    pw  = @view star.proj_west[idx, :]
    pn  = @view star.proj_north[idx, :]
    n   = length(idx)
    ρ   = [hypot(sum(@view pw[i, :])/4, sum(@view pn[i, :])/4) for i in 1:n]
    ρmax = maximum(ρ)
    ρmax > 0 || error("radflat_bins: all projected radii are zero")
    ρn = ρ ./ ρmax
    # bin edges equispaced in ρ; the outermost patch lands in the last bin, not one past it
    b  = [min(nbins, floor(Int, r*nbins) + 1) for r in ρn]
    pf = length(star.polyflux) == n ? abs.(star.polyflux) : abs.(star.polyflux[idx])
    ρc = [(k - 0.5)/nbins for k in 1:nbins]                 # bin centre radius
    w  = zeros(Float64, nbins)
    cnt = zeros(Int, nbins)
    @inbounds for i in 1:n
        w[b[i]]  += pf[i]
        cnt[b[i]] += 1
    end
    w .= (ρc .^ 2) .* w                                     # ρ_b² · flux_weight_b
    # Under-populated bins are the failure mode at coarse tessellations. Equal-ρ bins on a
    # projected disk hold ∝ ρ dρ patches, so the INNERMOST bin starves first: at HEALPix
    # nside=1 (24 visible patches) one bin comes out empty and the penalty is ~30 % wrong;
    # at nside=2 the innermost bin holds a single patch, so its "mean brightness" is one
    # number. An empty bin contributes zero rather than erroring — quietly, which is exactly
    # the kind of silent wrong answer worth a warning. nside ≥ 3 is stable to ~1 %.
    if any(iszero, cnt)
        @warn "radflat_bins: $(count(iszero, cnt)) of $nbins radial bins are EMPTY " *
              "($(cnt)); RADFLAT will be biased. Use a finer tessellation or fewer bins."
    elseif minimum(cnt) < 3
        @warn "radflat_bins: thinnest radial bin holds only $(minimum(cnt)) patch(es) " *
              "($(cnt)); RADFLAT will be noisy. Consider a finer tessellation or fewer bins."
    end
    return (idx = idx, bin = b, nbins = nbins, w = w, count = cnt, rho = ρn)
end

"""
    spheroid_radflat_fg(x, g, bins; scale=1e5, verbose=false) -> f

RADFLAT (J. Monnier, priv. comm. 2026-07-01): force the *azimuthally averaged* radial
brightness profile of the visible disk to be flat.

```
f = scale · Σ_b w_b · (I_b/I_mean − 1)²
```

with `I_b` the mean patch brightness in projected radial bin `b` and `I_mean` the mean over
the whole visible disk. `x` is the patch brightness **without** limb darkening — that is
deliberate, since the point is to stop the reconstruction from smuggling limb structure into
the surface map where it trades off against the LD coefficient.

# Why it exists

On a **non-rotating** star reconstructed from few epochs, the surface regularizer is weak
enough that dark or bright spots drift to the limb and become degenerate with the
limb-darkening parameter: χ² barely moves while the LD coefficient wanders. Monnier's RW Cep
test gave `α = 1.07 ± 0.19` uncontrolled against `α = 0.49 ± 0.02` at strong RADFLAT, with
χ² only slightly worse and the image essentially unchanged.

!!! note "Not for rotators"
    A rotating star observed over several epochs already breaks that degeneracy — spots move
    and the limb is resampled. RADFLAT would then be fighting real structure. Use it for the
    single-epoch / non-rotating case it was designed for.

# Gradient

With `n_b` patches in bin `b` and `N` visible patches, `∂I_b/∂x_j = [j∈b]/n_b` and
`∂I_mean/∂x_j = 1/N`, so

```
∂f/∂x_j = 2·scale·[ w_{b(j)}(r_{b(j)} − 1)/(n_{b(j)}·I_mean)
                    − (1/(N·I_mean))·Σ_b w_b (r_b − 1) r_b ],    r_b ≡ I_b/I_mean
```

the second term being a global coupling shared by every visible patch — it is what stops the
regularizer from simply rescaling the whole map.
"""
function spheroid_radflat_fg(x, g, bins; scale = 1e5, verbose = false, ϵ = 1e-12)
    nb = bins.nbins; b = bins.bin; w = bins.w; cnt = bins.count
    N  = length(x)
    N == length(b) ||
        error("spheroid_radflat_fg: x has $N entries but the binning was built for $(length(b)). " *
              "The regularizer's pixel subset must match the one passed to radflat_bins.")
    Imean = sum(x)/N
    abs(Imean) > ϵ || (g .= 0; return 0.0)
    S = zeros(Float64, nb)
    @inbounds for i in 1:N; S[b[i]] += x[i]; end
    r = [cnt[k] > 0 ? (S[k]/cnt[k])/Imean : 1.0 for k in 1:nb]   # empty bin contributes nothing
    f = scale * sum(w[k]*(r[k]-1)^2 for k in 1:nb)
    glob = sum(w[k]*(r[k]-1)*r[k] for k in 1:nb) / (N*Imean)
    @inbounds for i in 1:N
        k = b[i]
        loc = cnt[k] > 0 ? w[k]*(r[k]-1)/(cnt[k]*Imean) : 0.0
        g[i] = 2*scale*(loc - glob)
    end
    verbose && println(" RADFLAT: ", f)
    return f
end

function spheroid_regularization(x,g; printcolor = :black, regularizers=[], verbose=false)
  reg_f = 0.0;
  for ireg =1:length(regularizers)
      x_sub = x[regularizers[ireg][4]] # take the pixel subset if needed (example: only regularize visible pixels)
      temp_g = similar(x_sub);
      if regularizers[ireg][1] == "mem"
          reg_f += regularizers[ireg][2]*max_entropy_fg(x_sub, temp_g, verbose = verbose);
      elseif regularizers[ireg][1] == "tv2"
          reg_f += regularizers[ireg][2]*spheroid_total_variation2_fg(x_sub, temp_g, regularizers[ireg][3], verbose = verbose);
      elseif regularizers[ireg][1] == "tv"
          reg_f += regularizers[ireg][2]*spheroid_total_variation_fg(x_sub, temp_g, regularizers[ireg][3], verbose = verbose);
      elseif regularizers[ireg][1] == "mean"
          reg_f += regularizers[ireg][2]*spheroid_mean_fg(x_sub, temp_g, verbose = verbose);
      elseif regularizers[ireg][1] == "bias"
          reg_f += regularizers[ireg][2]*spheroid_harmon_bias_fg(x_sub, temp_g, regularizers[ireg][3], verbose = verbose );
      elseif regularizers[ireg][1] == "radflat"
          reg_f += regularizers[ireg][2]*spheroid_radflat_fg(x_sub, temp_g, regularizers[ireg][3], verbose = verbose);
      else
          # Previously this chain had no `else`: an unrecognised name contributed exactly
          # zero, silently, so a typo looked like a regularizer that simply did nothing.
          error("spheroid_regularization: unknown regularizer \"$(regularizers[ireg][1])\"; " *
                "known: mem, tv, tv2, mean, bias, radflat")
      end
      g[regularizers[ireg][4]] += regularizers[ireg][2]*temp_g
  end
  if verbose == true
      println("\n");
  end
return reg_f
end


using OptimPackNextGen

function image_reconstruct_oi(x_start, data, stars; epochs_weights =[], printcolor= [], verbose = true, lower=0, upper=Inf, maxiter = 100, regularizers =[])
  x_sol = [];
  crit_imaging = (x,g)->spheroid_crit_allepochs_fg(x, g, stars, data, regularizers=regularizers, epochs_weights=  epochs_weights, verbose = verbose);
  x_sol = OptimPackNextGen.vmlmb(crit_imaging, x_start, verb=verbose, lower=lower, upper=upper, maxiter=maxiter, blmvm=false, gtol=(0,1e-8));
  dummy = similar(x_sol);
  crit_opt = crit_imaging(x_sol,dummy);
  return x_sol
end

function image_reconstruct_oi_crit(x, data, stars; regularizers =[],  verbose = verbose)
  g = similar(x);
  crit = spheroid_crit_allepochs_fg(x, g, stars, data, regularizers=regularizers,epochs_weights=[],  verbose = verbose);
  return crit
end

function image_reconstruct_oi_chi2(x, data, stars;  verbose = verbose)
  g = similar(x);
  crit = spheroid_crit_allepochs_fg(x, g, stars, data, regularizers=[], epochs_weights= [], verbose = verbose);
  return crit
end

function image_reconstruct_oi_chi2_fg(x, data, stars;  verbose = verbose)
  g = similar(x);
  crit = spheroid_crit_allepochs_fg(x, g, stars, data, regularizers=[], epochs_weights= [], verbose = verbose);
  return crit,g
end

function multires_reconstruct_oi(data, star_params, tepochs; n_start=2, n_end=4, maxiter=500, reg_weight=1e-5, reg_type="tv2", verbose=true, kwargs...)
  tmap = nothing
  stars = nothing
  for n in n_start:n_end
    tessels = tessellation_healpix(n)
    stars = create_star_multiepochs(tessels, star_params, tepochs)
    if tmap === nothing
      tmap = parametric_temperature_map(star_params, stars[1])
    else
      tmap = vec(repeat(tmap, 1, 4)')  # upsample from previous level
    end
    setup_oi!(data, stars)
    regularizers = [[reg_type, reg_weight, tv_neighbors_healpix(n), 1:length(tmap)]]
    if verbose
      println("Multi-resolution: HEALPix level n=$n, npix=$(nside2npix(2^n))")
    end
    tmap = image_reconstruct_oi(tmap, data, stars; maxiter=maxiter, regularizers=regularizers, verbose=verbose, kwargs...)
  end
  return tmap, stars
end

"""
    rescale_temperature_tpole(tmap, stars, star_params)

Rescale a temperature map so that the temperature at the pole equals `star_params.tpole`.
`stars` can be a single `stellar_geometry` or a vector of them (multi-epoch).
Warns if the pole pixel is never visible across all epochs.
"""
function rescale_temperature_tpole(tmap, stars, star_params)
    star = stars isa AbstractVector ? stars[1] : stars
    colat = @view star.vertices_spherical[:, 5, 2]
    ipole = argmin(colat)
    # Check visibility
    if stars isa AbstractVector
        visible = any(ipole in s.index_quads_visible for s in stars)
    else
        visible = ipole in star.index_quads_visible
    end
    if !visible
        @warn "Pole pixel (index $ipole) is never visible — rescaling is based on an unconstrained pixel"
    end
    scale = star_params.tpole / tmap[ipole]
    return tmap .* scale
end

"""
    rescale_temperature_teff(tmap, stars, teff)

Rescale a temperature map so that the mean temperature of the visible pixels equals `teff`.
`stars` can be a single `stellar_geometry` or a vector of them (multi-epoch).
Visibility is the union across all epochs.
"""
function rescale_temperature_teff(tmap, stars, teff)
    vis = sometimes_visible(stars)
    if isempty(vis)
        @warn "No visible pixels — returning unmodified map"
        return copy(tmap)
    end
    mean_T = mean(tmap[vis])
    scale = teff / mean_T
    return tmap .* scale
end


# "Multitemporal" = dynamical reconstruction

# function oi_reconstruct_mutitemporal(x_start::Array{Float64,1}, data::Array{OIdata,1}, polyflux::Array{Array{Float64,1},1}, polyft::Array{Array{Complex{Float64},2},1}; epochs_weights =[], printcolor= [], verbose = true, maxiter = 100, regularizers =[])
#   crit_imaging = (x,g)->oi_multitemporal_fg(x, g, polyflux, polyft, data, regularizers=regularizers, epochs_weights=  epochs_weights);
#   x_sol = OptimPackNextGen.vmlmb(crit_imaging, x_start, verbose=verbose, lower=0, maxiter=maxiter, blmvm=false, gtol=(0,1e-8));
#   return reshape(x_sol, size(polyflux[1],1), length(data))
# end

# function oi_multitemporal_fg(x, g, polyflux, polyft, data; regularizers=[], epochs_weights = [], verbose=true)
#   # Explanation of the following: optimpack optimizes a vector, but we want an array of images
#   npix = size(polyflux[1],1);
#   chi2_f = 0.0;
#   g[:] .= 0.0;
#   #temp_g = deepcopy(g[:]);
#   #npix = size(x);
#   nepochs = length(data);
#   if epochs_weights == []
#     epochs_weights = ones(Float64, nepochs)/nepochs;
#   end

#   #singleepoch_g = zeros(Float64, npix);
#   #for i=1:nepochs # weighted sum -- should probably do the computation in parallel
#   #  chi2_f += epochs_weights[i]*spheroid_chi2_fg_alt(x[1:npix], singleepoch_g, polyflux[i], polyft[i], data[i], verbose=verbose);
#   #  g[1:npix] += epochs_weights[i]*singleepoch_g;
#   #end


#   for i=1:nepochs # weighted sum -- in the future, do the computation in parallel
#     tslice = 1+(i-1)*npix:i*npix; # temporal slice
#     subg = zeros(Float64, npix);
#     chi2_f += epochs_weights[i]*spheroid_chi2_fg_alt(x[tslice], subg, polyflux[i], polyft[i], data[i], verbose=verbose);
#     g[tslice] = epochs_weights[i]*subg;
#     #x[tslice] = x[1:npix];
#     #g[tslice] = g[1:npix];
#   end

#   # cross temporal regularization -- weight needs to be defined in the "regularizers" variable
#   if length(regularizers)>nepochs
#     if (regularizers[nepochs+1][1][1] == "temporal_tvsq")  & (regularizers[nepochs+1][1][2] > 0.0) & (nepochs>1)
#       y = reshape(x,(npix,nepochs))
#       temporalf = sum( (y[:,2:end]-y[:,1:end-1]).^2 )
#       tv_g = Array{Float64}(undef, npix,nepochs)
#       if nepochs>2
#          tv_g[:,1] = 2*(y[:,1] - y[:,2])
#          tv_g[:,2:end-1] = 4*y[:,2:end-1]-2*(y[:,1:end-2]+y[:,3:end])
#          tv_g[:,end] = 2*(y[:,end] - y[:,end-1])
#       else
#          tv_g[:,1] = 2*(y[:,1]-y[:,2]);
#          tv_g[:,2] = 2*(y[:,2]-y[:,1]);
#       end
#       chi2_f += regularizers[nepochs+1][1][2]*temporalf
#       g[:] += regularizers[nepochs+1][1][2]*vec(tv_g);
#       if verbose == true
#            printstyled("Temporal regularization: $temporalf\n", color=:yellow)
#       end
#     end
#   end
#   reg_f = spheroid_regularization(x, g, regularizers=regularizers[1:nepochs], verbose = verbose); # Note: adds to g
#   return chi2_f + reg_f;
# end
