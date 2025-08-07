
function setup_fullft!(uv, stars)
    nepochs = size(stars,1);
    T = eltype(stars[1].vertices_xyz);
    for i=1:nepochs
         stars[i].polyflux = setup_polyflux_single(stars[i])
         stars[i].polyft = setup_polyft_single(uv, stars[i]);
         indx = findall( (uv[2,:].==0) .&& (uv[1,:].==0))
         stars[i].polyft[indx,:]=ComplexF32.(stars[i].polyflux)
    end
end

function rfft_size(L::Int) # Used to find N for RFFT size = (N ÷ 2 + 1,  N)
    # Try the even case: N = 2m
    d1 = 1 + 2L
    r1 = sqrt(d1)
    if isinteger(r1)
        m = (-1 + r1) ÷ 2
        N = 2 * m
        return Int(N)
    end
    # Try the odd case: N = 2m + 1
    d2 = 8L + 1
    r2 = sqrt(d2)
    m = (-3 + r2) ÷ 4
    N = 2 * m + 1
    return Int(N)
end


function poly_to_img(x, stars)
    V = stars[1].polyft*x
    nx_samples = rfft_size(length(V))
    VV = reshape(V, (nx_samples÷2+1, nx_samples))
    model_img = FFTW.ifftshift(FFTW.irfft(VV,nx_samples))    
    model_img .*= model_img.>0
    return model_img
end

# function crit_fg(x,g,ft,ftT,data)
#     res = ft*x - data
#     f   = sum(abs2,res)
#     g[:] .= 2*real(ftT*res)
#     return f
# end

function crit_single_obj_fg(x,g,ft,ftT,psf,data)
    nx = size(psf,1)
    V = ft*x
    res = convolve_fo(psf,V) - data
    f   = sum(abs2,res)
    g[:] .= 2*real(ftT * vec(conj.(rfft(psf)) .* rfft(res)))
    return f
end

function crit_single_psf_fg(psf,g,V,data)
    res = convolve_fo(psf,V) .- data
    f   = sum(abs2,res)
    g[:] .= 2*real(V'*convolve_fo_tr(res,V))
    return f
end

function crit_obj_fg(x,g,stars,psfs,data)
  T = eltype(x)
  nepochs = length(data)
  chi2_t = zeros(T, nepochs);
  singleepoch_g = [zeros(T, length(x)) for i=1:nepochs];
  Threads.@threads for i=1:nepochs # weighted sum -- should probably do the computation in parallel
    chi2_t[i] = crit_single_obj_fg(x, singleepoch_g[i], stars[i].polyft, stars[i].polyft', psfs[:,:,i], data[i]);
  end
  f = sum(chi2_t)
  g[:] .= sum(singleepoch_g)
  #println("Total chi2: $f");
  return f;
end

function crit_psf_fg(psfs,g_psf,stars,tmap_visible,data)
  T = eltype(psfs)
  nepochs = length(data)
  chi2_t = zeros(T, nepochs);
  singleepoch_g = [similar(psfs[:,:,1]) for i=1:nepochs];
  for i=1:nepochs # weighted sum -- should probably do the computation in parallel
    chi2_t[i] = crit_single_psf_fg(psfs[:,:,i], singleepoch_g[i], stars[i].polyft*tmap_visible, data[i]);
  end
  f = sum(chi2_t)
  g_psf[:] .= [singleepoch_g[i][1] for i=1:nepochs]
  return f;
end

function convolve(object::Array{T,2}, psf::SubArray{T, 2, Array{T, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}) where T
    return irfft(rfft(object).*rfft(psf),size(object,1))
end
      
function convolve_tr(object::Array{T,2}, psf::Array{T,2}) where T
    return irfft(conj(rfft(object)).*rfft(psf),size(object,1))  
end

# # function convolve_fo(psf::Array{T,2}, fo::Array{Complex{T},2}) where T
# #     return irfft(rfft(psf).*fo,size(psf,1))
# # end
   
# function convolve_fo_tr(psf::Array{T,2}, fo::Array{Complex{T},2}) where T
#     return irfft(rfft(psf).*conj.(fo),size(psf,1))  
# end

# function convolve_fo(psf::Array{T,2}, fo::Array{Complex{T},1}) where T
#   return irfft(reshape(vec(rfft(psf)).*fo, nx÷2+1,nx), size(psf,1))
# end

# function convolve_fo_tr(psf::Array{T,2}, fo::Array{Complex{T},1}) where T
#   return irfft(reshape(vec(rfft(psf)).*conj(fo), nx÷2+1,nx), size(psf,1))
# end

# function convolve_fo_adj(y::Array{T,2}, fo::Array{Complex{T},1}) where T
#     nx = size(y, 1)
#     fft_y = rfft(y)                                      # size: (nx, nx÷2+1)
#     result_fft = vec(fft_y) .* conj.(fo)                # element-wise product with conj(fo)
#     return irfft(reshape(result_fft, size(fft_y)), nx)  # reshape back and inverse FFT
# end


function convolve_fo(psf::Array{T,2}, fo::Array{Complex{T},1}) where T
    nx = size(psf, 1)
    fft_psf = rfft(psf)                                 # size: (nx, nx÷2+1)
    result_fft = reshape(vec(fft_psf) .* fo, size(fft_psf))  # multiply element-wise in vector form, reshape back
    return irfft(result_fft, nx)                        # return real result of convolution
end

function grad_x(x::Vector{T}, psf::Array{T,2}, H::Array{Complex{T},2}, d::Array{T,2}) where T
    nx = size(psf, 1)
    y = convolve_fo(psf, H * x)
    res = y - d                   
    return 2*real(H' * vec(conj.(rfft(psf)) .* rfft(res)))
end