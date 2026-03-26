# Soft visibility via sigmoid — replaces hard normals[:,3] > 0 mask
# Provides smooth, differentiable limb weighting for gradient-based optimization.

"""
    sigmoid(x)
Numerically stable sigmoid function.
"""
function sigmoid(x::T) where T<:AbstractFloat
    if x >= zero(T)
        z = exp(-x)
        return one(T) / (one(T) + z)
    else
        z = exp(x)
        return z / (one(T) + z)
    end
end

"""
    dsigmoid(x)
Derivative of sigmoid: σ'(x) = σ(x) * (1 - σ(x))
"""
function dsigmoid(x::T) where T<:AbstractFloat
    s = sigmoid(x)
    return s * (one(T) - s)
end

"""
    soft_visibility(normals_z; κ=50.0)
Compute soft visibility weights from the z-component of surface normals.
Returns `(weights, sigmoid_args)` where:
- `weights[i] = sigmoid(κ * normals_z[i])` ∈ (0, 1)
- `sigmoid_args[i] = κ * normals_z[i]` (stored for gradient computation)

κ controls sharpness: larger κ → approaches hard mask.
"""
function soft_visibility(normals_z::AbstractVector{T}; κ::T=T(50)) where T<:AbstractFloat
    npix = length(normals_z)
    weights = Vector{T}(undef, npix)
    sig_args = Vector{T}(undef, npix)
    @inbounds for i in 1:npix
        sig_args[i] = κ * normals_z[i]
        weights[i] = sigmoid(sig_args[i])
    end
    return weights, sig_args
end
