# rasterize.jl — Fast polygon rasterization for stellar surface imaging
#
# Instead of computing the continuous polygon Fourier transform at every UV
# point (expensive: O(Npix x Nk x 4) trig ops), rasterize the polygons
# directly onto a real-space image grid via exact pixel-polygon intersection
# using Sutherland-Hodgman clipping.  The resulting image is the
# pixel-integrated polygon model; one FFT gives the Fourier-space model if
# needed.
#
# Adapted from planet_deconv for use with ROTIR's coordinate conventions
# (proj_west, proj_north in mas).
#
# This file only adds NEW functions -- it does not modify existing polyft code.

# --- Sutherland-Hodgman clip against each half-plane -------------------------
# Each clip reads (in_x, in_y, n) and writes to (out_x, out_y), returning the
# new vertex count.  Works for convex or concave polygons; output has at most
# n+1 vertices per clip.

@inline function _clip_left!(in_x, in_y, n::Int, xmin::T, out_x, out_y) where T
    m = 0
    @inbounds for i in 1:n
        j = i == n ? 1 : i + 1
        xi = in_x[i]; yi = in_y[i]
        xj = in_x[j]; yj = in_y[j]
        in_i = xi >= xmin
        in_j = xj >= xmin
        if in_j
            if !in_i
                t = (xmin - xi) / (xj - xi)
                m += 1; out_x[m] = xmin; out_y[m] = yi + t*(yj - yi)
            end
            m += 1; out_x[m] = xj; out_y[m] = yj
        elseif in_i
            t = (xmin - xi) / (xj - xi)
            m += 1; out_x[m] = xmin; out_y[m] = yi + t*(yj - yi)
        end
    end
    return m
end

@inline function _clip_right!(in_x, in_y, n::Int, xmax::T, out_x, out_y) where T
    m = 0
    @inbounds for i in 1:n
        j = i == n ? 1 : i + 1
        xi = in_x[i]; yi = in_y[i]
        xj = in_x[j]; yj = in_y[j]
        in_i = xi <= xmax
        in_j = xj <= xmax
        if in_j
            if !in_i
                t = (xmax - xi) / (xj - xi)
                m += 1; out_x[m] = xmax; out_y[m] = yi + t*(yj - yi)
            end
            m += 1; out_x[m] = xj; out_y[m] = yj
        elseif in_i
            t = (xmax - xi) / (xj - xi)
            m += 1; out_x[m] = xmax; out_y[m] = yi + t*(yj - yi)
        end
    end
    return m
end

@inline function _clip_bottom!(in_x, in_y, n::Int, ymin::T, out_x, out_y) where T
    m = 0
    @inbounds for i in 1:n
        j = i == n ? 1 : i + 1
        xi = in_x[i]; yi = in_y[i]
        xj = in_x[j]; yj = in_y[j]
        in_i = yi >= ymin
        in_j = yj >= ymin
        if in_j
            if !in_i
                t = (ymin - yi) / (yj - yi)
                m += 1; out_x[m] = xi + t*(xj - xi); out_y[m] = ymin
            end
            m += 1; out_x[m] = xj; out_y[m] = yj
        elseif in_i
            t = (ymin - yi) / (yj - yi)
            m += 1; out_x[m] = xi + t*(xj - xi); out_y[m] = ymin
        end
    end
    return m
end

@inline function _clip_top!(in_x, in_y, n::Int, ymax::T, out_x, out_y) where T
    m = 0
    @inbounds for i in 1:n
        j = i == n ? 1 : i + 1
        xi = in_x[i]; yi = in_y[i]
        xj = in_x[j]; yj = in_y[j]
        in_i = yi <= ymax
        in_j = yj <= ymax
        if in_j
            if !in_i
                t = (ymax - yi) / (yj - yi)
                m += 1; out_x[m] = xi + t*(xj - xi); out_y[m] = ymax
            end
            m += 1; out_x[m] = xj; out_y[m] = yj
        elseif in_i
            t = (ymax - yi) / (yj - yi)
            m += 1; out_x[m] = xi + t*(xj - xi); out_y[m] = ymax
        end
    end
    return m
end

"""
    quad_box_area(q1x,q1y,q2x,q2y,q3x,q3y,q4x,q4y, xmin,xmax,ymin,ymax, Ax,Ay,Bx,By)

Exact area of intersection between a (convex or concave) quadrilateral and an
axis-aligned box via Sutherland-Hodgman clipping + shoelace.  Scratch buffers
`Ax,Ay,Bx,By` must each hold at least 12 elements.
"""
@inline function quad_box_area(q1x::T, q1y::T, q2x::T, q2y::T,
                                q3x::T, q3y::T, q4x::T, q4y::T,
                                xmin::T, xmax::T, ymin::T, ymax::T,
                                Ax, Ay, Bx, By) where T
    @inbounds begin
        Ax[1] = q1x; Ay[1] = q1y
        Ax[2] = q2x; Ay[2] = q2y
        Ax[3] = q3x; Ay[3] = q3y
        Ax[4] = q4x; Ay[4] = q4y
    end
    n = 4
    n = _clip_left!(Ax, Ay, n, xmin, Bx, By);  n == 0 && return zero(T)
    n = _clip_right!(Bx, By, n, xmax, Ax, Ay); n == 0 && return zero(T)
    n = _clip_bottom!(Ax, Ay, n, ymin, Bx, By); n == 0 && return zero(T)
    n = _clip_top!(Bx, By, n, ymax, Ax, Ay);   n < 3  && return zero(T)
    area2 = zero(T)
    @inbounds for i in 1:n
        j = i == n ? 1 : i + 1
        area2 += Ax[i] * Ay[j] - Ax[j] * Ay[i]
    end
    return abs(area2) * T(0.5)
end

"""
    rasterize_polygon_image!(img, proj_west, proj_north, x_weighted, pixsize;
                             cx=nx/2+1, cy=ny/2+1, w_threshold=1e-10)

Build a centered real-space image by rasterizing each quadrilateral tile with
intensity `x_weighted[p]`.  Output is in the "data-centered" convention:
`proj_west=proj_north=0` maps to pixel `(cx, cy)`.

# Arguments
- `img`: `(ny, nx)` preallocated output matrix (will be zeroed and filled)
- `proj_west, proj_north`: `(Npix, 4)` polygon vertex coordinates (mas)
- `x_weighted`: `(Npix,)` intensity weights (e.g. `xmap .* vis_weights`)
- `pixsize`: mas per pixel
- `cx, cy`: pixel index (1-based) where the projection origin maps

Each pixel accumulates `sum_p x_weighted[p] * Area(tile_p intersect pixel)`.
Threaded over polygons (per-thread local images, reduced at the end).
"""
function rasterize_polygon_image!(img::Matrix{T}, proj_west::AbstractMatrix,
                                    proj_north::AbstractMatrix,
                                    x_weighted::AbstractVector,
                                    pixsize::Real;
                                    cx::Real = size(img,2) ÷ 2 + 1,
                                    cy::Real = size(img,1) ÷ 2 + 1,
                                    w_threshold::Real = 1e-10) where T
    ny, nx_loc = size(img)
    Npix = size(proj_west, 1)
    fill!(img, zero(T))

    invpx = T(1 / pixsize)
    cx_T = T(cx)
    cy_T = T(cy)
    wth  = T(w_threshold)
    px2  = T(pixsize) * T(pixsize)

    nt = Threads.maxthreadid()
    local_imgs = [zeros(T, ny, nx_loc) for _ in 1:nt]

    Threads.@threads for p in 1:Npix
        xw = T(x_weighted[p]) * px2
        if abs(xw) < wth * px2
            continue
        end

        q1x = T(proj_west[p,1]) * invpx + cx_T
        q2x = T(proj_west[p,2]) * invpx + cx_T
        q3x = T(proj_west[p,3]) * invpx + cx_T
        q4x = T(proj_west[p,4]) * invpx + cx_T
        q1y = T(proj_north[p,1]) * invpx + cy_T
        q2y = T(proj_north[p,2]) * invpx + cy_T
        q3y = T(proj_north[p,3]) * invpx + cy_T
        q4y = T(proj_north[p,4]) * invpx + cy_T

        xmn_f = min(q1x, q2x, q3x, q4x)
        xmx_f = max(q1x, q2x, q3x, q4x)
        ymn_f = min(q1y, q2y, q3y, q4y)
        ymx_f = max(q1y, q2y, q3y, q4y)

        ixmn = max(1, floor(Int, xmn_f + T(0.5)))
        ixmx = min(nx_loc, floor(Int, xmx_f + T(0.5)) + 1)
        iymn = max(1, floor(Int, ymn_f + T(0.5)))
        iymx = min(ny, floor(Int, ymx_f + T(0.5)) + 1)
        if ixmn > ixmx || iymn > iymx
            continue
        end

        tid = Threads.threadid()
        limg = local_imgs[tid]

        Ax = Vector{T}(undef, 12); Ay = Vector{T}(undef, 12)
        Bx = Vector{T}(undef, 12); By = Vector{T}(undef, 12)

        @inbounds for iy in iymn:iymx, ix in ixmn:ixmx
            pxmin = T(ix) - T(0.5)
            pxmax = T(ix) + T(0.5)
            pymin = T(iy) - T(0.5)
            pymax = T(iy) + T(0.5)
            area = quad_box_area(q1x,q1y, q2x,q2y, q3x,q3y, q4x,q4y,
                                  pxmin,pxmax,pymin,pymax,
                                  Ax,Ay,Bx,By)
            if area > zero(T)
                limg[iy, ix] += xw * area
            end
        end
    end

    @inbounds for t in 1:nt
        img .+= local_imgs[t]
    end
    return img
end

"""
    rasterize_polygon_image(proj_west, proj_north, x_weighted, pixsize, nx;
                            cx=nx/2+1, cy=nx/2+1, T=Float32)

Non-mutating convenience wrapper: allocates and returns an `(nx, nx)` image.

# Arguments
- `proj_west, proj_north`: `(Npix, 4)` polygon vertex coordinates (mas)
- `x_weighted`: `(Npix,)` intensity weights
- `pixsize`: mas per pixel
- `nx`: image dimension (pixels)
"""
function rasterize_polygon_image(proj_west, proj_north, x_weighted, pixsize, nx;
                                  cx::Real=nx ÷ 2 + 1, cy::Real=nx ÷ 2 + 1,
                                  T::Type=Float32)
    img = zeros(T, nx, nx)
    rasterize_polygon_image!(img, proj_west, proj_north, x_weighted, pixsize; cx=cx, cy=cy)
    return img
end

"""
    rasterize_adjoint!(grad, res_img, proj_west, proj_north, pixel_weight, pixsize;
                       cx=nx/2+1, cy=ny/2+1)

Adjoint (transpose) of `rasterize_polygon_image!` for gradient computation.
Given a residual image `res_img` (same layout as the rasterized output),
compute `grad[p] = pixel_weight[p] * pixsize^2 * sum_{iy,ix} res_img[iy,ix] * area(tile_p intersect pixel_{iy,ix})`.

This is the exact transpose of the forward rasterization: the forward spreads
polygon intensities into pixels weighted by overlap area; the adjoint gathers
pixel values back to polygons weighted by the same overlap area.

Threaded over polygons (per-thread local gradient, reduced at the end).
"""
function rasterize_adjoint!(grad::AbstractVector{T}, res_img::Matrix{T},
                             proj_west::AbstractMatrix, proj_north::AbstractMatrix,
                             pixel_weight::AbstractVector,
                             pixsize::Real;
                             cx::Real = size(res_img,2) ÷ 2 + 1,
                             cy::Real = size(res_img,1) ÷ 2 + 1) where T
    ny, nx_loc = size(res_img)
    Npix = size(proj_west, 1)
    fill!(grad, zero(T))

    invpx = T(1 / pixsize)
    cx_T = T(cx); cy_T = T(cy)
    px2  = T(pixsize) * T(pixsize)

    nt = Threads.maxthreadid()
    local_grads = [zeros(T, Npix) for _ in 1:nt]

    Threads.@threads for p in 1:Npix
        pw = T(pixel_weight[p]) * px2
        if abs(pw) < T(1e-10) * px2
            continue
        end

        q1x = T(proj_west[p,1]) * invpx + cx_T
        q2x = T(proj_west[p,2]) * invpx + cx_T
        q3x = T(proj_west[p,3]) * invpx + cx_T
        q4x = T(proj_west[p,4]) * invpx + cx_T
        q1y = T(proj_north[p,1]) * invpx + cy_T
        q2y = T(proj_north[p,2]) * invpx + cy_T
        q3y = T(proj_north[p,3]) * invpx + cy_T
        q4y = T(proj_north[p,4]) * invpx + cy_T

        xmn_f = min(q1x, q2x, q3x, q4x)
        xmx_f = max(q1x, q2x, q3x, q4x)
        ymn_f = min(q1y, q2y, q3y, q4y)
        ymx_f = max(q1y, q2y, q3y, q4y)

        ixmn = max(1, floor(Int, xmn_f + T(0.5)))
        ixmx = min(nx_loc, floor(Int, xmx_f + T(0.5)) + 1)
        iymn = max(1, floor(Int, ymn_f + T(0.5)))
        iymx = min(ny, floor(Int, ymx_f + T(0.5)) + 1)
        if ixmn > ixmx || iymn > iymx
            continue
        end

        tid = Threads.threadid()
        lg = local_grads[tid]

        Ax = Vector{T}(undef, 12); Ay = Vector{T}(undef, 12)
        Bx = Vector{T}(undef, 12); By = Vector{T}(undef, 12)

        acc = zero(T)
        @inbounds for iy in iymn:iymx, ix in ixmn:ixmx
            pxmin = T(ix) - T(0.5)
            pxmax = T(ix) + T(0.5)
            pymin = T(iy) - T(0.5)
            pymax = T(iy) + T(0.5)
            area = quad_box_area(q1x,q1y, q2x,q2y, q3x,q3y, q4x,q4y,
                                  pxmin,pxmax,pymin,pymax,
                                  Ax,Ay,Bx,By)
            if area > zero(T)
                acc += res_img[iy, ix] * area
            end
        end
        lg[p] = pw * acc
    end

    @inbounds for t in 1:nt
        grad .+= local_grads[t]
    end
    return grad
end
