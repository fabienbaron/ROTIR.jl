mutable struct OIdata{T}
    # Complex visibilities
    visamp::Array{T,1}
    visamp_err::Array{T,1}
    visphi::Array{T,1}
    visphi_err::Array{T,1}
    vis_baseline::Array{T,1}
    vis_mjd::Array{T,1}
    vis_lam::Array{T,1}
    vis_dlam::Array{T,1}
    vis_flag::Array{Bool,1}
    # V2
    v2::Array{T,1}
    v2_err::Array{T,1}
    v2_baseline::Array{T,1}
    v2_mjd::Array{T,1}
    mean_mjd::T
    v2_lam::Array{T,1}
    v2_dlam::Array{T,1}
    v2_flag::Array{Bool,1}
    # T3
    t3amp::Array{T,1}
    t3amp_err::Array{T,1}
    t3phi::Array{T,1}
    t3phi_err::Array{T,1}
    t3phi_vonmises_err::Array{T,1}
    t3phi_vonmises_chi2_offset::Array{T,1}
    t3_baseline::Array{T,1}
    t3_maxbaseline::Array{T,1}
    t3_mjd::Array{T,1}
    t3_lam::Array{T,1}
    t3_dlam::Array{T,1}
    t3_flag::Array{Bool,1}
    #T4
#    t4amp::Array{T,1}
#    t4amp_err::Array{T,1}
#    t4phi::Array{T,1}
#    t4phi_err::Array{T,1}
#    t4_baseline::Array{T,1}
#    t4_maxbaseline::Array{T,1}
#    t4_mjd::Array{T,1}
#    t4_lam::Array{T,1}
#    t4_dlam::Array{T,1}
#    t4_flag::Array{Bool,1}
    #OIFlux
    flux::Array{T,1}
    flux_err::Array{T,1}
    flux_mjd::Array{T,1}
    flux_lam::Array{T,1}
    flux_dlam::Array{T,1}
    flux_flag::Array{Bool,1}
    flux_sta_index::Array{Int64,1}
    #UV coverage
    uv::Array{T,2}
    uv_lam::Array{T,1}
    uv_dlam::Array{T,1}
    uv_mjd::Array{T,1}
    uv_baseline::Array{T,1}
    # Data product sizes
    nflux::Int64
    nvisamp::Int64
    nvisphi::Int64
    nv2::Int64
    nt3amp::Int64
    nt3phi::Int64
    # nt4amp::Int64
    # nt4phi::Int64
    nuv::Int64
    # Indexing logic
    indx_vis::Array{Int64,1}
    indx_v2::Array{Int64,1}
    indx_t3_1::Array{Int64,1}
    indx_t3_2::Array{Int64,1}
    indx_t3_3::Array{Int64,1}
    # indx_t4_1::Array{Int64,1}
    # indx_t4_2::Array{Int64,1}
    # indx_t4_3::Array{Int64,1}
    # indx_t4_4::Array{Int64,1}
    sta_name::Array{String,1}
    tel_name::Array{String,1}
    sta_index::Array{Int64,1}
    vis_sta_index::Array{Int64,2}
    v2_sta_index::Array{Int64,2}
    t3_sta_index::Array{Int64,2}
    # t4_sta_index::Array{Int64,2}
    filename::String
end

function dataF32(data)
return OIdata{Float32}(
    Float32.(data.visamp),
    Float32.(data.visamp_err),
    Float32.(data.visphi),
    Float32.(data.visphi_err),
    Float32.(data.vis_baseline),
    Float32.(data.vis_mjd),
    Float32.(data.vis_lam),
    Float32.(data.vis_dlam),
    copy(data.vis_flag),
    Float32.(data.v2),
    Float32.(data.v2_err),
    Float32.(data.v2_baseline),
    Float32.(data.v2_mjd),
    Float32(data.mean_mjd),
    Float32.(data.v2_lam),
    Float32.(data.v2_dlam),
    copy(data.v2_flag),
    # T3
    Float32.(data.t3amp),
    Float32.(data.t3amp_err),
    Float32.(data.t3phi),
    Float32.(data.t3phi_err),
    Float32.(data.t3phi_vonmises_err),
    Float32.(data.t3phi_vonmises_chi2_offset),
    Float32.(data.t3_baseline),
    Float32.(data.t3_maxbaseline),
    Float32.(data.t3_mjd),
    Float32.(data.t3_lam),
    Float32.(data.t3_dlam),
    copy(data.t3_flag),
    #OIFlux
    Float32.(data.flux),
    Float32.(data.flux_err),
    Float32.(data.flux_mjd),
    Float32.(data.flux_lam),
    Float32.(data.flux_dlam),
    copy(data.flux_flag),
    copy(data.flux_sta_index),
    copy(data.uv),
    Float32.(data.uv_lam),
    Float32.(data.uv_dlam),
    Float32.(data.uv_mjd),
    Float32.(data.uv_baseline),
    data.nflux,
    data.nvisamp,
    data.nvisphi,
    data.nv2,
    data.nt3amp,
    data.nt3phi,
    data.nuv,
    copy(data.indx_vis),
    copy(data.indx_v2),
    copy(data.indx_t3_1),
    copy(data.indx_t3_2),
    copy(data.indx_t3_3),
    copy(data.sta_name),
    copy(data.tel_name),
    copy(data.sta_index),
    copy(data.vis_sta_index),
    copy(data.v2_sta_index),
    copy(data.t3_sta_index),
    data.filename
)
end