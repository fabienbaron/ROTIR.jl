# Radial velocities against the orbit, PythonPlot backend. Included by ROTIRPythonPlotExt; the
# Makie version is `plot_rv_makie(bparams, rv::RVData)` in oiplot_spheroid_makie.jl. Both draw
# what `rv_curves` computes, so they cannot disagree about phases, offsets or residuals.

"""
    plot_rv(bparams, rv::RVData; K1, K2, γ = 0.0, offsets = Dict(), jitter = (a = 0, b = 0),
            residuals = true, xaxis = :phase, tref = nothing, figtitle) -> (fig, axes)

A table of radial velocities against the orbit, with the residuals `rv − model` below when
`residuals`. Each instrument has its own marker and is shifted by its offset; error bars are
√(σ² + s²) with `s` the star's jitter. `xaxis = :phase` or `:time`; see [`rv_curves`](@ref).
`axes` is `(ax,)` or `(ax, resax)`.
"""
function plot_rv(bparams, rv::RVData; K1::Real, K2::Real, γ::Real = 0.0,
                 offsets = Dict{String,Float64}(), jitter = (a = 0.0, b = 0.0),
                 residuals::Bool = true, xaxis::Symbol = :phase, tref = nothing,
                 figtitle = "Radial Velocities")
  set_oiplot_defaults()
  c = rv_curves(bparams, rv; K1 = K1, K2 = K2, gamma = γ, offsets = offsets,
                jitter = jitter, xaxis = xaxis, tref = tref)
  fig, axs = residuals ?
      subplots(2, 1, figsize = (10, 8), sharex = true,
               gridspec_kw = Dict("height_ratios" => [3, 1])) :
      subplots(1, 1, figsize = (10, 6))
  ax = residuals ? axs[0] : axs
  rax = residuals ? axs[1] : nothing
  ax.plot(c.xc, c.va, "b-", linewidth = 2, label = "Primary model")
  ax.plot(c.xc, c.vb, "r-", linewidth = 2, label = "Secondary model")
  ax.axhline(γ, color = "gray", linestyle = "--", linewidth = 0.8, alpha = 0.5)
  ins = rv_instruments(rv)
  marks = ("o", "s", "^", "D", "v", "*")
  for (s, col, name) in ((:a, "blue", "Primary"), (:b, "red", "Secondary")),
      (j, inst) in enumerate(ins)
    i = findall(k -> c.star[k] === s && c.instrument[k] == inst, eachindex(c.x))
    isempty(i) && continue
    mk = marks[mod1(j, length(marks))]
    lab = name * " data" * (isempty(inst) ? "" : " ($(inst))")
    ax.errorbar(c.x[i], c.y[i], yerr = c.sigma[i], fmt = mk, color = col, ms = 4,
                label = lab, zorder = 5)
    residuals && rax.errorbar(c.x[i], c.resid[i], yerr = c.sigma[i], fmt = mk, color = col,
                              ms = 4, zorder = 5)
  end
  ax.set_ylabel("Radial Velocity (km/s)")
  ax.set_title(figtitle)
  ax.legend()
  ax.grid(true, alpha = 0.3)
  xaxis === :phase && ax.set_xlim(0, 1)
  if residuals
    rax.axhline(0, color = "gray", linestyle = "--", linewidth = 0.8, alpha = 0.5)
    rax.set_ylabel("O − C (km/s)")
    rax.set_xlabel(c.xlabel)
    rax.grid(true, alpha = 0.3)
  else
    ax.set_xlabel(c.xlabel)
  end
  tight_layout()
  return fig, residuals ? (ax, rax) : (ax,)
end
