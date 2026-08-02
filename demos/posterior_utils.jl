# Shared plumbing for the sampler comparison: saving samples, corner plots, timing.
#
# Samples are written as plain CSV with `# key = value` metadata lines on top, so the
# comparison script needs no serialisation dependency and the files stay greppable.
#
# Corner plots go through the Python `corner` package via PyCall (it arrives with
# ultranest; otherwise `~/.julia/conda/3/x86_64/bin/pip install corner`). Plotting is
# always optional — a missing corner never costs you a run that took hours.

using DelimitedFiles, Printf, Statistics, PyCall
using Logging, TerminalLoggers, ProgressLogging

const RESULTS_DIR = get(ENV, "RESULTS_DIR", joinpath(@__DIR__, "results"))

# Pathfinder (and anything else built on Transducers/ProgressLogging) emits progress as
# log records rather than printing. A plain `julia script.jl` installs no sink for them, so
# they are silently discarded and the run looks frozen for minutes at a time. Installing a
# TerminalLogger makes those bars appear. Skipped when stdout is not a terminal — in a
# redirected log the escape codes are worse than no bar. PROGRESS=0 disables it.
# PROGRESS=auto (default) installs it only on a terminal; PROGRESS=1 forces it on (useful
# under tmux/nohup where stdout is a pipe); PROGRESS=0 disables it.
function install_progress_logger()
    mode = get(ENV, "PROGRESS", "auto")
    mode == "0" && return nothing
    (mode == "1" || isa(stdout, Base.TTY)) || return nothing
    global_logger(TerminalLogger(stderr, Logging.Info))
    return nothing
end
install_progress_logger()

"""
    progress_counter(total; every=max(total ÷ 20, 1), label="")

Thread-safe counter for hand-written parallel loops: call the returned closure once per
completed item. Prints `label i/total (pct%)` on one line.
"""
function progress_counter(total::Integer; every::Integer = max(total ÷ 20, 1), label = "")
    done = Threads.Atomic{Int}(0)
    t0 = time()
    return function tick()
        n = Threads.atomic_add!(done, 1) + 1
        if n % every == 0 || n == total
            el = time() - t0
            @printf("\r%s%d/%d (%.0f%%, %.0f s elapsed, ~%.0f s left)   ",
                    isempty(label) ? "" : label * " ", n, total, 100n / total, el,
                    el * (total - n) / max(n, 1))
            flush(stdout)
            n == total && println()
        end
        return n
    end
end

"""
    save_posterior(method, samples, labels; meta...)

Write `samples` (nsamples × nparams, physical units) to `results/<method>.csv`.
"""
function save_posterior(method::AbstractString, samples::AbstractMatrix,
                        labels::AbstractVector; meta...)
    mkpath(RESULTS_DIR)
    path = joinpath(RESULTS_DIR, method * ".csv")
    open(path, "w") do io
        println(io, "# method = ", method)
        println(io, "# nsamples = ", size(samples, 1))
        for (k, v) in pairs(meta)
            println(io, "# ", k, " = ", v)
        end
        println(io, join(labels, ","))
        writedlm(io, samples, ',')
    end
    @printf("wrote %s (%d samples)\n", path, size(samples, 1))
    return path
end

"""
    load_posterior(path) -> (samples, labels, meta)
"""
function load_posterior(path::AbstractString)
    meta = Dict{String,String}()
    labels = String[]
    rows = Vector{Vector{Float64}}()
    for line in eachline(path)
        if startswith(line, "#")
            kv = split(chopprefix(line, "#"), "=", limit = 2)
            length(kv) == 2 && (meta[strip(kv[1])] = strip(kv[2]))
        elseif isempty(labels)
            labels = String.(strip.(split(line, ",")))
        else
            push!(rows, parse.(Float64, split(line, ",")))
        end
    end
    return reduce(vcat, transpose.(rows)), labels, meta
end

"""
    summarise(samples, labels; title="")

Median and 16/84 percentile half-widths, printed and returned as (median, sigma).
"""
function summarise(samples::AbstractMatrix, labels::AbstractVector; title = "")
    med = [median(view(samples, :, j)) for j in eachindex(labels)]
    lo  = [quantile(view(samples, :, j), 0.16) for j in eachindex(labels)]
    hi  = [quantile(view(samples, :, j), 0.84) for j in eachindex(labels)]
    σ   = 0.5 .* (hi .- lo)
    isempty(title) || println(title)
    for j in eachindex(labels)
        @printf("  %-8s %10.5f  +%.5f -%.5f\n", labels[j], med[j], hi[j] - med[j],
                med[j] - lo[j])
    end
    return med, σ
end

"""
    corner_plot(samples, labels, filename; truths=nothing, title="")

Corner plot via Python `corner`. Returns the path, or `nothing` if corner is unavailable —
never throws, because these scripts run for hours and a plotting dependency is not worth
losing a run over.
"""
function corner_plot(samples::AbstractMatrix, labels::AbstractVector,
                     filename::AbstractString; truths = nothing, title = "")
    mkpath(RESULTS_DIR)
    path = joinpath(RESULTS_DIR, filename)
    try
        corner = pyimport("corner")
        plt = pyimport("matplotlib.pyplot")
        fig = corner.corner(Array{Float64}(samples); labels = collect(String, labels),
                            quantiles = [0.16, 0.5, 0.84], show_titles = true,
                            title_fmt = ".4f",
                            truths = truths === nothing ? nothing : collect(Float64, truths))
        isempty(title) || fig.suptitle(title)
        fig.savefig(path, dpi = 130, bbox_inches = "tight")
        plt.close(fig)
        println("wrote ", path)
        return path
    catch e
        @warn "corner plot skipped (pip install corner)" exception = e
        return nothing
    end
end

"Run `f`, returning `(result, wall_seconds)`."
function timed(f)
    t0 = time()
    r = f()
    return r, time() - t0
end
