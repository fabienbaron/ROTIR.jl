# =======================================================================================
# Tied parameters for orbit fitting
# =======================================================================================
# `orbit_fit_spec` already lets you choose which parameters are FREE. This adds the other
# half: parameters that are neither free nor fixed, but computed from the others.
#
# The motivating case is a disc lying in the orbital plane. Its projected major axis runs
# along the line of nodes, so its position angle is not an independent quantity — it is Ω.
# Fitting it freely both wastes a parameter and permits a solution in which the disc is
# oriented inconsistently with the orbit that was fitted alongside it:
#
#     spec = orbit_fit_spec(donor, disc; elements = el,
#                           free = [:a, :i, :Omega, :T0, :f, :c2_fwhm, :c2_ratio],
#                           ties = Dict(:c2_pa => "-Omega"))
#
# Expressions are ordinary Julia written against the parameter names in `spec.names`
# (`a, i, Omega, omega, e, P, T0, dP`, the component parameters prefixed `c1_`/`c2_`, and
# `f`). They are compiled ONCE, with RuntimeGeneratedFunctions, into a single call — a fit
# evaluates hundreds of thousands of likelihoods, and walking a dictionary of expressions at
# each one would cost more than the Kepler solve it is attached to.
#
# A tie may reference another tie; the evaluation order is resolved by topological sort, and
# a cycle is an error rather than a hang.
#
# CONVENTION WARNING. A tie is only as good as the convention it encodes. `EllipticalGaussian`
# takes `pa` as the angle whose axis is scaled by `ratio < 1`, i.e. the MINOR axis, measured
# from +RA toward +Dec — so the tie that puts the MAJOR axis on the line of nodes is
# `:c2_pa => "-Omega"`, not `"Omega"`. Check the component's own definition before writing
# one; a wrong tie is far more damaging than a free parameter, because it cannot be detected
# by looking at the fit quality alone.

using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

"""
    OrbitTies

Compiled tied-parameter rules for an [`OrbitFitSpec`](@ref).

`targets` indexes the tied entries of the parameter vector, and `fn(θ)` returns their values
as a tuple in that order. `exprs` keeps the original strings for display.

The type parameter `F` is load-bearing: with `fn::Function` (abstract) every call is a
dynamic dispatch and the returned tuple is boxed, which costs ~3 allocations per likelihood
evaluation. Concretely typed, resolving the ties allocates nothing beyond the parameter
vector copy.
"""
struct OrbitTies{F}
    targets::Vector{Int}
    names::Vector{Symbol}
    exprs::Vector{String}
    fn::F
end

Base.isempty(t::OrbitTies) = isempty(t.targets)

"Collect every symbol appearing in an expression (used to find parameter references)."
function _tie_symbols!(acc::Set{Symbol}, ex)
    if ex isa Symbol
        push!(acc, ex)
    elseif ex isa Expr
        # Skip the callee of a call: `sin` in `sin(Omega)` is a function, not a parameter.
        args = ex.head === :call ? ex.args[2:end] : ex.args
        for a in args; _tie_symbols!(acc, a); end
    end
    return acc
end

"""
    compile_ties(names, ties) -> OrbitTies

Compile `ties` (a `Symbol => expression-string` mapping) against the full parameter name
list `names`. Errors on an unknown target, an unknown reference, or a circular definition.
"""
function compile_ties(names::Vector{Symbol}, ties::AbstractDict)
    isempty(ties) && return OrbitTies(Int[], Symbol[], String[], _ -> ())
    nameset = Set(names)
    avail   = join(string.(names), ", ")

    targets = Symbol[]
    exprs   = Dict{Symbol,String}()
    for (k, v) in ties
        ks = Symbol(k)
        ks in nameset ||
            error("compile_ties: unknown tied parameter $ks; available: $avail")
        push!(targets, ks)
        exprs[ks] = String(v)
    end

    # Which tie depends on which other tie, so they can be evaluated in a valid order.
    parsed = Dict(k => Meta.parse(v) for (k, v) in exprs)
    refs   = Dict{Symbol,Vector{Symbol}}()
    for (k, ex) in parsed
        syms = _tie_symbols!(Set{Symbol}(), ex)
        for s in syms
            s == k && error("compile_ties: tie for $k refers to itself")
            (s in nameset || isdefined(@__MODULE__, s) || isdefined(Base, s)) ||
                error("compile_ties: tie for $k references unknown name `$s`; " *
                      "parameters are: $avail")
        end
        refs[k] = [s for s in syms if s in Set(targets) && s != k]
    end

    order  = Symbol[]
    indeg  = Dict(k => length(v) for (k, v) in refs)
    queue  = sort!([k for (k, d) in indeg if d == 0])
    while !isempty(queue)
        n = popfirst!(queue)
        push!(order, n)
        for (m, d) in refs
            if n in d
                indeg[m] -= 1
                indeg[m] == 0 && push!(queue, m)
            end
        end
    end
    length(order) == length(targets) ||
        error("compile_ties: circular tie among " *
              join(string.(setdiff(targets, order)), ", "))

    # Bind ONLY the parameters some tie actually reads. Binding all of them would leave the
    # compiled function carrying a dozen dead loads for a one-line tie; UltraNest calls this
    # hundreds of thousands of times, so the body is kept to exactly what is used.
    needed = Set{Symbol}()
    for ex in values(parsed); _tie_symbols!(needed, ex); end
    used = [i for i in eachindex(names) if names[i] in needed]
    binds = [:($(names[i]) = @inbounds θ[$i]) for i in used]
    tstmt = [:($k = $(parsed[k])) for k in order]
    # A TUPLE, not a vector: for the handful of ties a model realistically has this stays on
    # the stack, so resolving ties costs no allocation at all.
    body  = quote
        $(binds...)
        $(tstmt...)
        return $(Expr(:tuple, targets...))
    end
    fn = @RuntimeGeneratedFunction(:(θ -> $body))
    tidx = [findfirst(==(k), names)::Int for k in targets]
    return OrbitTies(tidx, targets, [exprs[k] for k in targets], fn)
end

"""
    apply_ties(t::OrbitTies, θ) -> θ

Return a copy of `θ` with every tied entry replaced by its computed value. `θ` is returned
unchanged (not copied) when there are no ties, so an untied fit pays nothing.
"""
function apply_ties(t::OrbitTies, θ)
    isempty(t.targets) && return θ
    out = collect(Float64, θ)
    vals = t.fn(out)
    @inbounds for (j, i) in enumerate(t.targets); out[i] = vals[j]; end
    return out
end

function Base.show(io::IO, t::OrbitTies)
    if isempty(t.targets)
        print(io, "OrbitTies (none)")
    else
        print(io, "OrbitTies with $(length(t.targets)) tie(s):")
        for (n, e) in zip(t.names, t.exprs); print(io, "\n  $n = $e"); end
    end
end
