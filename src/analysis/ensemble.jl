"""
    final_size(el::EventLog)

Number of seeded plus transmitted hosts in an event log.
"""
final_size(el::EventLog) = count(==(EK_Seeding), el.kind) + count(==(EK_Transmission), el.kind)


"""
    event_count(el::EventLog, kind)

Count events of a single kind in an event log.

`kind` may be an [`EventKind`](@ref), a user-facing symbol such as
`:transmission`, `:activation`, `:removal`, `:fossilized_sampling`, or
`:serial_sampling`, or the corresponding string.
"""
event_count(el::EventLog, kind::EventKind) = count(==(kind), el.kind)
event_count(el::EventLog, kind::Symbol) = event_count(el, event_kind(kind))
event_count(el::EventLog, kind::AbstractString) = event_count(el, event_kind(kind))


"""
    final_time(el::EventLog)

Time of the last event in an event log, or `0.0` for an empty log.
"""
final_time(el::EventLog) = isempty(el.time) ? 0.0 : el.time[end]


"""
    EnsembleSummary

Compact per-replicate summary returned by [`run_ensemble`](@ref).

The summary keeps primitive vectors for common outbreak-level statistics and
stores raw logs only when requested. This keeps the default representation
small enough for large ensembles while still allowing derived views to use full
event logs when needed.
"""
struct EnsembleSummary
    nrep::Int
    final_size::Vector{Int}
    total_events::Vector{Int}
    final_time::Vector{Float64}
    transmissions::Vector{Int}
    activations::Vector{Int}
    removals::Vector{Int}
    fossilized_samples::Vector{Int}
    serial_samples::Vector{Int}
    logs::Union{Nothing,Vector{EventLog}}
end

"""
    EnsembleReplicateSummary

One replicate row from an [`EnsembleSummary`](@ref).
"""
struct EnsembleReplicateSummary
    replicate::Int
    final_size::Int
    total_events::Int
    final_time::Float64
    transmissions::Int
    activations::Int
    removals::Int
    fossilized_samples::Int
    serial_samples::Int
    log::Union{Nothing,EventLog}
end


Base.length(summary::EnsembleSummary) = summary.nrep
Base.isempty(summary::EnsembleSummary) = length(summary) == 0
Base.firstindex(::EnsembleSummary) = 1
Base.lastindex(summary::EnsembleSummary) = length(summary)
Base.eltype(::Type{EnsembleSummary}) = EnsembleReplicateSummary


function Base.getindex(summary::EnsembleSummary, i::Integer)
    return EnsembleReplicateSummary(
        Int(i),
        summary.final_size[i],
        summary.total_events[i],
        summary.final_time[i],
        summary.transmissions[i],
        summary.activations[i],
        summary.removals[i],
        summary.fossilized_samples[i],
        summary.serial_samples[i],
        summary.logs === nothing ? nothing : summary.logs[i],
    )
end


function Base.iterate(summary::EnsembleSummary, state::Int=firstindex(summary))
    state > lastindex(summary) && return nothing
    return (summary[state], state + 1)
end


function Base.show(io::IO, summary::EnsembleSummary)
    retained = summary.logs === nothing ? "not retained" : "retained"
    print(io, "EnsembleSummary(", summary.nrep, " replicates")
    print(io, ", logs ", retained, ")")
end


function Base.show(io::IO, row::EnsembleReplicateSummary)
    retained = row.log === nothing ? "log not retained" : "log retained"
    print(io, "EnsembleReplicateSummary(rep=", row.replicate,
          ", final_size=", row.final_size,
          ", total_events=", row.total_events,
          ", final_time=", row.final_time,
          ", transmissions=", row.transmissions,
          ", activations=", row.activations,
          ", removals=", row.removals,
          ", fossilized_samples=", row.fossilized_samples,
          ", serial_samples=", row.serial_samples,
          ", ", retained,
          ")")
end


"""
    mean_final_size(summary::EnsembleSummary)

Mean final outbreak size across ensemble replicates.
"""
mean_final_size(summary::EnsembleSummary) = _mean(summary.final_size)


"""
    mean_final_time(summary::EnsembleSummary)

Mean final event time across ensemble replicates.
"""
mean_final_time(summary::EnsembleSummary) = _mean(summary.final_time)


"""
    attack_rate(summary::EnsembleSummary, population_size)

Mean final outbreak size divided by a finite population size.
"""
function attack_rate(summary::EnsembleSummary, population_size::Integer)
    population_size > 0 || throw(ArgumentError("population_size must be positive"))
    return mean_final_size(summary) / population_size
end


_mean(xs::AbstractVector) = isempty(xs) ? NaN : sum(xs) / length(xs)


function _summarize_log!(summary::EnsembleSummary, rep::Int, el::EventLog)
    summary.final_size[rep] = final_size(el)
    summary.total_events[rep] = length(el)
    summary.final_time[rep] = final_time(el)
    summary.transmissions[rep] = event_count(el, EK_Transmission)
    summary.activations[rep] = event_count(el, EK_Activation)
    summary.removals[rep] = event_count(el, EK_Removal)
    summary.fossilized_samples[rep] = event_count(el, EK_FossilizedSampling)
    summary.serial_samples[rep] = event_count(el, EK_SerialSampling)
    return summary
end


function _empty_ensemble_summary(nrep::Int, retain_logs::Bool)
    nrep >= 0 || throw(ArgumentError("nrep must be non-negative"))
    return EnsembleSummary(
        nrep,
        Vector{Int}(undef, nrep),
        Vector{Int}(undef, nrep),
        Vector{Float64}(undef, nrep),
        Vector{Int}(undef, nrep),
        Vector{Int}(undef, nrep),
        Vector{Int}(undef, nrep),
        Vector{Int}(undef, nrep),
        Vector{Int}(undef, nrep),
        retain_logs ? Vector{EventLog}(undef, nrep) : nothing,
    )
end


"""
    run_ensemble(simulator, nrep; rng=Random.default_rng(), retain_logs=false, validate=true)

Run `nrep` independent replicates from a `simulator(rng)` returning
an [`EventLog`](@ref), and return a compact [`EnsembleSummary`](@ref).

Each replicate receives a fresh child RNG seeded from `rng`. The ensemble layer
never reseeds Julia's global RNG. Set `retain_logs=true` to keep raw event logs
for later per-log derived views such as trajectories or host summaries; by
default only compact per-replicate statistics are retained.
"""
function run_ensemble(simulator, nrep::Integer; rng::AbstractRNG=Random.default_rng(), retain_logs::Bool=false, validate::Bool=true)
    n = Int(nrep)
    summary = _empty_ensemble_summary(n, retain_logs)

    for rep in 1:n
        rep_rng = Random.MersenneTwister(rand(rng, UInt))
        el = simulator(rep_rng)
        el isa EventLog || throw(ArgumentError("simulator must return an EventLog"))
        validate && validate_event_log(el)
        _summarize_log!(summary, rep, el)
        retain_logs && (summary.logs[rep] = el)
    end

    return summary
end
