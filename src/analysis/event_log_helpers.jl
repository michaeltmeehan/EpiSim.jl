"""
    event_indices(el::EventLog, kind)

Return event-log row indices whose event kind is `kind`.

`kind` may be an [`EventKind`](@ref), a user-facing symbol such as
`:transmission`, or the corresponding string.
"""
function event_indices(el::EventLog, kind::EventKind)
    indices = Int[]
    @inbounds for i in eachindex(el.kind)
        el.kind[i] == kind && push!(indices, i)
    end
    return indices
end
event_indices(el::EventLog, kind::Symbol) = event_indices(el, event_kind(kind))
event_indices(el::EventLog, kind::AbstractString) = event_indices(el, event_kind(kind))


"""
    event_times(el::EventLog, kind)

Return event times for events whose kind is `kind`, preserving log order.

`kind` may be an [`EventKind`](@ref), a user-facing symbol such as
`:transmission`, or the corresponding string.
"""
function event_times(el::EventLog, kind::EventKind)
    times = Float64[]
    @inbounds for i in eachindex(el.kind)
        el.kind[i] == kind && push!(times, el.time[i])
    end
    return times
end
event_times(el::EventLog, kind::Symbol) = event_times(el, event_kind(kind))
event_times(el::EventLog, kind::AbstractString) = event_times(el, event_kind(kind))


"""
    first_event_time(el::EventLog, kind)

Return the first time for event kind `kind`, or `nothing` if absent.
"""
function first_event_time(el::EventLog, kind::EventKind)
    @inbounds for i in eachindex(el.kind)
        el.kind[i] == kind && return el.time[i]
    end
    return nothing
end
first_event_time(el::EventLog, kind::Symbol) = first_event_time(el, event_kind(kind))
first_event_time(el::EventLog, kind::AbstractString) = first_event_time(el, event_kind(kind))


"""
    last_event_time(el::EventLog, kind)

Return the last time for event kind `kind`, or `nothing` if absent.
"""
function last_event_time(el::EventLog, kind::EventKind)
    @inbounds for i in lastindex(el.kind):-1:firstindex(el.kind)
        el.kind[i] == kind && return el.time[i]
    end
    return nothing
end
last_event_time(el::EventLog, kind::Symbol) = last_event_time(el, event_kind(kind))
last_event_time(el::EventLog, kind::AbstractString) = last_event_time(el, event_kind(kind))


"""
    has_event_kind(el::EventLog, kind)

Return whether `el` contains at least one event of `kind`.
"""
has_event_kind(el::EventLog, kind::EventKind) = first_event_time(el, kind) !== nothing
has_event_kind(el::EventLog, kind::Symbol) = has_event_kind(el, event_kind(kind))
has_event_kind(el::EventLog, kind::AbstractString) = has_event_kind(el, event_kind(kind))


"""
    total_events(el::EventLog)

Return the number of events in an event log.
"""
total_events(el::EventLog) = length(el)


"""
    observed_hosts(el::EventLog; include_sources=true)

Return sorted semantically valid host IDs observed in `host`, and optionally
positive `infector` source IDs. Sentinel values such as `0` and negative
placeholders are never returned.
"""
function observed_hosts(el::EventLog; include_sources::Bool=true)
    hosts = Set{Int}()
    @inbounds for i in eachindex(el.host)
        h = el.host[i]
        h > 0 && push!(hosts, h)
        if include_sources
            source = el.infector[i]
            source > 0 && push!(hosts, source)
        end
    end
    return sort!(collect(hosts))
end


"""
    distinct_host_count(el::EventLog; include_sources=true)

Count distinct semantically valid positive host IDs observed in `host`, and
optionally positive `infector` source IDs.
"""
distinct_host_count(el::EventLog; include_sources::Bool=true) =
    length(observed_hosts(el; include_sources=include_sources))


"""
    HostEventSummary

Compact per-observed-host event participation summary for one [`EventLog`](@ref).

`host_id` is sorted ascending. All count vectors are aligned to `host_id`.
`transmissions_caused` counts appearances as `infector` on transmission events.
Sampling counts include both fossilized and serial sampling, while `removals`
counts `EK_Removal` plus `EK_SerialSampling` because serial sampling removes the
host from the infectious state.
"""
struct HostEventSummary
    host_id::Vector{Int}
    transmissions_caused::Vector{Int}
    samples::Vector{Int}
    removals::Vector{Int}
    activations::Vector{Int}
end

"""
    HostEventRecord

One row from a [`HostEventSummary`](@ref), indexed by sorted summary position.
"""
struct HostEventRecord
    host_id::Int
    transmissions_caused::Int
    samples::Int
    removals::Int
    activations::Int
end


Base.length(summary::HostEventSummary) = length(summary.host_id)
Base.isempty(summary::HostEventSummary) = isempty(summary.host_id)
Base.firstindex(::HostEventSummary) = 1
Base.lastindex(summary::HostEventSummary) = length(summary)
Base.eltype(::Type{HostEventSummary}) = HostEventRecord


function Base.getindex(summary::HostEventSummary, i::Integer)
    return HostEventRecord(
        summary.host_id[i],
        summary.transmissions_caused[i],
        summary.samples[i],
        summary.removals[i],
        summary.activations[i],
    )
end


function Base.iterate(summary::HostEventSummary, state::Int=firstindex(summary))
    state > lastindex(summary) && return nothing
    return (summary[state], state + 1)
end


function Base.show(io::IO, summary::HostEventSummary)
    print(io, "HostEventSummary(", length(summary), " observed hosts")
    if isempty(summary)
        print(io, ", empty")
    else
        print(io, ", host id range ", first(summary.host_id), " to ", last(summary.host_id))
    end
    print(io, ")")
end


function Base.show(io::IO, record::HostEventRecord)
    print(io, "HostEventRecord(host_id=", record.host_id,
          ", transmissions_caused=", record.transmissions_caused,
          ", samples=", record.samples,
          ", removals=", record.removals,
          ", activations=", record.activations,
          ")")
end


"""
    host_event_summary(el::EventLog)

Summarize host participation for one event log over observed host IDs.

The returned [`HostEventSummary`](@ref) is indexed by sorted observed host ID,
not by dense population index. Observed hosts include positive IDs from both
`host` and `infector` columns. Non-positive sentinel values are ignored.
"""
function host_event_summary(el::EventLog)
    ids = observed_hosts(el; include_sources=true)
    positions = Dict{Int,Int}(id => i for (i, id) in pairs(ids))
    n = length(ids)

    transmissions_caused = zeros(Int, n)
    samples = zeros(Int, n)
    removals = zeros(Int, n)
    activations = zeros(Int, n)

    @inbounds for i in eachindex(el.kind)
        kind = el.kind[i]
        host = el.host[i]

        if kind == EK_Transmission
            source = el.infector[i]
            source > 0 && haskey(positions, source) && (transmissions_caused[positions[source]] += 1)
        elseif kind == EK_Activation
            host > 0 && haskey(positions, host) && (activations[positions[host]] += 1)
        elseif kind == EK_Removal
            host > 0 && haskey(positions, host) && (removals[positions[host]] += 1)
        elseif kind == EK_FossilizedSampling
            host > 0 && haskey(positions, host) && (samples[positions[host]] += 1)
        elseif kind == EK_SerialSampling
            if host > 0 && haskey(positions, host)
                samples[positions[host]] += 1
                removals[positions[host]] += 1
            end
        end
    end

    return HostEventSummary(ids, transmissions_caused, samples, removals, activations)
end
