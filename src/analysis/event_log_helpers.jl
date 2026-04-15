"""
    event_indices(el::EventLog, kind::EventKind)

Return event-log row indices whose event kind is `kind`.
"""
function event_indices(el::EventLog, kind::EventKind)
    indices = Int[]
    @inbounds for i in eachindex(el.kind)
        el.kind[i] == kind && push!(indices, i)
    end
    return indices
end


"""
    event_times(el::EventLog, kind::EventKind)

Return event times for events whose kind is `kind`, preserving log order.
"""
function event_times(el::EventLog, kind::EventKind)
    times = Float64[]
    @inbounds for i in eachindex(el.kind)
        el.kind[i] == kind && push!(times, el.time[i])
    end
    return times
end


"""
    first_event_time(el::EventLog, kind::EventKind)

Return the first time for event kind `kind`, or `nothing` if absent.
"""
function first_event_time(el::EventLog, kind::EventKind)
    @inbounds for i in eachindex(el.kind)
        el.kind[i] == kind && return el.time[i]
    end
    return nothing
end


"""
    last_event_time(el::EventLog, kind::EventKind)

Return the last time for event kind `kind`, or `nothing` if absent.
"""
function last_event_time(el::EventLog, kind::EventKind)
    @inbounds for i in lastindex(el.kind):-1:firstindex(el.kind)
        el.kind[i] == kind && return el.time[i]
    end
    return nothing
end


"""
    has_event_kind(el::EventLog, kind::EventKind)

Return whether `el` contains at least one event of `kind`.
"""
has_event_kind(el::EventLog, kind::EventKind) = first_event_time(el, kind) !== nothing


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


Base.length(summary::HostEventSummary) = length(summary.host_id)


function Base.show(io::IO, summary::HostEventSummary)
    print(io, "HostEventSummary(", length(summary), " hosts)")
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
