
"""
    EventKind

Kind of event recorded in an [`EventLog`](@ref).

- `EK_Seeding`: host is present in the initial exposed/infected set at time zero.
- `EK_Transmission`: host was newly infected/exposed by `infector`.
- `EK_Activation`: host moved from exposed to infectious.
- `EK_Removal`: host left the infectious state without being sampled.
- `EK_FossilizedSampling`: host was sampled but remains infectious.
- `EK_SerialSampling`: host was sampled and removed from the infectious state.
- `EK_None`: internal sentinel; it is not expected in a completed event log.
"""
@enum EventKind::UInt8 begin
    EK_None = 0
    EK_Seeding = 1
    EK_Transmission = 2
    EK_Removal = 3
    EK_Activation = 4
    EK_FossilizedSampling = 5
    EK_SerialSampling = 6
end


const EVENTKIND_LABELS = Dict(
    EK_None              => "None",
    EK_Seeding           => "Seeding",
    EK_Transmission      => "Transmission",
    EK_Removal           => "Removal",
    EK_Activation        => "Activation",
    EK_FossilizedSampling => "FossilizedSampling",
    EK_SerialSampling    => "SerialSampling"
)

const EVENTKIND_CANONICAL_NAMES = (
    :seeding,
    :transmission,
    :activation,
    :removal,
    :fossilized_sampling,
    :fossilised_sampling,
    :serial_sampling,
)

const EVENTKIND_ALIASES = Dict(
    :seeding => EK_Seeding,
    :seedings => EK_Seeding,
    :transmission => EK_Transmission,
    :transmissions => EK_Transmission,
    :activation => EK_Activation,
    :activations => EK_Activation,
    :removal => EK_Removal,
    :removals => EK_Removal,
    :fossilized_sampling => EK_FossilizedSampling,
    :fossilized_samplings => EK_FossilizedSampling,
    :fossilised_sampling => EK_FossilizedSampling,
    :fossilised_samplings => EK_FossilizedSampling,
    :serial_sampling => EK_SerialSampling,
    :serial_samplings => EK_SerialSampling,
)


"""
    event_kind(kind)

Return the [`EventKind`](@ref) for a user-facing event-kind name.

Existing enum values are returned unchanged. Symbols and strings such as
`:transmission`, `:activation`, `:removal`, `:fossilized_sampling`,
`:fossilised_sampling`, and `:serial_sampling` are accepted for public
event-log helper APIs. Plural forms are accepted consistently for the same
canonical names.
"""
event_kind(kind::EventKind) = kind


event_kind(kind::Symbol) = _event_kind_from_name(String(kind), ":$kind")
event_kind(kind::AbstractString) = _event_kind_from_name(kind, repr(kind))


function _event_kind_from_name(kind::AbstractString, original)
    normalized = Symbol(lowercase(strip(kind)))
    return get(EVENTKIND_ALIASES, normalized) do
        throw(ArgumentError("unknown event kind $original; accepted names are " *
                            join((":" * String(name) for name in EVENTKIND_CANONICAL_NAMES), ", ") *
                            " (plural forms are also accepted)"))
    end
end


function Base.show(io::IO, x::EventKind)
    print(io, EVENTKIND_LABELS[x])
end


"""
    EventLog

Columnar event log produced by EpiSim simulation engines.

The four vectors must have equal length and are interpreted row-wise:

- `time[i]`: event time. Completed simulation logs are expected to be finite,
  non-negative, and monotonically non-decreasing.
- `host[i]`: one-based host identifier affected by the event.
- `infector[i]`: transmission source host for `EK_Transmission`; otherwise `0`.
- `kind[i]`: event kind.

Downstream consumers may assume validated logs use stable host identifiers within
one log, never use `0` as a host id, and use `infector == 0` as the no-source
sentinel. Logs do not encode complete model parameters, initial exposed versus
infectious status, unobserved non-events, or a tree-native representation.
"""
struct EventLog
    time::Vector{Float64}
    host::Vector{Int}
    infector::Vector{Int}
    kind::Vector{EventKind}
end

"""
    EventRecord

One row from an [`EventLog`](@ref), returned by indexing or iterating the log.
"""
struct EventRecord
    time::Float64
    host::Int
    infector::Int
    kind::EventKind
end


Base.length(el::EventLog) = length(el.time)
Base.isempty(el::EventLog) = isempty(el.time)
Base.firstindex(::EventLog) = 1
Base.lastindex(el::EventLog) = length(el)
Base.eltype(::Type{EventLog}) = EventRecord


function Base.getindex(el::EventLog, i::Integer)
    return EventRecord(el.time[i], el.host[i], el.infector[i], el.kind[i])
end


function Base.iterate(el::EventLog, state::Int=firstindex(el))
    state > lastindex(el) && return nothing
    return (el[state], state + 1)
end


function Base.show(io::IO, el::EventLog)
    print(io, "EventLog(", length(el), " events")
    if isempty(el)
        print(io, ", empty")
    else
        print(io, ", time range ", el.time[1], " to ", el.time[end])
    end
    print(io, ")")
end


function Base.show(io::IO, record::EventRecord)
    print(io, "EventRecord(time=", _compact_time(record.time),
          ", host=", record.host)
    if record.kind == EK_Transmission
        print(io, ", infector=", record.infector)
    end
    print(io, ", kind=")
    show(io, record.kind)
    print(io, ")")
end


_compact_time(t::Float64) = string(round(t; digits=3))


function EventLog(n::Int)
    return EventLog(zeros(Float64, n), collect(1:n), zeros(Int, n), fill(EK_Seeding, n))
end


function update_event_log!(el::EventLog,
                              t::Float64,
                              host::Int,
                              infector::Int,
                              kind::EventKind)
    push!(el.time, t)
    push!(el.host, host)
    push!(el.infector, infector)
    push!(el.kind, kind)
end


function _eventlog_validation_error(message::AbstractString, throw::Bool)
    throw && error(message)
    return false
end


"""
    validate_event_log(el::EventLog; population_size=nothing, throw=true)

Check the semantic contract for an EpiSim event log.

Validation covers equal column lengths, finite non-negative monotone times,
positive host ids, optional population bounds, event-specific source rules, and
duplicate infection/seeding of the same host. `EK_None` is rejected because it is
an internal sentinel rather than a completed event. When `throw=false`, the
function returns `false` instead of raising an error.

This is intentionally a log-level validator. It does not reconstruct hidden
state such as whether a seeded host began exposed or infectious.
"""
function validate_event_log(el::EventLog; population_size::Union{Nothing,Int}=nothing, throw::Bool=true)
    n = length(el.time)
    if length(el.host) != n || length(el.infector) != n || length(el.kind) != n
        return _eventlog_validation_error("EventLog columns must have equal lengths", throw)
    end

    if population_size !== nothing && population_size < 0
        return _eventlog_validation_error("population_size must be non-negative", throw)
    end

    infected_hosts = Set{Int}()
    previous_time = -Inf

    for i in 1:n
        t = el.time[i]
        host = el.host[i]
        source = el.infector[i]
        kind = el.kind[i]

        if !isfinite(t) || t < 0.0
            return _eventlog_validation_error("Event $i has invalid time $t", throw)
        end
        if t < previous_time
            return _eventlog_validation_error("Event $i occurs at $t before previous time $previous_time", throw)
        end
        previous_time = t

        if host < 1
            return _eventlog_validation_error("Event $i has invalid host id $host", throw)
        end
        if population_size !== nothing && host > population_size
            return _eventlog_validation_error("Event $i host id $host exceeds population_size $population_size", throw)
        end

        if kind == EK_None
            return _eventlog_validation_error("Event $i uses EK_None, which is an internal sentinel", throw)
        elseif kind == EK_Transmission
            if source < 1
                return _eventlog_validation_error("Transmission event $i must have a positive infector id", throw)
            end
            if source == host
                return _eventlog_validation_error("Transmission event $i has host equal to infector", throw)
            end
            if population_size !== nothing && source > population_size
                return _eventlog_validation_error("Transmission event $i source id $source exceeds population_size $population_size", throw)
            end
            if source ∉ infected_hosts
                return _eventlog_validation_error("Transmission event $i source id $source has not appeared as seeded or infected", throw)
            end
            if host ∈ infected_hosts
                return _eventlog_validation_error("Transmission event $i reinfects host $host", throw)
            end
            push!(infected_hosts, host)
        else
            if source != 0
                return _eventlog_validation_error("Non-transmission event $i must have infector/source id 0", throw)
            end
            if kind == EK_Seeding
                if !iszero(t)
                    return _eventlog_validation_error("Seeding event $i occurs at nonzero time $t", throw)
                end
                if host ∈ infected_hosts
                    return _eventlog_validation_error("Seeding event $i duplicates infected host $host", throw)
                end
                push!(infected_hosts, host)
            elseif host ∉ infected_hosts
                return _eventlog_validation_error("Event $i refers to host $host before seeding or transmission", throw)
            end
        end
    end

    return true
end


# Write iterators for time, host, infector, kind
eachtime(el::EventLog) = (el.time[i] for i in 1:length(el))
eachhost(el::EventLog) = (el.host[i] for i in 1:length(el))
eachinfector(el::EventLog) = (el.infector[i] for i in 1:length(el))
eachkind(el::EventLog) = (el.kind[i] for i in 1:length(el))
