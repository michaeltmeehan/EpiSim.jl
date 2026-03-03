# ================================
# Canonical EventLog definition
# ================================

export EventLog, push_event!

"""
Struct-of-arrays event log.

All simulation engines must emit events into this structure.
"""
mutable struct EventLog
    t::Vector{Time}
    kind::Vector{EventKind}
    host::Vector{HostID}
    src::Vector{HostID}
end

EventLog() = EventLog(Time[], EventKind[], HostID[], HostID[])

function push_event!(
    log::EventLog,
    t::Time,
    kind::EventKind,
    host::HostID,
    src::HostID = 0,
)
    push!(log.t, t)
    push!(log.kind, kind)
    push!(log.host, host)
    push!(log.src, src)
    return nothing
end


"""
Lightweight row view of an EventLog entry.
"""
struct Event
    t::Time
    kind::EventKind
    host::HostID
    src::HostID
end


Base.length(log::EventLog) = length(log.t)

function Base.iterate(log::EventLog, i=1)
    i > length(log.t) && return nothing
    return (
        Event(log.t[i], log.kind[i], log.host[i], log.src[i]),
        i + 1
    )
end


Base.getindex(log::EventLog, i::Int) =
    Event(log.t[i], log.kind[i], log.host[i], log.src[i])


Base.size(log::EventLog) = (length(log),)


Base.iterate(r::Base.Iterators.Reverse{EventLog}, i=length(r.itr)) =
    i < 1 ? nothing :
    (Event(r.itr.t[i], r.itr.kind[i], r.itr.host[i], r.itr.src[i]), i - 1)

Base.reverse(log::EventLog) = Base.Iterators.Reverse(log)