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