"""
    StateCountTrajectory

Compact event-time compartment-count trajectory derived from one [`EventLog`](@ref).

The vectors are aligned by index. Entry `1` is the caller-provided initial state
at time `0.0`, before processing any event-log rows. Each subsequent entry is
the state immediately after processing one event. Simultaneous events are kept
as separate entries with repeated times, preserving raw event-log order.
"""
struct StateCountTrajectory
    time::Vector{Float64}
    S::Vector{Int}
    E::Vector{Int}
    I::Vector{Int}
    R::Vector{Int}
end

"""
    StateCountPoint

One event-time SEIR count point from a [`StateCountTrajectory`](@ref).
"""
struct StateCountPoint
    time::Float64
    S::Int
    E::Int
    I::Int
    R::Int
end


Base.length(traj::StateCountTrajectory) = length(traj.time)
Base.isempty(traj::StateCountTrajectory) = isempty(traj.time)
Base.firstindex(::StateCountTrajectory) = 1
Base.lastindex(traj::StateCountTrajectory) = length(traj)
Base.eltype(::Type{StateCountTrajectory}) = StateCountPoint


function Base.getindex(traj::StateCountTrajectory, i::Integer)
    return StateCountPoint(traj.time[i], traj.S[i], traj.E[i], traj.I[i], traj.R[i])
end


function Base.iterate(traj::StateCountTrajectory, state::Int=firstindex(traj))
    state > lastindex(traj) && return nothing
    return (traj[state], state + 1)
end


function Base.show(io::IO, traj::StateCountTrajectory)
    print(io, "StateCountTrajectory(", length(traj), " points")
    if isempty(traj)
        print(io, ", empty")
    else
        print(io, ", time range ", first(traj.time), " to ", last(traj.time))
    end
    print(io, ")")
end


function Base.show(io::IO, point::StateCountPoint)
    print(io, "StateCountPoint(time=", point.time,
          ", S=", point.S,
          ", E=", point.E,
          ", I=", point.I,
          ", R=", point.R,
          ")")
end


"""
    event_time_state_counts(el::EventLog; S0, E0, I0, R0=0, validate=true)

Recover a stepwise SEIR count trajectory at raw event times from one event log.

`S0`, `E0`, `I0`, and `R0` are the state counts at time zero before walking the
log. `EK_SerialSampling` is sampling with removal and is interpreted as an
`I -> R` transition. `EK_FossilizedSampling` is sampling without removal and
does not change compartment counts. `EK_Seeding` does not change counts during
recovery because seeded hosts must already be included in `E0` or `I0` according
to the simulation setup. The returned trajectory includes the initial state as
its first entry and then one entry after each event-log row.
"""
function event_time_state_counts(el::EventLog; S0::Integer, E0::Integer, I0::Integer, R0::Integer=0, validate::Bool=true)
    validate && validate_event_log(el)

    S = Int(S0)
    E = Int(E0)
    I = Int(I0)
    R = Int(R0)
    _check_nonnegative_state(S, E, I, R)

    n = length(el)
    time = Vector{Float64}(undef, n + 1)
    out_S = Vector{Int}(undef, n + 1)
    out_E = Vector{Int}(undef, n + 1)
    out_I = Vector{Int}(undef, n + 1)
    out_R = Vector{Int}(undef, n + 1)

    time[1] = 0.0
    out_S[1] = S
    out_E[1] = E
    out_I[1] = I
    out_R[1] = R

    @inbounds for i in 1:n
        kind = el.kind[i]

        if kind == EK_Transmission
            S -= 1
            E += 1
        elseif kind == EK_Activation
            E -= 1
            I += 1
        elseif kind == EK_Removal || kind == EK_SerialSampling
            I -= 1
            R += 1
        elseif kind == EK_Seeding || kind == EK_FossilizedSampling
            nothing
        else
            throw(ArgumentError("unsupported event kind $kind at row $i"))
        end

        _check_nonnegative_state(S, E, I, R)
        j = i + 1
        time[j] = el.time[i]
        out_S[j] = S
        out_E[j] = E
        out_I[j] = I
        out_R[j] = R
    end

    return StateCountTrajectory(time, out_S, out_E, out_I, out_R)
end


function _check_nonnegative_state(S::Int, E::Int, I::Int, R::Int)
    if S < 0 || E < 0 || I < 0 || R < 0
        throw(ArgumentError("state counts must remain non-negative"))
    end
    return nothing
end
