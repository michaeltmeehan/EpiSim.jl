# ================================
# Core type definitions
# ================================

export HostID, Time, EventKind

const HostID = Int
const Time   = Float64

@enum EventKind::UInt8 begin
    Seeding
    Transmission
    Activation
    Removal
    SerialSampling
    FossilisedSampling
end


mutable struct SimulationState
    t::Float64
    cumulative_infected::Int
    cumulative_sampled::Int
end


function update_state!(state::SimulationState, t::Time, ev::EventKind)
    state.t = t
    if ev == Transmission
        state.cumulative_infected += 1
    elseif ev == SerialSampling || ev == FossilisedSampling
        state.cumulative_sampled += 1
    end
end