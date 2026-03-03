abstract type AbstractStoppingCriterion end

struct NoStopping <: AbstractStoppingCriterion end

struct StopWhenTimeReached <: AbstractStoppingCriterion
    tmax::Float64
end

struct StopWhenCumulativeInfected <: AbstractStoppingCriterion
    threshold::Int
end

struct StopWhenCumulativeSampled <: AbstractStoppingCriterion
    threshold::Int
end

struct CompositeStopping <: AbstractStoppingCriterion
    criteria::Vector{AbstractStoppingCriterion}
end


should_stop(::NoStopping, state::SimulationState) = false

should_stop(c::StopWhenTimeReached, state) =
    state.t ≥ c.tmax

should_stop(c::StopWhenCumulativeInfected, state) =
    state.cumulative_infected ≥ c.threshold

should_stop(c::StopWhenCumulativeSampled, state) =
    state.cumulative_sampled ≥ c.threshold

should_stop(c::CompositeStopping, state) =
    any(crit -> should_stop(crit, state), c.criteria)