struct BDParameters <: AbstractParameters
    birth_rate::Float64
    death_rate::Float64
    sampling_rate::Float64
end

# Convenience keyword constructor
BDParameters(; birth_rate=2.0, death_rate=0.9, sampling_rate=0.1) =
    BDParameters(birth_rate, death_rate, sampling_rate)


const BD_EVENT_TYPES = [Transmission, Recovery, Sampling]


mutable struct AgenticBDState <: AgenticState
    t::Float64
    I::Int
    currently_infected::Vector{Int}
    n_sampled::Int
    n_cumulative::Int
end

function AgenticBDState(; t=0., I=1)
    currently_infected = collect(1:I)
    n_sampled = 0
    n_cumulative = I
    return AgenticBDState(t, I, currently_infected, n_sampled, n_cumulative)
end


mutable struct AggregateBDState <: AggregateState
    t::Float64
    I::Int
end

AggregateBDState(; t=0., I=1) = AggregateBDState(t, I)


BDState = Union{AgenticBDState, AggregateBDState}


function BDModel(; birth_rate=2.0, death_rate=0.9, sampling_rate=0.1, I=1, agentic=true)
    par = BDParameters(; birth_rate, death_rate, sampling_rate)
    state = agentic ?
        AgenticBDState(; I=I) :
        AggregateBDState(; I=I)
    return Model(par, BD_EVENT_TYPES, state)
end


@inline function update_event_rates!(event_rates::Vector{Float64}, 
                                     par::BDParameters, 
                                     state::BDState)
    λ = par.birth_rate
    μ = par.death_rate
    ψ = par.sampling_rate

    I = state.I

    event_rates[1] = λ * I    # Birth rate
    event_rates[2] = μ * I    # Death rate
    event_rates[3] = ψ * I # Sampling rate
end


@inline function update_state!(rng::AbstractRNG,
                               parameters::BDParameters,
                               state::AggregateBDState, 
                               ::Type{Transmission})::Transmission
    state.I += 1
    return Transmission(nothing, nothing, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               parameters::BDParameters,
                               state::AgenticBDState, 
                               ::Type{Transmission})::Transmission
    # Update state for Transmission event
    state.I += 1
    state.n_cumulative += 1
    infectee = state.n_cumulative         # Label infected individuals sequentially
    infector = sample(rng, state.currently_infected)
    push!(state.currently_infected, infectee)
    return Transmission(infector, infectee, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               parameters::BDParameters,
                               state::AggregateBDState, 
                               ::Type{Recovery})::Recovery
    state.I -= 1
    return Recovery(nothing, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               parameters::BDParameters,
                               state::AgenticBDState, 
                               ::Type{Recovery})::Recovery
    # Update state for Recovery event
    state.I -= 1
    recovered = pop_random!(rng, state.currently_infected)
    return Recovery(recovered, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               parameters::BDParameters,
                               state::AggregateBDState, 
                               ::Type{Sampling})::Sampling
    state.I -= 1
    return Sampling(nothing, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               parameters::BDParameters,
                               state::AgenticBDState, 
                               ::Type{Sampling})::Sampling
    # Update state for Sampling event
    state.I -= 1
    sampled = pop_random!(rng, state.currently_infected)
    state.n_sampled += 1
    return Sampling(sampled, state.t)
end