struct SIRParameters <: AbstractParameters
    N::Int
    transmission_rate::Float64
    recovery_rate::Float64
    sampling_rate::Float64
end

# Convenience keyword constructor
SIRParameters(; N::Int=10_000, transmission_rate::Float64=2.0, recovery_rate::Float64=0.9, sampling_rate::Float64=0.1) =
    SIRParameters(N, transmission_rate, recovery_rate, sampling_rate)


const SIR_EVENT_TYPES = [Transmission, Recovery, Sampling]


mutable struct AgenticSIRState <: AgenticState
    t::Float64
    S::Int
    I::Int
    currently_infected::Vector{Int}
    n_sampled::Int
    n_cumulative::Int
end


function AgenticSIRState(; t::Float64=0., 
                           S::Int=9999, 
                           I::Int=1)

    currently_infected = collect(1:I) 
    n_sampled = 0 
    n_cumulative = I 
    return AgenticSIRState(t, S, I, currently_infected, n_sampled, n_cumulative)
end


mutable struct AggregateSIRState <: AggregateState
    t::Float64
    S::Int
    I::Int
end


AggregateSIRState(; t::Float64=0., S::Int=9999, I::Int=1) = AggregateSIRState(t, S, I)


SIRState = Union{AgenticSIRState, AggregateSIRState}


function SIRModel(; N::Int=10_000, 
                    transmission_rate::Float64=2.0, 
                    recovery_rate::Float64=0.9, 
                    sampling_rate::Float64=0.1, 
                    I::Int=1, 
                    agentic::Bool=true)
    par = SIRParameters(; N, transmission_rate, recovery_rate, sampling_rate)
    state = agentic ?
        AgenticSIRState(; S=N-I, I=I) :
        AggregateSIRState(; S=N-I, I=I)
    return Model(par, SIR_EVENT_TYPES, state)
end


capture(state::SIRState) = (; t=state.t, S=state.S, I=state.I)


@inline function update_event_rates!(event_rates::Vector{Float64}, 
                                     par::SIRParameters, 
                                     state::SIRState)
    N = par.N
    β = par.transmission_rate
    α = par.recovery_rate
    ψ = par.sampling_rate

    S = state.S
    I = state.I

    event_rates[1] = β * S * I / N    # Transmission rate
    event_rates[2] = α * I            # Recovery rate
    event_rates[3] = ψ * I            # Sampling rate
end


# TODO: Perhaps still return an event log aggregate models
@inline function update_state!(rng::AbstractRNG,
                               par::SIRParameters,
                               state::AggregateSIRState, 
                               ::Type{Transmission})::Transmission
    state.S -= 1
    state.I += 1
    return Transmission(nothing, nothing, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               par::SIRParameters,
                               state::AgenticSIRState, 
                               ::Type{Transmission})::Transmission
    # Update state for Transmission event
    state.S -= 1
    state.I += 1
    state.n_cumulative += 1
    infectee = state.n_cumulative         # Label infected individuals sequentially
    infector = sample(rng, state.currently_infected)
    push!(state.currently_infected, infectee)
    return Transmission(infector, infectee, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               par::SIRParameters,
                               state::AggregateSIRState, 
                               ::Type{Recovery})::Recovery
    state.I -= 1
    return Recovery(nothing, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               par::SIRParameters,
                               state::AgenticSIRState, 
                               ::Type{Recovery})::Recovery
    # Update state for Recovery event
    state.I -= 1
    recovered = pop_random!(rng, state.currently_infected)
    return Recovery(recovered, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               par::SIRParameters,
                               state::AggregateSIRState, 
                               ::Type{Sampling})::Sampling
    state.I -= 1
    return Sampling(nothing, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               par::SIRParameters,
                               state::AgenticSIRState, 
                               ::Type{Sampling})::Sampling
    # Update state for Sampling event
    state.I -= 1
    sampled = pop_random!(rng, state.currently_infected)
    state.n_sampled += 1
    return Sampling(sampled, state.t)
end