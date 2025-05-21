struct SEIRParameters <: AbstractParameters
    N::Int
    transmission_rate::Float64
    activation_rate::Float64
    recovery_rate::Float64
    sampling_rate::Float64
end

# Convenience keyword constructor
SEIRParameters(; N::Int=10_000, transmission_rate::Float64=2.0, activation_rate::Float64=1.0, recovery_rate::Float64=0.9, sampling_rate::Float64=0.1) =
    SEIRParameters(N, transmission_rate, activation_rate, recovery_rate, sampling_rate)


const SEIR_EVENT_TYPES = [Transmission, Activation, Recovery, Sampling]


mutable struct AgenticSEIRState <: AgenticState
    t::Float64
    S::Int
    E::Int
    I::Int
    currently_exposed::Vector{Int}
    currently_infected::Vector{Int}
    n_sampled::Int
    n_cumulative::Int
end


function AgenticSEIRState(; t::Float64=0., 
                            S::Int=9999, 
                            E::Int=0, 
                            I::Int=1) 
    currently_exposed = collect(1:E) 
    currently_infected = collect(E+1:E+I)
    n_sampled = 0
    n_cumulative = E + I
    return AgenticSEIRState(t, S, E, I, currently_exposed, currently_infected, n_sampled, n_cumulative)
end


mutable struct AggregateSEIRState <: AggregateState
    t::Float64
    S::Int
    E::Int
    I::Int
end


AggregateSEIRState(; t::Float64=0., S::Int=9999, E::Int=0, I::Int=1) = AggregateSEIRState(t, S, E, I)


const SEIRState = Union{AgenticSEIRState, AggregateSEIRState}


function SEIRModel(; N::Int=10_000, 
                     transmission_rate::Float64=2.0,
                     activation_rate::Float64=1.0,
                     recovery_rate::Float64=0.9, 
                     sampling_rate::Float64=0.1,
                     E::Int=0,
                     I::Int=1, 
                     agentic::Bool=true)
    par = SEIRParameters(; N, transmission_rate, activation_rate, recovery_rate, sampling_rate)
    state = agentic ?
        AgenticSEIRState(; S=N-E-I, E=E, I=I) :
        AggregateSEIRState(; S=N-E-I, E=E, I=I)
    return Model(par, SEIR_EVENT_TYPES, state)
end


capture(state::SEIRState) = (; t=state.t, S=state.S, E=state.E, I=state.I)


@inline function update_event_rates!(event_rates::Vector{Float64}, 
                                     par::SEIRParameters, 
                                     state::SEIRState)
    N = par.N
    β = par.transmission_rate
    ν = par.activation_rate
    α = par.recovery_rate
    ψ = par.sampling_rate

    S = state.S
    E = state.E
    I = state.I

    event_rates[1] = β * S * I / N    # Transmission rate
    event_rates[2] = ν * E            # Activation rate
    event_rates[3] = α * I            # Recovery rate
    event_rates[4] = ψ * I            # Sampling rate
end


@inline function update_state!(rng::AbstractRNG,
                               par::SEIRParameters,
                               state::AggregateSEIRState, 
                               ::Type{Transmission})::Transmission
    state.S -= 1
    state.E += 1
    return Transmission(nothing, nothing, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               par::SEIRParameters,
                               state::AgenticSEIRState, 
                               ::Type{Transmission})::Transmission
    # Update state for Transmission event
    state.S -= 1
    state.E += 1
    state.n_cumulative += 1
    infectee = state.n_cumulative         # Label infected individuals sequentially
    infector = sample(rng, state.currently_infected)
    push!(state.currently_exposed, infectee)
    return Transmission(infector, infectee, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               par::SEIRParameters,
                               state::AggregateSEIRState, 
                               ::Type{Recovery})::Recovery
    state.I -= 1
    return Recovery(nothing, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               par::SEIRParameters,
                               state::AgenticSEIRState, 
                               ::Type{Recovery})::Recovery
    # Update state for Recovery event
    state.I -= 1
    recovered = pop_random!(rng, state.currently_infected)
    return Recovery(recovered, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               par::SEIRParameters,
                               state::AggregateSEIRState, 
                               ::Type{Sampling})::Sampling
    state.I -= 1
    return Sampling(nothing, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               par::SEIRParameters,
                               state::AgenticSEIRState, 
                               ::Type{Sampling})::Sampling
    # Update state for Sampling event
    state.I -= 1
    sampled = pop_random!(rng, state.currently_infected)
    state.n_sampled += 1
    return Sampling(sampled, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               par::SEIRParameters,
                               state::AggregateSEIRState, 
                               ::Type{Activation})::Activation
    # Update state for Activation event
    state.E -= 1
    state.I += 1
    return Activation(nothing, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               par::SEIRParameters,
                               state::AgenticSEIRState, 
                               ::Type{Activation})::Activation
    # Update state for Activation event
    state.E -= 1
    state.I += 1
    activated = pop_random!(rng, state.currently_exposed)
    push!(state.currently_infected, activated)
    return Activation(activated, state.t)
end