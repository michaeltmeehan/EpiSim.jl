struct SIRParameters <: AbstractEpiParameters
    N::Int
    transmission_rate::Float64
    recovery_rate::Float64
    sampling_rate::Float64
end

# Convenience keyword constructor
SIRParameters(; N::Int=10_000, transmission_rate::Float64=2.0, recovery_rate::Float64=0.9, sampling_rate::Float64=0.1) =
    SIRParameters(N, transmission_rate, recovery_rate, sampling_rate)

calc_R0(par::SIRParameters) = par.transmission_rate / (par.recovery_rate + par.sampling_rate)

calc_infectious_period(par::SIRParameters) = 1.0 / (par.recovery_rate + par.sampling_rate)

calc_sampling_fraction(par::SIRParameters) = par.sampling_rate / (par.recovery_rate + par.sampling_rate)

calc_extinction_probability(par::SIRParameters) = 1. / calc_R0(par)


# TODO: Generalize this function to p::AbstractEpiParameters
function summarize(p::SIRParameters)
    R₀ = calc_R0(p)
    infectious_period = calc_infectious_period(p)
    return (;R₀ = R₀, infectious_period = infectious_period)
end


# TODO: This can probably be abstracted to Model{P<:AbstracctEpiParameters, S<:AbstractEpiState}
struct SIRModel{S<:AbstractEpiState} <: AbstractEpiModel
    parameters::SIRParameters
    event_types::Vector{DataType}
    initial_state::S
end


@forward SIRModel.parameters summarize

# TODO: These can probably be generalized to calc_R0(model::AbstractEpiModel)
calc_R0(model::SIRModel) = calc_R0(model.parameters)
calc_infectious_period(model::SIRModel) = calc_infectious_period(model.parameters)
calc_extinction_probability(model::SIRModel) = calc_extinction_probability(model.parameters)^model.initial_state.I


const SIR_EVENT_TYPES = [Transmission, Recovery, Sampling]


# TODO: Implement AbstractSIRState and have AgenticSIRState and AggregateSIRState inherit from it
mutable struct AgenticSIRState <: AgenticState
    t::Float64
    S::Int
    I::Int
    currently_infected::Vector{Int}
    n_sampled::Int
    n_cumulative::Int
end

AgenticSIRState(; t::Float64=0., S::Int=9999, I::Int=1, currently_infected::Vector{Int}=collect(1:I), n_sampled::Int=0, n_cumulative::Int=0) =
    AgenticSIRState(t, S, I, currently_infected, n_sampled, n_cumulative)


mutable struct AggregateSIRState <: AggregateState
    t::Float64
    S::Int
    I::Int
end

AggregateSIRState(; t::Float64=0., S::Int=9999, I::Int=1) = AggregateSIRState(t, S, I)


function SIRModel(; N::Int=10_000, transmission_rate::Float64=2.0, recovery_rate::Float64=0.9, sampling_rate::Float64=0.1, I::Int=1, agentic::Bool=true)
    par = SIRParameters(; N, transmission_rate, recovery_rate, sampling_rate)
    state = agentic ?
        AgenticSIRState(; S=N-I, I=I) :
        AggregateSIRState(; S=N-I, I=I)
    return SIRModel(par, SIR_EVENT_TYPES, state)
end



capture(state::Union{AgenticSIRState, AggregateSIRState}) = (; t=state.t, S=state.S, I=state.I)

get_default_stop_condition(model::SIRModel{AgenticSIRState}) = s -> s.n_sampled >= 100

get_default_stop_condition(model::SIRModel{AggregateSIRState}) = s -> s.I == 0 || s.t >= 100.0


function initialize_event_log(state::AgenticSIRState)::Vector{AbstractEpiEvent}
    event_log = Vector{AbstractEpiEvent}()
    for i in 1:state.I
        push!(event_log, Seed(i, 0.0))
    end
    return event_log
end

# TODO: Change parameters to par throughout
@inline function update_event_rates!(event_rates::Vector{Float64}, 
                                     par::SIRParameters, 
                                     state::Union{AgenticSIRState, AggregateSIRState})
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
                               ::Type{Transmission})::Nothing
    state.S -= 1
    state.I += 1
    return
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
                               ::Type{Recovery})::Nothing
    state.I -= 1
    return
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
                               ::Type{Sampling})::Nothing
    state.I -= 1
    return nothing
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