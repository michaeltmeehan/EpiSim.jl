const SIR_EVENT_TYPES = [Transmission, Recovery, Sampling]

const SIREvent = Union{Seed, SIR_EVENT_TYPES...}

get_event_types(model::SIRModel) = SIR_EVENT_TYPES


mutable struct SIRState <: AbstractEpiState
    t::Float64
    S::Int
    I::Int
    currently_infected::Vector{Int}
    n_sampled::Int
    n_cumulative::Int
end


EpiState(model::SIRModel) = SIRState(0.0, model.N - 1, 1, [1], 0, 1)


default_stop_condition(model::SIRModel) = s -> isempty(s.currently_infected)


function initialize_event_log(state::SIRState)::Vector{SIREvent}
    event_log = Vector{SIREvent}()
    for i in 1:state.I
        push!(event_log, Seed(i, 0.0))
    end
    return event_log
end


@inline function update_event_rates!(event_rates::Vector{Float64}, 
                                     model::SIRModel, 
                                     state::SIRState)
    β = model.transmission_rate
    γ = model.recovery_rate
    ψ = model.sampling_rate
    N = model.N

    S = state.S
    I = state.I

    event_rates[1] = β * I * S / N  # Infection rate
    event_rates[2] = γ * I        # Recovery rate
    event_rates[3] = ψ * I        # Sampling rate
end


function update_state!(rng::AbstractRNG,
                       state::SIRState, 
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


function update_state!(rng::AbstractRNG,
                       state::SIRState, 
                       ::Type{Recovery})::Recovery
    # Update state for Recovery event
    state.I -= 1
    recovered = pop_random!(rng, state.currently_infected)
    return Recovery(recovered, state.t)
end


function update_state!(rng::AbstractRNG,
                       state::SIRState, 
                       ::Type{Sampling})::Sampling
    # Update state for Sampling event
    state.I -= 1
    sampled = pop_random!(rng, state.currently_infected)
    state.n_sampled += 1
    return Sampling(sampled, state.t)
end


function update_event_rates!(event_rates::Vector{Float64},
                             model::AbstractEpiModel,
                             state::AbstractEpiState)
    error("update_event_rates! not implemented for $(typeof(model))")
end


function get_default_stop_condition(model::AbstractEpiModel)
    return s -> false
end



function sample_event_type(rand_number::Float64, 
                           event_types::Vector{DataType}, 
                           event_rates::Vector{Float64}, 
                           total_event_rate::Float64)
    cumulative_probability = event_rates[1] / total_event_rate
    event_index = 1
    while rand_number > cumulative_probability && event_index < length(event_rates)
        event_index += 1
        cumulative_probability += event_rates[event_index] / total_event_rate
    end
    return event_types[event_index]
end


function gillespie(rng::AbstractRNG,
                   model::AbstractEpiModel;
                   state::AbstractEpiState=EpiState(model),
                   stop_condition::Function=default_stop_condition(model))

    event_types = get_event_types(model)
    event_log = initialize_event_log(state)
    event_rates = Vector{Float64}(undef, length(event_types))

    # Main simulation loop
    while !stop_condition(state)

        update_event_rates!(event_rates, model, state)

        total_event_rate = sum(event_rates)

        # Generate random number for time step and event selection
        rand_number = rand(rng)

        # Update time based on the total event rate
        state.t -= log(rand_number) / total_event_rate

        # Sample an event type based on the event rates
        event_type = sample_event_type(rand_number, event_types, event_rates, total_event_rate)

        # Update model state and extract concrete event record (e.g., Transmission(1, 2, 1.0))
        event = update_state!(rng, state, event_type)

        # Update event log
        push!(event_log, event)

    end
    return event_log
end