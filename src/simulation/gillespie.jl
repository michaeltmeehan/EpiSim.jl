using ..Models


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
                   model::M;
                   state::S=EpiState(model),
                   stop_condition::Function=get_default_stop_condition(model)) where {M <: AbstractEpiModel, S <: AbstractEpiState}

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
        event = update_state!(rng, model, state, event_type)

        # Update event log
        push!(event_log, event)

    end
    return event_log
end