using ..Models


abstract type AbstractOutbreak end
abstract type AbstractEnsemble end

struct Outbreak{M<:AbstractEpiModel} <: AbstractOutbreak
    model::M
    state_log::DataFrame
    event_log::Union{Nothing, Vector{AbstractEpiEvent}}
end


struct Ensemble{M<:AbstractEpiModel} <: AbstractEnsemble
    model::M
    replicates::Vector{Outbreak{M}}
    seeds::Vector{Int}
end


function sample_event_type(rng::AbstractRNG,
                           event_types::Vector{DataType}, 
                           event_rates::Vector{Float64}, 
                           total_event_rate::Float64)
    r = rand(rng)
    cumulative_probability = event_rates[1] / total_event_rate
    event_index = 1
    while r > cumulative_probability && event_index < length(event_rates)
        event_index += 1
        cumulative_probability += event_rates[event_index] / total_event_rate
    end
    return event_types[event_index]
end


function simulate(rng::AbstractRNG,
                  model::AbstractEpiModel;
                  stop_condition::Function=get_default_stop_condition(model))

    state = deepcopy(model.initial_state)
    event_types = model.event_types
    event_log = isagentic(model) ? initialize_event_log(state) : nothing
    event_rates = Vector{Float64}(undef, length(event_types))

    state_log = [capture(state)]

    # Main simulation loop
    while !stop_condition(state)

        update_event_rates!(event_rates, model.parameters, state)

        total_event_rate = sum(event_rates)

        total_event_rate == 0.0 && break

        # Generate random number for time step and event selection
        # rand_number = rand(rng)

        # Update time based on the total event rate
        state.t -= log(rand(rng)) / total_event_rate

        # Sample an event type based on the event rates
        event_type = sample_event_type(rng, event_types, event_rates, total_event_rate)

        # Update model state and extract concrete event record (e.g., Transmission(1, 2, 1.0))
        event = update_state!(rng, model.parameters, state, event_type)

        # Update event log
        !isnothing(event) && push!(event_log, event)

        # Update state log
        push!(state_log, capture(state))
    end
    return Outbreak(model, DataFrame(state_log), event_log)
end


function simulate(model::AbstractEpiModel;
                  stop_condition::Function=get_default_stop_condition(model))
    return simulate(Random.GLOBAL_RNG, model; stop_condition=stop_condition)
end


function simulate(rng::AbstractRNG,
                  model::M,
                  n::Int;
                  stop_condition::Function = get_default_stop_condition(model)) where M <: AbstractEpiModel

    if n == 1
        return simulate(rng, model; stop_condition=stop_condition)
    else
        seeds = rand(rng, UInt32(1):UInt32(2^31-1), n)
        replicates = Vector{Outbreak{M}}(undef, n)
        Threads.@threads for i in 1:n
            replicates[i] = simulate(MersenneTwister(seeds[i]), model; stop_condition=stop_condition)
        end
        return Ensemble{M}(model, replicates, seeds)
    end
end


function simulate(model::M, n::Int;
                  stop_condition::Function = get_default_stop_condition(model)) where M <: AbstractEpiModel
    return simulate(Random.GLOBAL_RNG, model, n; stop_condition=stop_condition)
end
