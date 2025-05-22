using ..Models


abstract type AbstractOutbreak end
abstract type AbstractEnsemble end

struct Outbreak{M<:AbstractModel} <: AbstractOutbreak
    model::M
    state_log::DataFrame
    event_log::Vector{<:AbstractEpiEvent}
end


struct Ensemble{M<:AbstractModel} <: AbstractEnsemble
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
                  model::AbstractModel;
                  stop_condition::Function=get_default_stop_condition(model),
                  max_iter::Int=100_000)

    state = deepcopy(model.initial_state)
    event_types = model.event_types
    event_rates = Vector{Float64}(undef, length(event_types))

    initial_state_log = capture(state)
    state_log = Vector{typeof(initial_state_log)}(undef, max_iter)
    state_log[1] = initial_state_log

    initial_event_log = initialize_event_log(state)
    event_log = Vector{Union{event_types..., [typeof(event) for event in initial_event_log]...}}(undef, max_iter)
    for (i, event) in enumerate(initial_event_log)
        event_log[i] = event
    end
    
    n = 1
    m = length(initial_event_log)
    # event_log = initialize_event_log(state)
    # state_log = [capture(state)]

    # Main simulation loop
    while !stop_condition(state) && n < max_iter && m < max_iter

        update_event_rates!(event_rates, model.par, state)

        total_event_rate = sum(event_rates)

        total_event_rate == 0.0 && break

        # Update time based on the total event rate
        state.t -= log(rand(rng)) / total_event_rate

        # Sample an event type based on the event rates
        event_type = sample_event_type(rng, event_types, event_rates, total_event_rate)

        # Update model state and extract concrete event record (e.g., Transmission(1, 2, 1.0))
        event = update_state!(rng, model.par, state, event_type)

        # Update iteration counters
        n += 1
        m += 1

        # Update event log
        event_log[m] = event
        # push!(event_log, event)

        # Update state log
        state_log[n] = capture(state)
        # push!(state_log, capture(state))
    end
    return Outbreak(model, DataFrame(state_log[1:n]), event_log[1:m])
    # return Outbreak(model, DataFrame(state_log), event_log)
end


function simulate(model::AbstractModel;
                  stop_condition::Function=get_default_stop_condition(model),
                  max_iter::Int=100_000)
    return simulate(Random.GLOBAL_RNG, model; stop_condition=stop_condition)
end


function simulate(rng::AbstractRNG,
                  model::M,
                  n::Int;
                  stop_condition::Function = get_default_stop_condition(model),
                  max_iter::Int=100_000) where M <: AbstractModel

    if n == 1
        return simulate(rng, model; stop_condition=stop_condition, max_iter=max_iter)
    else
        seeds = rand(rng, UInt32(1):UInt32(2^31-1), n)
        replicates = Vector{Outbreak{M}}(undef, n)
        Threads.@threads for i in 1:n
            replicates[i] = simulate(MersenneTwister(seeds[i]), model; stop_condition=stop_condition, max_iter=max_iter)
        end
        return Ensemble{M}(model, replicates, seeds)
    end
end


function simulate(model::AbstractModel, 
                  n::Int;
                  stop_condition::Function = get_default_stop_condition(model),
                  max_iter::Int=100_000)
    return simulate(Random.GLOBAL_RNG, model, n; stop_condition=stop_condition, max_iter=max_iter)
end
