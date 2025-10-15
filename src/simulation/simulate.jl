using ..Models

abstract type AbstractSimulation end
abstract type AbstractEnsemble end


struct Simulation{M<:AbstractModel} <: AbstractSimulation
    model::M
    state_log::DataFrame
    event_log::Vector{AbstractEvent}
end


eachstate(sim::Simulation) = eachrow(sim.state_log)
eachevent(sim::Simulation) = sim.event_log


function Base.show(io::IO, sim::Simulation)
    duration = sim.state_log.t[end] - sim.state_log.t[1]
    n_steps = nrow(sim.state_log)
    n_events = length(sim.event_log)
    counts = event_counts(sim)

    println(io, "Simulation of ", typeof(sim.model))
    println(io, "  Duration   :  ", round(duration, digits=3))
    println(io, "  Time steps :  ", n_steps)
    println(io, "  Events     :  ", n_events)
    for (etype, count) in sort(collect(counts); by=last, rev=true)
        println(io, "    ", rpad(nameof(etype), 14), ": ", count)
    end
end


# TODO: Consider formalizing the above as custom iterators
# struct EventIterator
#     sim::Simulation
# end

# Base.iterate(it::EventIterator) = iterate(it.sim.event_log)
# Base.iterate(it::EventIterator, state) = iterate(it.sim.event_log, state)

# eachevent(sim::Simulation) = EventIterator(sim)


struct Ensemble{M<:AbstractModel} <: AbstractEnsemble
    model::M
    simulations::Vector{Simulation{M}}
    seeds::Vector{Int}
end


eachsim(ens::Ensemble) = ens.simulations

Base.length(ens::Ensemble) = length(ens.simulations)
Base.getindex(ens::Ensemble, i::Int) = ens.simulations[i]
Base.getindex(ens::Ensemble, i::AbstractVector{<:Int}) = ens.simulations[i]

Base.iterate(ens::Ensemble) = iterate(ens.simulations)
Base.iterate(ens::Ensemble, state) = iterate(ens.simulations, state)
Base.eltype(ens::Ensemble{M}) where M = Simulation{M}

function ensemble_summary(ens::Ensemble)
    durations = [sim.state_log.t[end] for sim in eachsim(ens)]
    n_events  = [length(sim.event_log) for sim in eachsim(ens)]

    return (
        n = length(ens),
        duration_mean = mean(durations),
        duration_range = extrema(durations),
        event_mean = mean(n_events),
        event_range = extrema(n_events)
    )
end


function Base.show(io::IO, ens::Ensemble)
    s = ensemble_summary(ens)

    println(io, "Ensemble of ", s.n, " simulations (", typeof(ens.model), ")")
    println(io, "  Duration:   mean = $(round(s.duration_mean, digits=3))")
    println(io, "              range = ", round.(s.duration_range, digits=3))
    println(io, "  Events:     mean = $(round(s.event_mean, digits=1))")
    println(io, "              range = ", s.event_range)
end


function sample_event_type(rng::AbstractRNG,
                           event_types::Vector{DataType}, 
                           event_rates::Vector{Float64}, 
                           total_event_rate::Float64)
    r = rand(rng) * total_event_rate
    cumulative_weight = event_rates[1]
    event_index = 1
    while r > cumulative_weight && event_index < length(event_rates)
        event_index += 1
        cumulative_weight += event_rates[event_index]
    end
    return event_types[event_index]
end


function simulate(rng::AbstractRNG,
                  model::AbstractModel;
                  stop_condition::Function=get_default_stop_condition(model))

    state = deepcopy(model.initial_state)
    event_types = model.event_types
    event_rates = Vector{Float64}(undef, length(event_types))
    event_log = Vector{AbstractEvent}()
    push!(event_log, initialize_event_log(state)...)
    state_log = [capture(state)]

    # Main simulation loop
    while !stop_condition(state)

        update_event_rates!(event_rates, model.par, state)

        total_event_rate = sum(event_rates)

        total_event_rate == 0.0 && break

        # Update time based on the total event rate
        state.t -= log(rand(rng)) / total_event_rate

        # Sample an event type based on the event rates
        event_type = sample_event_type(rng, event_types, event_rates, total_event_rate)

        # Update model state and extract concrete event record (e.g., Transmission(1, 2, 1.0))
        event = update_state!(rng, model.par, state, event_type)

        # Update event log
        push!(event_log, event)

        # Update state log
        push!(state_log, capture(state))
    end

    # Check if final state violates stop condition
    if stop_condition(state)
        pop!(state_log)  # Remove last state if it violates stop condition
        pop!(event_log)  # Remove last event if it violates stop condition
    end
    return Simulation(model, DataFrame(state_log), event_log)
end


function simulate(model::AbstractModel;
                  stop_condition::Function=get_default_stop_condition(model))
    return simulate(Random.GLOBAL_RNG, model; stop_condition=stop_condition)
end


function simulate(rng::AbstractRNG,
                  model::M,
                  n::Int;
                  stop_condition::Function = get_default_stop_condition(model)) where M <: AbstractModel

    if n == 1
        return simulate(rng, model; stop_condition=stop_condition)
    else
        seeds = rand(rng, UInt32(1):UInt32(2^31-1), n)
        simulations = Vector{Simulation{M}}(undef, n)
        Threads.@threads for i in 1:n
            simulations[i] = simulate(MersenneTwister(seeds[i]), model; stop_condition=stop_condition)
        end
        return Ensemble{M}(model, simulations, seeds)
    end
end


function simulate(model::AbstractModel, 
                  n::Int;
                  stop_condition::Function = get_default_stop_condition(model))
    return simulate(Random.GLOBAL_RNG, model, n; stop_condition=stop_condition)
end
