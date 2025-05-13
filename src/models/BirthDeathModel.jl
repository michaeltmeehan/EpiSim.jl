struct BirthDeathParameters <: AbstractEpiParameters
    birth_rate::Float64
    death_rate::Float64
    sampling_rate::Float64
end


struct BirthDeathModel{S<:AbstractEpiState} <: AbstractEpiModel
    parameters::BirthDeathParameters
    event_types::Vector{DataType}
    initial_state::S
end


get_default_stop_condition(model::BirthDeathModel{AgenticBirthDeathState}) = s -> s.n_sampled >= 100

get_default_stop_condition(model::BirthDeathModel{AggregateBirthDeathState}) = s -> s.I == 0 || s.t >= 100.0 || s.I >= 10_000

const BIRTHDEATH_EVENT_TYPES = [Transmission, Recovery, Sampling]

const BirthDeathEvent = Union{Seed, BIRTHDEATH_EVENT_TYPES...}


mutable struct AgenticBirthDeathState <: AgenticState
    t::Float64
    I::Int
    currently_infected::Vector{Int}
    n_sampled::Int
    n_cumulative::Int
end


mutable struct AggregateBirthDeathState <: AggregateState
    t::Float64
    I::Int
end


capture(state::Union{AgenticBirthDeathState, AggregateBirthDeathState}) = (; t=state.t, I=state.I)



function initialize_event_log(state::AgenticBirthDeathState)::Vector{BirthDeathEvent}
    event_log = Vector{BirthDeathEvent}()
    for i in 1:state.I
        push!(event_log, Seed(i, 0.0))
    end
    return event_log
end


@inline function update_event_rates!(event_rates::Vector{Float64}, 
                                     parameters::BirthDeathParameters, 
                                     state::Union{AgenticBirthDeathState, AggregateBirthDeathState})
    birth_rate = parameters.birth_rate
    death_rate = parameters.death_rate
    sampling_rate = parameters.sampling_rate

    I = state.I

    event_rates[1] = birth_rate * I    # Birth rate
    event_rates[2] = death_rate * I    # Death rate
    event_rates[3] = sampling_rate * I # Sampling rate
end


@inline function update_state!(rng::AbstractRNG,
                               parameters::BirthDeathParameters,
                               state::AggregateBirthDeathState, 
                               ::Type{Transmission})::Nothing
    state.I += 1
    return
end


@inline function update_state!(rng::AbstractRNG,
                               parameters::BirthDeathParameters,
                               state::AgenticBirthDeathState, 
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
                               parameters::BirthDeathParameters,
                               state::AggregateBirthDeathState, 
                               ::Type{Recovery})::Nothing
    state.I -= 1
    return
end


@inline function update_state!(rng::AbstractRNG,
                               parameters::BirthDeathParameters,
                               state::AgenticBirthDeathState, 
                               ::Type{Recovery})::Recovery
    # Update state for Recovery event
    state.I -= 1
    recovered = pop_random!(rng, state.currently_infected)
    return Recovery(recovered, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               parameters::BirthDeathParameters,
                               state::AggregateBirthDeathState, 
                               ::Type{Sampling})::Nothing
    state.I -= 1
    return nothing
end


@inline function update_state!(rng::AbstractRNG,
                               parameters::BirthDeathParameters,
                               state::AgenticBirthDeathState, 
                               ::Type{Sampling})::Sampling
    # Update state for Sampling event
    state.I -= 1
    sampled = pop_random!(rng, state.currently_infected)
    state.n_sampled += 1
    return Sampling(sampled, state.t)
end



function update_event_rates!(event_rates::Vector{Float64}, parameters::BirthDeathParameters)
    λ = parameters.birth_rate
    μ = parameters.death_rate
    ψ = parameters.sampling_rate

    event_rates[1] = λ  # Birth rate
    event_rates[2] = μ  # Death rate
    event_rates[3] = ψ  # Sampling rate
end


function simulate_events(rng::AbstractRNG,
                         model::BirthDeathModel;
                         N_max::Int=10_000, 
                         t_max::Float64=100.0, 
                         S_max::Int=100, 
                         I_init::Int=1)

    I = I_init
    n_cumulative = I_init
    n_sampled = 0

    currently_infected = collect(1:I_init)

    events = Vector{BirthDeathEvent}()

    for seed_host in 1:I_init
        push!(events, Seed(seed_host, 0.0))
    end

    # Pre-calculate event rates
    event_rates = Vector{Float64}(undef, 3)
    update_event_rates!(event_rates, model.parameters)
    total_event_rate = sum(event_rates)
    
    t = 0.0

    while !isempty(currently_infected) && n_cumulative < N_max && n_sampled < S_max
        rand_number = rand(rng)
        t -= log(rand_number) / (total_event_rate * I)
        t > t_max && break

        if rand_number ≤ event_rates[1] / total_event_rate
            # Birth event
            I += 1
            n_cumulative += 1
            infectee = n_cumulative         # Label infected individuals sequentially
            infector = sample(rng, currently_infected)
            push!(currently_infected, infectee)
            transmission!(events, infector, infectee, t)
        elseif rand_number ≤ (event_rates[1] + event_rates[2]) / total_event_rate
            # Death event
            I -= 1
            recovered = pop_random!(rng, currently_infected)
            recovery!(events, recovered, t)
        else
            # Sampling event
            I -= 1
            n_sampled += 1
            sampled = pop_random!(rng, currently_infected)
            sampling!(events, sampled, t)
        end
    end
    return events
end


function simulate_events(model::BirthDeathModel;
                           N_max::Int=10_000, 
                           t_max::Float64=100.0, 
                           S_max::Int=100, 
                           I_init::Int=1)
    return simulate_events(Random.GLOBAL_RNG, model; N_max=N_max, t_max=t_max, S_max=S_max, I_init=I_init)
end