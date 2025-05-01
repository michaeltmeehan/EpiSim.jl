struct BirthDeathModel <: AbstractEpiModel
    birth_rate::Float64
    death_rate::Float64
    sampling_rate::Float64
end

const BirthDeath_EVENT_TYPES = [Transmission, Recovery, Sampling]

const BirthDeathEvent = Union{Seed, BirthDeath_EVENT_TYPES...}

get_event_types(model::BirthDeathModel) = BirthDeath_EVENT_TYPES


mutable struct BirthDeathState <: AbstractEpiState
    t::Float64
    I::Int
    currently_infected::Vector{Int}
    n_sampled::Int
    n_cumulative::Int
end


function validate_state(state::BirthDeathState)
    VALIDATE_STATE || return
    @assert state.I >= 0 "Infected individuals cannot be negative."
    @assert state.n_sampled >= 0 "Sampled individuals cannot be negative."
    @assert state.n_cumulative >= 0 "Cumulative infected individuals cannot be negative."
    @assert state.t >= 0 "Time cannot be negative."
    # for i in state.currently_infected
    #     @assert i ≥ 0 "Currently infected individuals cannot be negative."
    #     @assert i ≤ state.n_cumulative "Currently infected individuals cannot exceed cumulative count."
    # end
    # @assert all(state.currently_infected .≥ 0) "Currently infected individuals cannot be negative."
    # @assert all(state.currently_infected .<= state.n_cumulative) "Currently infected individuals cannot exceed cumulative infected individuals."
    @assert state.I == length(state.currently_infected) "Number of currently infected individuals does not match the length of the currently infected vector."
end


EpiState(model::BirthDeathModel) = BirthDeathState(0.0, 1, [1], 0, 1)


get_default_stop_condition(model::BirthDeathModel) = s -> isempty(s.currently_infected) || s.n_cumulative >= 10_000 || s.n_sampled >= 100 || s.t >= 100.0


function initialize_event_log(state::BirthDeathState)::Vector{BirthDeathEvent}
    event_log = Vector{BirthDeathEvent}()
    for i in 1:state.I
        push!(event_log, Seed(i, 0.0))
    end
    return event_log
end


@inline function update_event_rates!(event_rates::Vector{Float64}, 
                                     model::BirthDeathModel, 
                                     state::BirthDeathState)
    birth_rate = model.birth_rate
    death_rate = model.death_rate
    sampling_rate = model.sampling_rate

    I = state.I

    event_rates[1] = birth_rate * I    # Birth rate
    event_rates[2] = death_rate * I    # Death rate
    event_rates[3] = sampling_rate * I # Sampling rate
end


@inline function update_state!(rng::AbstractRNG,
                               model::BirthDeathModel,
                               state::BirthDeathState, 
                               ::Type{Transmission})::Transmission
    # Update state for Transmission event
    state.I += 1
    state.n_cumulative += 1
    infectee = state.n_cumulative         # Label infected individuals sequentially
    infector = sample(rng, state.currently_infected)
    push!(state.currently_infected, infectee)
    validate_state(state)
    return Transmission(infector, infectee, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               model::BirthDeathModel,
                               state::BirthDeathState, 
                               ::Type{Recovery})::Recovery
    # Update state for Recovery event
    state.I -= 1
    recovered = pop_random!(rng, state.currently_infected)
    validate_state(state)
    return Recovery(recovered, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               model::BirthDeathModel,
                               state::BirthDeathState, 
                               ::Type{Sampling})::Sampling
    # Update state for Sampling event
    state.I -= 1
    sampled = pop_random!(rng, state.currently_infected)
    state.n_sampled += 1
    validate_state(state)
    return Sampling(sampled, state.t)
end



function update_event_rates!(event_rates::Vector{Float64}, model::BirthDeathModel)
    λ = model.birth_rate
    μ = model.death_rate
    ψ = model.sampling_rate

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
    update_event_rates!(event_rates, model)
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