struct SIRModel <: AbstractEpiModel
    transmission_rate::Float64
    recovery_rate::Float64
    sampling_rate::Float64
    N::Int
end


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


function validate_state(state::SIRState)
    VALIDATE_STATE || return
    @assert state.S >= 0 "Susceptible individuals cannot be negative."
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


EpiState(model::SIRModel) = SIRState(0.0, model.N - 1, 1, [1], 0, 1)


get_default_stop_condition(model::SIRModel) = s -> isempty(s.currently_infected)


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


@inline function update_state!(rng::AbstractRNG,
                       state::SIRState, 
                       ::Type{Transmission})::Transmission
    # Update state for Transmission event
    state.S -= 1
    state.I += 1
    state.n_cumulative += 1
    infectee = state.n_cumulative         # Label infected individuals sequentially
    infector = sample(rng, state.currently_infected)
    push!(state.currently_infected, infectee)
    validate_state(state)
    return Transmission(infector, infectee, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                       state::SIRState, 
                       ::Type{Recovery})::Recovery
    # Update state for Recovery event
    state.I -= 1
    recovered = pop_random!(rng, state.currently_infected)
    validate_state(state)
    return Recovery(recovered, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                       state::SIRState, 
                       ::Type{Sampling})::Sampling
    # Update state for Sampling event
    state.I -= 1
    sampled = pop_random!(rng, state.currently_infected)
    state.n_sampled += 1
    validate_state(state)
    return Sampling(sampled, state.t)
end


@inline function update_event_rates!(event_rates::Vector{Float64}, 
                                     model::SIRModel, 
                                     S::Int, 
                                     I::Int)
    β = model.transmission_rate
    γ = model.recovery_rate
    ψ = model.sampling_rate
    N = model.N

    event_rates[1] = β * I * S / N  # Infection rate
    event_rates[2] = γ * I        # Recovery rate
    event_rates[3] = ψ * I        # Sampling rate
end


function simulate_outbreak(rng::AbstractRNG,
                           model::SIRModel; 
                           S_init::Int=9999, 
                           I_init::Int=1,
                           S_max::Int=100)

    # Initialize population parameters
    S = S_init
    I = I_init

    # Initialize the simulation parameters
    n_cumulative = I_init
    n_sampled = 0
    currently_infected = collect(1:I_init)

    events = Vector{SIREvent}()
    event_rates = Vector{Float64}(undef, 3)
    
    t = 0.0
    
    # Add initial infections as seeds
    for i in 1:I_init
        push!(events, Seed(i, 0.0))
    end

    while !isempty(currently_infected) && n_sampled < S_max

        update_event_rates!(event_rates, model, S, I)
        total_event_rate = sum(event_rates)

        rand_number = rand(rng)
        t -= log(rand_number) / total_event_rate

        if rand_number ≤ event_rates[1] / total_event_rate
            # Infection event
            S -= 1
            I += 1
            n_cumulative += 1
            infectee = n_cumulative
            
            # Get a random infector from the infected pool
            infector = sample(rng, currently_infected)
            transmission!(events, infector, infectee, t)
            push!(currently_infected, infectee)
        elseif rand_number ≤ (event_rates[1] + event_rates[2]) / total_event_rate
            # Recovery event
            I -= 1
            recovered = pop_random!(rng, currently_infected)
            recovery!(events, recovered, t)
        else
            # Sampling event
            I -= 1
            sampled = pop_random!(rng, currently_infected)
            sampling!(events, sampled, t)
            n_sampled += 1
        end
    end
    return events
end


function simulate_outbreak(model::SIRModel;
                           S_init::Int=9999,
                           I_init::Int=1,
                           S_max::Int=100)
    return simulate_outbreak(Random.GLOBAL_RNG, model; S_init=S_init, I_init=I_init, S_max=S_max)
end