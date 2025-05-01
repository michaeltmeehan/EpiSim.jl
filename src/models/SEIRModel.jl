struct SEIRModel <: AbstractEpiModel
    transmission_rate::Float64
    activation_rate::Float64
    recovery_rate::Float64
    sampling_rate::Float64
    N::Int
end


const SEIR_EVENT_TYPES = [Transmission, Activation, Recovery, Sampling]

const SEIREvent = Union{Seed, SEIR_EVENT_TYPES...}

get_event_types(model::SEIRModel) = SEIR_EVENT_TYPES


mutable struct SEIRState <: AbstractEpiState
    t::Float64
    S::Int
    E::Int
    I::Int
    currently_exposed::Vector{Int}
    currently_infected::Vector{Int}
    n_sampled::Int
    n_cumulative::Int
end


EpiState(model::SEIRModel) = SEIRState(0.0, model.N - 1, 0, 1, [], [1], 0, 1)


get_default_stop_condition(model::SEIRModel) = s -> isempty(s.currently_infected)


function initialize_event_log(state::SEIRState)::Vector{SEIREvent}
    event_log = Vector{SEIREvent}()
    for i in 1:state.I
        push!(event_log, Seed(i, 0.0))
    end
    return event_log
end


@inline function update_event_rates!(event_rates::Vector{Float64}, 
                                     model::SEIRModel, 
                                     state::SEIRState)
    β = model.transmission_rate
    α = model.activation_rate
    γ = model.recovery_rate
    ψ = model.sampling_rate
    N = model.N

    S = state.S
    E = state.E
    I = state.I

    event_rates[1] = β * I * S / N  # Infection rate
    event_rates[2] = α * E        # Activation rate
    event_rates[3] = γ * I        # Recovery rate
    event_rates[4] = ψ * I        # Sampling rate
end


@inline function update_state!(rng::AbstractRNG,
                               model::SEIRModel,
                               state::SEIRState, 
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
                               model::SEIRModel,
                               state::SEIRState, 
                               ::Type{Activation})::Activation
    # Update state for Activation event
    state.E -= 1
    state.I += 1
    activated = pop_random!(rng, state.currently_exposed)
    push!(state.currently_infected, activated)
    return Activation(activated, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               model::SEIRModel,
                               state::SEIRState, 
                               ::Type{Recovery})::Recovery
    # Update state for Recovery event
    state.I -= 1
    recovered = pop_random!(rng, state.currently_infected)
    return Recovery(recovered, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               model::SEIRModel,
                               state::SEIRState, 
                               ::Type{Sampling})::Sampling
    # Update state for Sampling event
    state.I -= 1
    sampled = pop_random!(rng, state.currently_infected)
    state.n_sampled += 1
    return Sampling(sampled, state.t)
end


function update_event_rates!(event_rates::Vector{Float64}, model::SEIRModel, S::Int, E::Int, I::Int)
    β = model.transmission_rate
    ν = model.activation_rate
    γ = model.recovery_rate
    ψ = model.sampling_rate
    N = model.N

    event_rates[1] = β * I * S / N  # Infection rate
    event_rates[2] = ν * E        # Activation rate
    event_rates[3] = γ * I        # Recovery rate
    event_rates[4] = ψ * I        # Sampling rate
end


function simulate_events(rng::AbstractRNG,
                           model::SEIRModel;
                           S_init::Int=9999, 
                           E_init::Int=0,
                           I_init::Int=1,
                           S_max::Int=100)

    # Initialize population parameters
    S = S_init
    E = E_init
    I = I_init
     
    # Initialize the simulation parameters
    n_cumulative = E_init + I_init
    n_sampled = 0
    currently_exposed = Vector{Int}(undef, E_init)
    if E_init > 0
        currently_exposed[1:E_init] .= 1:E_init
    end
    
    currently_infected = Vector{Int}(undef, I_init)
    if I_init > 0
        currently_infected[1:I_init] .= (E_init + 1):(E_init + I_init)
    end

    events = Vector{AbstractEpiEvent}()
    event_rates = Vector{Float64}(undef, 4)
    
    t = 0.0

    # Add initial infections as seeds
    for i in 1:E_init
        push!(events, Seed(i, 0.0))
    end
    
    for i in (E_init + 1):(E_init + I_init)
        push!(events, Seed(i, 0.0))
    end

    while !isempty(currently_infected) && n_sampled < S_max

        update_event_rates!(event_rates, model, S, E, I)
        total_event_rate = sum(event_rates)

        rand_number = rand(rng)
        t -= log(rand_number) / total_event_rate

        if rand_number ≤ event_rates[1] / total_event_rate
            # Infection event (S -> E)
            S -= 1
            E += 1
            n_cumulative += 1
            infectee = n_cumulative         # Label infected individuals sequentially
            
            # Get a random infector from the infected pool
            infector = sample(rng, currently_infected)
            transmission!(events, infector, infectee, t)
            push!(currently_exposed, infectee)
        elseif rand_number ≤ (event_rates[1] + event_rates[2]) / total_event_rate
            # Activation event (E -> I)
            E -= 1
            I += 1
            activated = pop_random!(rng, currently_exposed)
            activation!(events, activated, t)
            push!(currently_infected, activated)
        elseif rand_number ≤ (event_rates[1] + event_rates[2] + event_rates[3]) / total_event_rate
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


function simulate_events(model::SEIRModel;
                           S_init::Int=9999,
                           E_init::Int=0,
                           I_init::Int=1,
                           S_max::Int=100)
    return simulate_events(Random.GLOBAL_RNG, model; S_init=S_init, E_init=E_init, I_init=I_init, S_max=S_max)
end