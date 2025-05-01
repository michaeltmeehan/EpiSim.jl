struct MultiTypeBirthDeathModel <: AbstractEpiModel
    birth_rate::Matrix{Float64} # Birth rate matrix (n_types x n_types): element (i,j) is the rate individuals in state i give birth to individuals in state j
    death_rate::Vector{Float64}
    sampling_rate::Vector{Float64}
    n_types::Int
    Λ::Vector{Float64} # Total birth rate for each type
end


function MultiTypeBirthDeathModel(birth_rate::Matrix{Float64}, 
                                  death_rate::Vector{Float64}, 
                                  sampling_rate::Vector{Float64})::MultiTypeBirthDeathModel
    typeof(death_rate)
    n_types = size(birth_rate, 1)
    Λ = sum(birth_rate, dims=2)[:]  # Total birth rate for each type
    return MultiTypeBirthDeathModel(birth_rate, death_rate, sampling_rate, n_types, Λ)
end


const MultiTypeBirthDeath_EVENT_TYPES = [Transmission, Recovery, Sampling]

const MultiTypeBirthDeathEvent = Union{Seed, MultiTypeBirthDeath_EVENT_TYPES...}

get_event_types(model::MultiTypeBirthDeathModel) = MultiTypeBirthDeath_EVENT_TYPES


mutable struct MultiTypeBirthDeathState <: AbstractEpiState
    t::Float64
    I::Vector{Int}
    currently_infected::Vector{Vector{Int}}
    n_sampled::Int
    n_cumulative::Int
end


EpiState(model::MultiTypeBirthDeathModel) = MultiTypeBirthDeathState(0.0, [1, 0], [[1], []], 0, 1)

get_default_stop_condition(model::MultiTypeBirthDeathModel) = s -> all(isempty.(s.currently_infected)) || s.n_cumulative >= 10_000 || s.n_sampled >= 100 || s.t >= 100.0


function initialize_event_log(state::MultiTypeBirthDeathState)::Vector{MultiTypeBirthDeathEvent}
    event_log = Vector{MultiTypeBirthDeathEvent}()
    for i in 1:state.n_cumulative
        push!(event_log, Seed(i, 0.0))
    end
    return event_log
end


@inline function update_event_rates!(event_rates::Vector{Float64}, 
                                     model::MultiTypeBirthDeathModel, 
                                     state::MultiTypeBirthDeathState)
    Λ = model.Λ
    μ = model.death_rate
    ψ = model.sampling_rate
    I = state.I

    event_rates[1] = Λ ⋅ I # Infection rate
    event_rates[2] = μ ⋅ I # Recovery rate
    event_rates[3] = ψ ⋅ I # Sampling rate
end


@inline function update_state!(rng::AbstractRNG,
                               model::MultiTypeBirthDeathModel,
                               state::MultiTypeBirthDeathState, 
                               ::Type{Transmission})::Transmission

    # Transmission event
    transmission_weights = model.Λ .* state.I
    parent_type = wsample(rng, 1:model.n_types, transmission_weights)
    child_type = wsample(rng, 1:model.n_types, model.Λ[parent_type, :])
    state.I[child_type] += 1
    state.n_cumulative += 1
    infectee = state.n_cumulative         # Label infected individuals sequentially
    infector = sample(rng, state.currently_infected[parent_type])
    push!(state.currently_infected[child_type], infectee)
    Transmission(infector, infectee, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               model::MultiTypeBirthDeathModel,
                               state::MultiTypeBirthDeathState, 
                               ::Type{Recovery})::Recovery
    # Recovery event
    recovery_weights = model.death_rate .* state.I
    recovery_type = wsample(rng, 1:model.n_types, recovery_weights)
    state.I[recovery_type] -= 1
    recovered = pop_random!(rng, state.currently_infected[recovery_type])
    return Recovery(recovered, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               model::MultiTypeBirthDeathModel,
                               state::MultiTypeBirthDeathState, 
                               ::Type{Sampling})::Sampling
    # Sampling event
    sampling_weights = model.sampling_rate .* state.I
    sampled_type = wsample(rng, 1:model.n_types, sampling_weights)
    state.I[sampled_type] -= 1
    state.n_sampled += 1
    sampled = pop_random!(rng, state.currently_infected[sampled_type])
    return Sampling(sampled, state.t)
end


function update_event_rates!(event_rates::Vector{Float64}, model::MultiTypeBirthDeathModel, I::Vector{Int})
    Λ = model.Λ
    μ = model.death_rate
    ψ = model.sampling_rate

    event_rates[1] = Λ ⋅ I  # Birth rate
    event_rates[2] = μ ⋅ I  # Death rate
    event_rates[3] = ψ ⋅ I  # Sampling rate
end


function simulate_events(rng::AbstractRNG,
                           model::MultiTypeBirthDeathModel;
                           I_init::Vector{Int}=[1, 0],
                           N_max::Int=10_000,
                           S_max::Int=100)

    I = I_init
    n_cumulative = sum(I_init)
    n_sampled = 0

    currently_infected = [collect(1:I_init[type]) for type in eachindex(I_init)]

    events = Vector{AbstractEpiEvent}()

    # Pre-calculate event rates
    event_rates = Vector{Float64}(undef, 3)
    
    t = 0.0
    
    # Add initial infections as seeds
    for type in eachindex(I_init)
        for i in 1:I_init[type]
            push!(events, Seed(i, 0.0))
        end
    end

    while !all(isempty.(currently_infected)) && n_cumulative < N_max && n_sampled < S_max
        
        update_event_rates!(event_rates, model, I)
        total_event_rate = sum(event_rates)
        
        rand_number = rand(rng)
        t -= log(rand_number) / total_event_rate

        if rand_number ≤ event_rates[1] / total_event_rate
            # Birth event
            birth_weights = model.Λ .* I
            parent_type = wsample(rng, 1:model.n_types, birth_weights)
            child_type = wsample(rng, 1:model.n_types, model.birth_rate[parent_type, :])
            I[child_type] += 1
            n_cumulative += 1
            infectee = n_cumulative         # Label infected individuals sequentially
            infector = sample(rng, currently_infected[parent_type])
            push!(currently_infected[child_type], infectee)
            transmission!(events, infector, infectee, t)
        elseif rand_number ≤ (event_rates[1] + event_rates[2]) / total_event_rate
            # Death event
            death_type = wsample(rng, 1:model.n_types, I)
            I[death_type] -= 1
            recovered = pop_random!(rng, currently_infected[death_type])
            recovery!(events, recovered, t)
        else
            # Sampling event
            sampled_type = wsample(rng, 1:model.n_types, I)
            I[sampled_type] -= 1
            n_sampled += 1
            sampled = pop_random!(rng, currently_infected[sampled_type])
            sampling!(events, sampled, t)
        end
    end
    return events
end


function simulate_events(model::MultiTypeBirthDeathModel; 
                           I_init::Vector{Int}=[1, 0],
                           N_max::Int=10_000, 
                           S_max::Int=100)
    return simulate_events(Random.GLOBAL_RNG, model, I_init=I_init, N_max=N_max, S_max=S_max)
end