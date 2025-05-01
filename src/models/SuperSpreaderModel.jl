struct SuperSpreaderModel <: AbstractEpiModel
    total_R0::Float64
    relative_transmissibility::Float64
    superspreader_fraction::Float64
    recovery_rate::Vector{Float64}
    sampling_rate::Vector{Float64}
    transmission_rate::Vector{Float64}
end


function SuperSpreaderModel(total_R0::Float64,
                            relative_transmissibility::Float64,
                            superspreader_fraction::Float64,
                            recovery_rate::Vector{Float64},
                            sampling_rate::Vector{Float64})
    transmission_rate = [relative_transmissibility, 1.] .* (recovery_rate .+ sampling_rate) .* total_R0 ./ (superspreader_fraction + (1. - superspreader_fraction) .* relative_transmissibility)
    return SuperSpreaderModel(total_R0, relative_transmissibility, superspreader_fraction, recovery_rate, sampling_rate, transmission_rate)
end


const SuperSpreader_EVENT_TYPES = [Transmission, Recovery, Sampling]

const SuperSpreaderEvent = Union{Seed, SuperSpreader_EVENT_TYPES...}

get_event_types(model::SuperSpreaderModel) = SuperSpreader_EVENT_TYPES


mutable struct SuperSpreaderState <: AbstractEpiState
    t::Float64
    I::Vector{Int}
    currently_infected::Vector{Vector{Int}}
    n_sampled::Int
    n_cumulative::Int
end


struct SuperSpreaderStateSlice <: AbstractEpiStateSlice
    t::Float64
    I_regular::Int
    I_super::Int
end


slice(state::SuperSpreaderState)::SuperSpreaderStateSlice = SuperSpreaderStateSlice(state.t, state.I[1], state.I[2])


EpiState(model::SuperSpreaderModel) = SuperSpreaderState(0.0, [0, 1], [[], [1]], 0, 1)

get_default_stop_condition(model::SuperSpreaderModel) = s -> all(isempty.(s.currently_infected)) || s.n_cumulative >= 10_000 || s.n_sampled >= 100 || s.t >= 100.0


function initialize_event_log(state::SuperSpreaderState)::Vector{SuperSpreaderEvent}
    event_log = Vector{SuperSpreaderEvent}()
    for i in 1:state.n_cumulative
        push!(event_log, Seed(i, 0.0))
    end
    return event_log
end


@inline function update_event_rates!(event_rates::Vector{Float64}, 
                                     model::SuperSpreaderModel, 
                                     state::SuperSpreaderState)
    λ = model.transmission_rate
    μ = model.recovery_rate
    ψ = model.sampling_rate
    I = state.I

    event_rates[1] = λ ⋅ I # Infection rate
    event_rates[2] = μ ⋅ I # Recovery rate
    event_rates[3] = ψ ⋅ I # Sampling rate
end


@inline function update_state!(rng::AbstractRNG,
                               model::SuperSpreaderModel,
                               state::SuperSpreaderState, 
                               ::Type{Transmission})::Transmission

    # Transmission event
    transmission_weights = model.transmission_rate .* state.I
    parent_type = wsample(rng, 1:2, transmission_weights)
    child_type = wsample(rng, 1:2, model.transmission_rate[parent_type, :])
    state.I[child_type] += 1
    state.n_cumulative += 1
    infectee = state.n_cumulative         # Label infected individuals sequentially
    infector = sample(rng, state.currently_infected[parent_type])
    push!(state.currently_infected[child_type], infectee)
    Transmission(infector, infectee, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               model::SuperSpreaderModel,
                               state::SuperSpreaderState, 
                               ::Type{Recovery})::Recovery
    # Recovery event
    recovery_weights = model.recovery_rate .* state.I
    recovery_type = wsample(rng, 1:2, recovery_weights)
    state.I[recovery_type] -= 1
    recovered = pop_random!(rng, state.currently_infected[recovery_type])
    return Recovery(recovered, state.t)
end


@inline function update_state!(rng::AbstractRNG,
                               model::SuperSpreaderModel,
                               state::SuperSpreaderState, 
                               ::Type{Sampling})::Sampling
    # Sampling event
    sampling_weights = model.sampling_rate .* state.I
    sampled_type = wsample(rng, 1:2, sampling_weights)
    state.I[sampled_type] -= 1
    state.n_sampled += 1
    sampled = pop_random!(rng, state.currently_infected[sampled_type])
    return Sampling(sampled, state.t)
end


@inline function update_event_rates!(event_rates::Vector{Float64},
                             model::SuperSpreaderModel,
                             I::Vector{Int})
    # Update event rates for the SuperSpreaderModel

    λ = model.transmission_rate
    μ = model.recovery_rate
    ψ = model.sampling_rate

    event_rates[1] = λ ⋅ I # Infection rate
    event_rates[2] = μ ⋅ I # Recovery rate
    event_rates[3] = ψ ⋅ I # Sampling rate
end


function simulate_events(rng::AbstractRNG,
                           model::SuperSpreaderModel;
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
            # Transmission event
            transmission_weights = model.transmission_rate .* I
            parent_type = wsample(rng, 1:2, transmission_weights)
            child_type = wsample(rng, 1:2, model.transmission_rate[parent_type, :])
            I[child_type] += 1
            n_cumulative += 1
            infectee = n_cumulative         # Label infected individuals sequentially
            infector = sample(rng, currently_infected[parent_type])
            push!(currently_infected[child_type], infectee)
            transmission!(events, infector, infectee, t)
        elseif rand_number ≤ (event_rates[1] + event_rates[2]) / total_event_rate
            # Recovery event
            recovery_weights = model.recovery_rate .* I
            recovery_type = wsample(rng, 1:2, recovery_weights)
            I[recovery_type] -= 1
            recovered = pop_random!(rng, currently_infected[recovery_type])
            recovery!(events, recovered, t)
        else
            # Sampling event
            sampling_weights = model.sampling_rate .* I
            sampled_type = wsample(rng, 1:2, sampling_weights)
            I[sampled_type] -= 1
            n_sampled += 1
            sampled = pop_random!(rng, currently_infected[sampled_type])
            sampling!(events, sampled, t)
        end
    end
    return events
end


function simulate_events(model::SuperSpreaderModel; 
                           I_init::Vector{Int}=[1, 0],
                           N_max::Int=10_000, 
                           S_max::Int=100)
    return simulate_events(Random.GLOBAL_RNG, model, I_init=I_init, N_max=N_max, S_max=S_max)
end