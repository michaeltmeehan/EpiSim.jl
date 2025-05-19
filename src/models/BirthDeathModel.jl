struct BirthDeathParameters <: AbstractEpiParameters
    birth_rate::Float64
    death_rate::Float64
    sampling_rate::Float64
end

# Convenience keyword constructor
BirthDeathParameters(; birth_rate=2.0, death_rate=0.9, sampling_rate=0.1) =
    BirthDeathParameters(birth_rate, death_rate, sampling_rate)

calc_R0(par::BirthDeathParameters) = par.birth_rate / (par.death_rate + par.sampling_rate)

calc_infectious_period(par::BirthDeathParameters) = 1.0 / (par.death_rate + par.sampling_rate)

calc_sampling_fraction(par::BirthDeathParameters) = par.sampling_rate / (par.death_rate + par.sampling_rate)

calc_extinction_probability(par::BirthDeathParameters) = 1. / calc_R0(par)


function summarize(p::BirthDeathParameters)
    R₀ = p.birth_rate / (p.death_rate + p.sampling_rate)
    infectious_period = 1 / (p.death_rate + p.sampling_rate)
    return (;R₀ = R₀, infectious_period = infectious_period)
end


struct BirthDeathModel{S<:AbstractEpiState} <: AbstractEpiModel
    parameters::BirthDeathParameters
    event_types::Vector{DataType}
    initial_state::S
end


@forward BirthDeathModel.parameters summarize

calc_R0(model::BirthDeathModel) = calc_R0(model.parameters)
calc_infectious_period(model::BirthDeathModel) = calc_infectious_period(model.parameters)
calc_extinction_probability(model::BirthDeathModel) = calc_extinction_probability(model.parameters)^model.initial_state.I


const BIRTHDEATH_EVENT_TYPES = [Transmission, Recovery, Sampling]


mutable struct AgenticBirthDeathState <: AgenticState
    t::Float64
    I::Int
    currently_infected::Vector{Int}
    n_sampled::Int
    n_cumulative::Int
end

AgenticBirthDeathState(; t=0., I=1, currently_infected=collect(1:I), n_sampled=0, n_cumulative=0) =
    AgenticBirthDeathState(t, I, currently_infected, n_sampled, n_cumulative)


mutable struct AggregateBirthDeathState <: AggregateState
    t::Float64
    I::Int
end

AggregateBirthDeathState(; t=0., I=1) = AggregateBirthDeathState(t, I)


function BirthDeathModel(; birth_rate=2.0, death_rate=0.9, sampling_rate=0.1, I=1, agentic=true)
    parms = BirthDeathParameters(; birth_rate, death_rate, sampling_rate)
    state = agentic ?
        AgenticBirthDeathState(; I=I) :
        AggregateBirthDeathState(; I=I)
    return BirthDeathModel(parms, BIRTHDEATH_EVENT_TYPES, state)
end



capture(state::Union{AgenticBirthDeathState, AggregateBirthDeathState}) = (; t=state.t, I=state.I)

get_default_stop_condition(model::BirthDeathModel{AgenticBirthDeathState}) = s -> s.n_sampled >= 100

get_default_stop_condition(model::BirthDeathModel{AggregateBirthDeathState}) = s -> s.I == 0 || s.t >= 100.0 || s.I >= 10_000


function initialize_event_log(state::AgenticBirthDeathState)::Vector{AbstractEpiEvent}
    event_log = Vector{AbstractEpiEvent}()
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

    events = Vector{Union{Seed, BIRTHDEATH_EVENT_TYPES...}}()

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