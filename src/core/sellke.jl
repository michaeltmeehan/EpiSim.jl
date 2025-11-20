_dirac(x::Float64) = DiscreteNonParametric([x], [1.0])
_to_dist(x::Union{Float64, Distribution}, exp_default=false) = x isa Distribution ? x : (exp_default ? Exponential(x) : _dirac(x))


function sellke(I₀::Int, 
                S₀::Int, 
                β::Union{Float64, Distribution}, 
                τ::Union{Float64, Distribution})

    @assert I₀ ≥ 0 "Number of initial infected I₀ must be at least 0"
    @assert S₀ ≥ 0 "Number of susceptibles S₀ must be at least 0"

    # Initialize number of infected, susceptible and recovered
    I = I₀; S = S₀; R = 0; N = I₀ + S₀

    Λ = 0.0     # Cumulative infection pressure
    t = 0.0     # Current time

    # Draw resistances for susceptibles (Exp(1))
    resistances = BinaryMinHeap{Float64}(randexp(S₀))

    # Convert β and τ if necessary
    dβ = _to_dist(β)
    dτ = _to_dist(τ, true)

    # Draw transmission rate and recovery times for initial infecteds
    recovery_times = t .+ rand(dτ, I₀)
    transmission_rates = rand(dβ, I₀)
    recoveries = BinaryMinHeap{Tuple{Float64,Float64}}(collect(zip(recovery_times, transmission_rates)))
    slope = sum(transmission_rates) / N

    while I > 0

        # Work out whether infection or recovery occurs next
        if !isempty(resistances) && top(resistances) ≤ Λ + slope * (top(recoveries)[1] - t) # Infection event

            # Update number of infected and susceptibles
            I += 1; S -= 1

            # Update time
            t += (top(resistances) - Λ) / slope

            # Update cumulative infection pressure (and increment heap)
            Λ = pop!(resistances)

            # Draw recovery time and transmission rate for newly infected
            τ = rand(dτ)
            βᵢ = rand(dβ)

            # Add to recoveries heap
            push!(recoveries, (t + τ, β))

            # Update slope
            slope += (βᵢ / N)

        else    # Recovery event

            # Update number of infected and recovered
            I -= 1; R += 1

            # Update cumulative infection pressure
            Λ += slope * (top(recoveries)[1] - t)

            # Update time (and increment heap)
            t, βᵢ = pop!(recoveries)

            # Decrease slope
            slope -= (βᵢ / N)
        end
    end
    return (S=S, I=I, R=R)
end


# mutable struct Agent
#     id::Int
#     state::State
#     resistance::Float64
#     transmission_rate::Float64
#     incubation_period::Float64
#     infectious_period::Float64
#     sampling_period::Float64
#     remove_on_sample::Bool
# end


# Base.isless(a1::Agent, a2::Agent) = a1.resistance < a2.resistance


abstract type Agent end


mutable struct Susceptible <: Agent
    id::Int
    resistance::Float64
end

Base.isless(s1::Susceptible, s2::Susceptible) = s1.resistance < s2.resistance


struct Traits
    transmission_rate::Float64
    incubation_period::Float64
    infectious_period::Float64
    sampling_period::Float64
    remove_on_sample::Bool
end


mutable struct Infected <: Agent
    id::Int
    traits::Traits
    events::Vector{Event}
    head::Int
end

isdone(infected::Infected) = infected.head > length(infected.events)
next_event(infected::Infected) = !isdone(infected) ? infected.events[infected.head] : nothing
next_time(agent::Infected) = agent.events[agent.head].time
Base.isless(i1::Infected, i2::Infected) = next_time(i1) < next_time(i2)


isfossil(infected::Infected) = infected.traits.remove_on_sample


struct TraitDists
    dβ::Distribution
    dτₑ::Distribution
    dτᵢ::Distribution
    dτₛ::Distribution
end


function TraitDists(β::Union{Float64, Distribution},
                    τₑ::Union{Float64, Distribution},
                    τᵢ::Union{Float64, Distribution},
                    τₛ::Union{Float64, Distribution})
    return TraitDists(_to_dist(β), _to_dist(τₑ, true), _to_dist(τᵢ, true), _to_dist(τₛ, true))
end


Base.rand(td::TraitDists, r::Float64) = Traits(rand(td.dβ), rand(td.dτₑ), rand(td.dτᵢ), rand(td.dτₛ), rand() < r)


# TODO: In principle, infecteds could be sampled multiple times before recovery / removal (run loop until either sampling results in recovery or sampling time > recovery time)
function make_infected(t::Float64, id::Int, td::TraitDists, r::Float64; exposed::Bool)
    tr = rand(td, r)    # Generate traits

    # Absolute times
    t_EI  = exposed ? t + tr.incubation_period : t            # E→I (activation) or already infectious
    t_rec = t_EI + tr.infectious_period
    t_smp = t_EI + tr.sampling_period

    events = Event[]
    if exposed
        push!(events, Activation(t_EI, id))
    end

    if tr.remove_on_sample && t_smp < t_rec
        # sampling causes removal; no recovery event needed
        push!(events, SerialSampling(t_smp, id))
    else
        # optional non-removal sampling (only if it occurs before recovery)
        if t_smp < t_rec
            push!(events, FossilizedSampling(t_smp, id))
        end
        push!(events, Recovery(t_rec, id))
    end

    return Infected(id, tr, events, 1)
end


function sellke(S₀::Int, 
                E₀::Int,
                I₀::Int, 
                β::Union{Float64, Distribution},    # Transmission rate
                τₑ::Union{Float64, Distribution},   # Incubation period
                τᵢ::Union{Float64, Distribution},   # Infectious period
                τₛ::Union{Float64, Distribution},    # Sampling time
                r::Float64                           # Probability of removal upon sampling
                )

    # Initialize number of exposed, infected, susceptible and recovered
    S, E, I, R = S₀, E₀, I₀, 0
    
    # Total population size
    N = E₀ + I₀ + S₀

    # Initialize time and cumulative infection pressure
    t = Λ = 0.0

    # Package parameters as distributions of traits
    td = TraitDists(β, τₑ, τᵢ, τₛ)

    # Initialize vectors of exposed and infective individuals
    exposed = [make_infected(t, id, td, r; exposed=true) for id in 1:E₀]
    infectives = [make_infected(t, id, td, r; exposed=false) for id in (E₀ + 1):(E₀ + I₀)]

    # Initialize ordered heap of susceptible individuals
    susceptibles = BinaryMinHeap{Susceptible}([Susceptible(id, randexp()) for id in (E₀ + I₀ + 1):N])

    # Initialize ordered heap of infected individuals (i.e., those who are exposed or infective)
    infecteds = BinaryMinHeap{Infected}(vcat(exposed, infectives))

    # Initialize list of transmission rates for sampling infectors
    transmission_rates = fill(0.0, N)
    for infective in infectives
        transmission_rates[infective.id] = infective.traits.transmission_rate
    end

    # Calculate initial slope (i.e., cumulative transmission rate)
    dΛ = sum(transmission_rates) / N

    # Initialize event log
    events = Event[Seeding(0.0, id) for id in 1:(E₀ + I₀)]

    # Initialize state log
    states = fill(State(t, S, E, I, R), E₀ + I₀)

    while E + I > 0

        # Look up time for next infected event
        t_next = next_time(top(infecteds))

        # Work out whether infection or temporal event occurs next
        if !isempty(susceptibles) && top(susceptibles).resistance ≤ Λ + dΛ * (t_next - t) # Infection event

            # Remove susceptible from heap
            susceptible = pop!(susceptibles)

            # Create new infected individual
            new_infectee = make_infected(t, susceptible.id, td, r; exposed=true)

            # Add new infected to heap
            push!(infecteds, new_infectee)

            # Update time
            t += (susceptible.resistance - Λ) / dΛ

            # Update number of exposed and susceptibles
            E += 1; S -= 1

            # Update cumulative infection pressure
            Λ = susceptible.resistance

            # Update event log
            push!(events, Transmission(t, new_infectee.id, wsampleindex(transmission_rates)))

        else    # Temporal event (activation, recovery, sampling)
            # Withdraw infected from heap
            infected = pop!(infecteds)

            # Retrieve next event for this infected
            event = next_event(infected)

            # Update event log
            push!(events, event)

            # Update cumulative infection pressure and time
            Λ += dΛ * (t_next - t)
            t = t_next

            # Advance head to next event
            infected.head += 1

            if event isa Activation   # Activation event
                # Update number of exposed and infected
                E -= 1; I += 1

                # Activate their transmission rate
                transmission_rates[infected.id] = infected.traits.transmission_rate

                # Increment slope
                dΛ += infected.traits.transmission_rate / N

            elseif event isa SerialSampling || event isa Recovery   # Sampling with removal or recovery event
                # Update number of infected and recovered
                I -= 1; R += 1

                # Remove their transmission rate
                transmission_rates[infected.id] = 0.0

                # Decrease slope
                dΛ -= (infected.traits.transmission_rate / N)

            end

            # If events remain for this infected, push back onto heap
            !isdone(infected) && push!(infecteds, infected)

        end

        # Push new state to log
        push!(states, State(t, S, E, I, R))
    end
    return Simulation(states, events, 0)
end