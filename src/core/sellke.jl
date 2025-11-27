_dirac(x::Float64) = DiscreteNonParametric([x], [1.0])
_to_dist(x::Union{Float64, Distribution}, exp_default=false) = x isa Distribution ? x : (exp_default ? Exponential(x) : _dirac(x))


@enum StateKind::UInt8 begin
    SK_None = 0
    SK_Susceptible = 1
    SK_Exposed = 2
    SK_Infected = 3
    SK_Removed = 4
end


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


mutable struct Population
    kind::Vector{StateKind}   # 0: None, 1: Susceptible, 2: Exposed, 3: Infected, 4: Removed
    resistance::Vector{Float64}   # only for Susceptible; Inf otherwise
    transmission_rate::Vector{Float64}   # only for Infected; 0.0 otherwise
    incubation_period::Vector{Float64}   # only for Exposed; 0.0 otherwise
    infectious_period::Vector{Float64}   # only for Infected; 0.0 otherwise
    sampling_period::Vector{Float64}   # only for Infected; 0.0 otherwise
    t_infection::Vector{Float64}   # only for Exposed and Infected; NaN otherwise
    t_activation::Vector{Float64}   # only for Exposed; NaN otherwise
    t_recovery::Vector{Float64}   # only for Infected; NaN
    t_sampling::Vector{Float64}   # only for Infected; NaN otherwise
    next_event::Vector{EventKind}
    next_event_time::Vector{Float64}
end


function Population(N::Int)
    return Population(fill(SK_Susceptible, N),
                      randexp(N),
                      zeros(Float64, N),
                      zeros(Float64, N),
                      zeros(Float64, N),
                      zeros(Float64, N),
                      fill(NaN, N),
                      fill(NaN, N),
                      fill(NaN, N),
                      fill(NaN, N),
                      fill(EK_None, N),
                      fill(Inf, N)
                      )
end


function schedule_host!(pop::Population,
                        i::Int,
                        t::Float64,
                        td::TraitDists,
                        r::Float64;
                        exposed::Bool)

    pop.transmission_rate[i] = rand(td.dβ)
    pop.incubation_period[i] = rand(td.dτₑ)
    pop.infectious_period[i] = rand(td.dτᵢ)
    pop.sampling_period[i] = rand(td.dτₛ)

    pop.t_infection[i] = t
    if exposed
        pop.kind[i] = SK_Exposed
        pop.next_event[i] = EK_Activation
        pop.next_event_time[i] = pop.t_activation[i] = t + pop.incubation_period[i]
    else
        pop.kind[i] = SK_Infected
        pop.t_activation[i] = t
        if pop.infectious_period[i] < pop.sampling_period[i]
            pop.next_event[i] = EK_Recovery
            pop.next_event_time[i] = t + pop.infectious_period[i]
        elseif rand() < r
            pop.next_event[i] = EK_SerialSampling
            pop.next_event_time[i] = t + pop.sampling_period[i]
        else
            pop.next_event[i] = EK_FossilizedSampling
            pop.next_event_time[i] = t + pop.sampling_period[i]
        end
    end
end


function advance_host!(pop::Population, id::Int, td::TraitDists, r::Float64)
    event_kind = pop.next_event[id]
    t = pop.next_event_time[id]
    if event_kind == EK_Activation
        pop.kind[id] = SK_Infected
        pop.t_activation[id] = t
        if pop.infectious_period[id] < pop.sampling_period[id]
            pop.next_event[id] = EK_Recovery
            pop.next_event_time[id] = pop.t_activation[id] + pop.infectious_period[id]
        elseif rand() < r
            pop.next_event[id] = EK_SerialSampling
            pop.next_event_time[id] = pop.t_activation[id] + pop.sampling_period[id]
        else
            pop.next_event[id] = EK_FossilizedSampling
            pop.next_event_time[id] = pop.t_activation[id] + pop.sampling_period[id]
        end
    elseif event_kind == EK_SerialSampling || event_kind == EK_Recovery
        pop.kind[id] = SK_Removed
        pop.transmission_rate[id] = 0.0
        pop.next_event[id] = EK_None
        pop.next_event_time[id] = Inf
        if event_kind == EK_SerialSampling
            pop.t_sampling[id] = t
        else
            pop.t_recovery[id] = t
        end
    elseif event_kind == EK_FossilizedSampling
        pop.t_sampling[id] = t
        # Draw a new sampling interval (gap) and update time since activation
        Δτₛ = rand(td.dτₛ)
        pop.sampling_period[id] += Δτₛ  # now: time since activation for next sample

        if pop.sampling_period[id] < pop.infectious_period[id]
            # still within infectious window
            if rand() < r
                pop.next_event[id] = EK_SerialSampling
            else
                pop.next_event[id] = EK_FossilizedSampling
            end
            pop.next_event_time[id] =
                pop.t_activation[id] + pop.sampling_period[id]
        else
            pop.next_event[id]      = EK_Recovery
            pop.next_event_time[id] =
                pop.t_activation[id] + pop.infectious_period[id]
        end
    elseif event_kind == EK_None
        # Do nothing
    else
        error("Unknown event kind: $event_kind")
    end
end


# susceptibility heap: (resistance, host_index)
const SusKey = Tuple{Float64,Int}
# infected heap: (next_event_time, host_index)
const InfKey = Tuple{Float64,Int}


# TODO: Add likelihood calculation
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

    # Initialize population
    pop = Population(N)

    # Assign resistances to susceptibles
    susceptibles = BinaryMinHeap{SusKey}([(pop.resistance[id], id) for id in (E₀ + I₀ + 1):N])

    # Initialize exposed individuals
    for id in 1:E₀
        schedule_host!(pop, id, t, td, r; exposed=true)
    end

    # Initialize infected individuals
    for id in (E₀ + 1):(E₀ + I₀)
        schedule_host!(pop, id, t, td, r; exposed=false)
    end

    infecteds = BinaryMinHeap{InfKey}([(pop.next_event_time[id], id) for id in 1:(E₀ + I₀)])

    # Calculate initial slope (i.e., cumulative transmission rate)
    dΛ = sum(pop.transmission_rate[(E₀+1):(E₀ + I₀)]) / N

    # Initialize event log components with seeding events
    times = zeros(Float64, E₀+I₀)
    hosts = collect(1:(E₀+I₀))
    infectors = zeros(Int, E₀+I₀)
    kinds = fill(EK_Seeding, E₀+I₀)    

    while E + I > 0

        # Look up time for next infected event
        t_next = top(infecteds)[1]

        # Work out whether infection or temporal event occurs next
        if !isempty(susceptibles) && top(susceptibles)[1] ≤ Λ + dΛ * (t_next - t) # Infection event

            # Remove susceptible from heap
            resistance, id = pop!(susceptibles)

            # Update time
            t += (resistance - Λ) / dΛ

            # Update number of exposed and susceptibles
            E += 1; S -= 1

            # Update cumulative infection pressure
            Λ = resistance

            # Create new exposed individual
            schedule_host!(pop, id, t, td, r; exposed=true)

            # Add new infected to heap
            push!(infecteds, (pop.next_event_time[id], id))

            # Update event log
            push!(times, t)
            push!(hosts, id)
            push!(infectors, wsampleindex(pop.transmission_rate))   # TODO: assign infector ID
            push!(kinds, EK_Transmission)

        else    # Temporal event (activation, recovery, sampling)
            # Withdraw infected from heap
            t_next, id = pop!(infecteds)

            # Retrieve scheduled event for this infected
            scheduled_event = pop.next_event[id]

            # Update event log
            push!(times, t_next)
            push!(hosts, id)
            push!(infectors, 0)
            push!(kinds, scheduled_event)

            # Update cumulative infection pressure and time
            Λ += dΛ * (t_next - t)
            t = t_next

            if scheduled_event == EK_Activation   # Activation event
                # Update number of exposed and infected
                E -= 1; I += 1

                # Increment slope
                dΛ += pop.transmission_rate[id] / N

            elseif scheduled_event == EK_SerialSampling || scheduled_event == EK_Recovery   # Sampling with removal or recovery event
                # Update number of infected and recovered
                I -= 1; R += 1

                # Decrease slope
                dΛ -= (pop.transmission_rate[id] / N)

            end

            # Advance host state
            advance_host!(pop, id, td, r)

            # If events remain for this infected, push back onto heap
            pop.next_event[id] != EK_None && push!(infecteds, (pop.next_event_time[id], id))

        end
    end
    return EventLog(times, hosts, infectors, kinds)
end