_dirac(x::Float64) = DiscreteNonParametric([x], [1.0])
_to_dist(x::Union{Float64, Distribution}, exp_default=false) = x isa Distribution ? x : (exp_default ? Exponential(x) : _dirac(x))


@enum StateKind::UInt8 begin
    SK_None = 0
    SK_Susceptible = 1
    SK_Exposed = 2
    SK_Infected = 3
    SK_Removed = 4
end
# State kinds:
# SK_Susceptible: in susceptibles heap, resistance[i] < ∞, β[i] = 0
# SK_Exposed: in infecteds heap, β[i] = 0, next_event = EK_Activation
# SK_Infected: in infecteds heap, β[i] > 0, next_event ∈ {EK_FossilizedSampling, EK_SerialSampling, EK_Recovery}
# SK_Removed: never in any heap, β[i] = 0, next_event = EK_None


const STATEKIND_LABELS = Dict(
    SK_None              => "None",
    SK_Susceptible       => "Susceptible",
    SK_Exposed          => "Exposed",
    SK_Infected         => "Infected",
    SK_Removed          => "Removed"
)


function Base.show(io::IO, x::StateKind)
    print(io, STATEKIND_LABELS[x])
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
    β::Vector{Float64}   # transmission rate; only for Infected; 0.0 otherwise
    τₑ::Vector{Float64}   # incubation period; only for Exposed; 0.0 otherwise
    τᵢ::Vector{Float64}   # infectious period; only for Infected; 0.0 otherwise
    τₛ::Vector{Float64}   # time since activation to next sample; only for Infected; 0.0 otherwise
    t_infection::Vector{Float64}   # only for Exposed and Infected; NaN otherwise
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
                      fill(EK_None, N),
                      fill(Inf, N)
                      )
end


mutable struct ActiveList
    ids::Vector{Int}
    β::Vector{Float64}
    pos_in_active::Vector{Int}  # position of each host in active_ids and active_β; 0 if not active
    total_β::Float64
end


function ActiveList(N::Int)
    return ActiveList(Int[], Float64[], zeros(Int, N), 0.0)
end


function add_active!(actives::ActiveList,
                     id::Int,
                     β::Float64)

    push!(actives.ids, id)
    push!(actives.β, β)
    actives.pos_in_active[id] = length(actives.ids)
    actives.total_β += β
    return
end


function remove_active!(actives::ActiveList,
                        id::Int)

    actives.total_β -= actives.β[actives.pos_in_active[id]]
    idx = actives.pos_in_active[id]
    @assert idx != 0 "remove_active! called for non-active host $id"
    last_id = actives.ids[end]
    last_β = actives.β[end]

    # Swap with last element
    actives.ids[idx] = last_id
    actives.β[idx] = last_β
    actives.pos_in_active[last_id] = idx

    # Pop last element
    pop!(actives.ids)
    pop!(actives.β)
    actives.pos_in_active[id] = 0
    return
end


function expose_host!(pop::Population,
                      id::Int,
                      t::Float64,
                      td::TraitDists)

    @unpack kind, resistance, β, τₑ, τᵢ, τₛ, 
            t_infection, next_event, next_event_time = pop

    t_infection[id] = t
    kind[id] = SK_Exposed
    τₑ[id] = rand(td.dτₑ)
    next_event[id] = EK_Activation
    next_event_time[id] = t + τₑ[id]
end


function schedule_sampling(τₑ::Float64,
                           τᵢ::Float64,
                           τₛ::Float64,
                           t_infection::Float64,
                           r::Float64)
    t_activation = t_infection + τₑ
    if τₛ < τᵢ
        # Sampling within infectious period
        if rand() < r
            next_event = EK_SerialSampling
        else
            next_event = EK_FossilizedSampling
        end
        next_event_time = t_activation + τₛ
    else
        # Sampling after infectious period
        next_event = EK_Recovery
        next_event_time = t_activation + τᵢ
    end
    return next_event, next_event_time
end


function activate_host!(pop::Population,
                        id::Int,
                        td::TraitDists,
                        r::Float64,
                        active::ActiveList)

    @unpack kind, resistance, β, τₑ, τᵢ, τₛ, 
            t_infection, next_event, next_event_time = pop

    kind[id] = SK_Infected
    β[id] = rand(td.dβ)
    τᵢ[id] = rand(td.dτᵢ)
    τₛ[id] = rand(td.dτₛ)
    next_event[id], next_event_time[id] = schedule_sampling(τₑ[id], τᵢ[id], τₛ[id], t_infection[id], r)

    add_active!(active, id, β[id])
    return
end


function fossilize_host!(pop::Population,
                        id::Int,
                        td::TraitDists,
                        r::Float64)

    @unpack kind, resistance, β, τₑ, τᵢ, τₛ, 
            t_infection, next_event, next_event_time = pop

    # Draw a new sampling interval (gap) and update time since activation
    Δτₛ = rand(td.dτₛ)
    τₛ[id] += Δτₛ  # now: time since activation for next sample
    next_event[id], next_event_time[id] = schedule_sampling(τₑ[id], τᵢ[id], τₛ[id], t_infection[id], r)
end


function deactivate_host!(pop::Population,
                        id::Int,
                        actives::ActiveList)

    @unpack kind, resistance, β, τₑ, τᵢ, τₛ, 
            t_infection, next_event, next_event_time = pop

    kind[id] = SK_Removed
    actives.total_β -= β[id]
    β[id] = 0.0
    next_event[id] = EK_None
    next_event_time[id] = Inf

    remove_active!(actives, id)
    return
end


struct SusKey
    resistance::Float64
    id::Int
end

Base.isless(a::SusKey, b::SusKey) = a.resistance < b.resistance

struct InfKey
    next_event_time::Float64
    id::Int
end

Base.isless(a::InfKey, b::InfKey) = a.next_event_time < b.next_event_time


function sample_infector(actives::ActiveList)
    @assert actives.total_β > 0.0 "sample_infector called with total_β = 0"
    u = rand() * actives.total_β
    cumulative = 0.0
    @inbounds for idx in eachindex(actives.β)
        cumulative += actives.β[idx]
        if cumulative ≥ u
            return actives.ids[idx]
        end
    end
    return actives.ids[end]  # Fallback
end


# TODO: Add likelihood calculation
# TODO: Implement in-place allocations of population
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
    susceptibles = BinaryMinHeap{SusKey}([SusKey(pop.resistance[id], id) for id in (E₀ + I₀ + 1):N])

    # Initialize exposed individuals
    for id in 1:E₀
        expose_host!(pop, id, t, td)
    end

    # Initialize list of active infectors
    actives = ActiveList(N)

    # Initialize infected individuals
    for id in (E₀ + 1):(E₀ + I₀)
        pop.t_infection[id] = t
        activate_host!(pop, id, td, r, actives)
    end

    # Initialize infected event heap
    infecteds = BinaryMinHeap{InfKey}([InfKey(pop.next_event_time[id], id) for id in 1:(E₀ + I₀)])

    # Calculate initial slope (i.e., cumulative transmission rate)
    dΛ = actives.total_β / N

    # Initialize event log components with seeding events
    el = EventLog(E + I)

    while E + I > 0

        # Look up time for next infected event
        t_next = top(infecteds).next_event_time

        @assert t_next ≥ t "Next infected event time $t_next is earlier than current time $t"

        # Work out whether infection or temporal event occurs next
        if !isempty(susceptibles) && top(susceptibles).resistance ≤ Λ + dΛ * (t_next - t) # Infection event

            @assert dΛ > 0.0 "Infection event triggered with zero force of infection"

            # Remove susceptible from heap
            @unpack resistance, id = pop!(susceptibles)

            # Update time
            t += (resistance - Λ) / dΛ

            # Update number of exposed and susceptibles
            E += 1; S -= 1

            # Update cumulative infection pressure
            Λ = resistance

            # Create new exposed individual
            expose_host!(pop, id, t, td)

            # Add new infected to heap
            push!(infecteds, InfKey(pop.next_event_time[id], id))

            update_event_log!(el, t, id, sample_infector(actives), EK_Transmission)

        else    # Temporal event (activation, recovery, sampling)
            # Withdraw infected from heap
            @unpack next_event_time, id = pop!(infecteds)

            # Retrieve scheduled event for this infected
            scheduled_event = pop.next_event[id]

            update_event_log!(el, t_next, id, 0, scheduled_event)

            # Update cumulative infection pressure and time
            Λ += dΛ * (t_next - t)
            t = t_next

            if scheduled_event == EK_Activation   # Activation event
                # Update number of exposed and infected
                E -= 1; I += 1

                # Activate host
                activate_host!(pop, id, td, r, actives)

                # Increment slope
                dΛ = actives.total_β / N

            elseif scheduled_event == EK_SerialSampling || scheduled_event == EK_Recovery   # Sampling with removal or recovery event
                # Update number of infected and recovered
                I -= 1; R += 1
                deactivate_host!(pop, id, actives)
                dΛ = actives.total_β / N
            elseif scheduled_event == EK_FossilizedSampling
                # No change in counts, no change in slope
                fossilize_host!(pop, id, td, r)
            end

            # If events remain for this infected, push back onto heap
            pop.next_event[id] != EK_None && push!(infecteds, InfKey(pop.next_event_time[id], id))

        end
    end
    return el
end