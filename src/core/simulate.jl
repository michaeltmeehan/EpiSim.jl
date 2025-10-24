
abstract type AbstractSimulation end


struct Simulation{M<:AbstractModel,S<:AbstractState,E<:AbstractEvent} <: AbstractSimulation
    model::M
    t::Vector{Float64}
    states::Vector{S}
    events::Vector{E}
    seed::Int
end


# Length/indices
Base.length(sim::Simulation) = length(sim.t)
Base.firstindex(sim::Simulation) = 1
Base.lastindex(sim::Simulation) = length(sim)
Base.IndexStyle(::Type{<:Simulation}) = IndexLinear()

# Random access: sim[i] -> SimStep(...)
Base.getindex(sim::Simulation{M,S,E}, i::Int) where {M,S<:AbstractState,E<:AbstractEvent} = (sim.t[i], sim.states[i], sim.events[i])

# Iteration: for step in sim ...
Base.iterate(sim::Simulation) = _iterate(sim, 1)
Base.iterate(sim::Simulation, i::Int) = _iterate(sim, i)

# internal helper
@inline function _iterate(sim::Simulation, i::Int)
    i > length(sim) && return nothing
    return (sim[i], i + 1)
end


# What a Simulation yields when iterated
Base.eltype(::Type{<:Simulation{<:Any,S,E}}) where {S<:AbstractState,E<:AbstractEvent} = Tuple{Float64,S,E}

# Convenience views (no copying)
eacht(sim::Simulation)      = sim.t
eachstate(sim::Simulation)  = sim.states
eachevent(sim::Simulation)  = sim.events


Base.summary(sim::Simulation{M,S,E}) where {M,S,E} = begin
    n = length(sim)
    span = n == 0 ? "empty" : "t ∈ [$(sim.t[1]), $(sim.t[end])]"
    "Simulation{$(M)} with $n steps ($span)"
end


abstract type AbstractEnsemble end


struct Ensemble{M<:AbstractModel, S<:AbstractState, E<:AbstractEvent} <: AbstractEnsemble
    sims::Vector{Simulation{M,S,E}}
end


Base.length(ens::Ensemble) = length(ens.sims)
Base.firstindex(::Ensemble) = 1
Base.lastindex(ens::Ensemble) = length(ens)
Base.IndexStyle(::Type{<:Ensemble}) = IndexLinear()
Base.eltype(::Type{<:Ensemble{M,S,E}}) where {M,S,E} = Simulation{M,S,E}
Base.getindex(ens::Ensemble, i::Int) = ens.sims[i]
Base.iterate(ens::Ensemble) = iterate(ens.sims)
Base.iterate(ens::Ensemble, st) = iterate(ens.sims, st)

# mutation passthroughs
Base.push!(ens::Ensemble, sim::Simulation{M,S,E}) where {M,S,E} = (push!(ens.sims, sim); ens)
Base.append!(ens::Ensemble, sims::Vector{Simulation{M,S,E}}) where {M,S,E} = (append!(ens.sims, sims); ens)

# convenience
eachsim(ens::Ensemble) = ens.sims
pairs(ens::Ensemble) = enumerate(ens)


Base.summary(ens::Ensemble{M,S,E}) where {M,S,E} = "Ensemble{$(M)} with $(length(ens)) simulations"


function Base.show(io::IO, ::MIME"text/plain", ens::Ensemble{M,S,E}) where {M,S,E}
    println(io, summary(ens))
    n = min(length(ens), 5)
    for i in 1:n
        sim = ens[i]
        span = isempty(sim.t) ? "empty" : "t∈[$(sim.t[1]), $(sim.t[end])]"
        println(io, " [$i] ", typeof(sim), " (", length(sim), " steps, ", span, ", seed=", string(sim.seed), ")")
    end
    length(ens) > n && println(io, " …")
end




function simulate(model::M,
                  tspan::Tuple{Float64,Float64},
                  state0::S;
                  seed::Int=rand(Int64),
                  stop::Function=((x::S) -> false)) where {M<:AbstractModel,S<:AbstractState}

    t, tf = tspan
    @assert t < tf "tspan must have t0 < tf"

    rng = MersenneTwister(seed)
    
    state = deepcopy(state0)

    EVENTS = get_events(model)

    # Store trajectory
    states = [capture(state)]
    times = Float64[t]
    events = Union{EVENTS...}[Seeding(0)]  # Initialize with seed event

    # Initialize rates vector
    λ = Vector{Float64}(undef, length(EVENTS))

    while t < tf && !stop(state)
        # Update rates
        update_rates!(λ, model, state)

        λ_tot = sum(λ)
        λ_tot == 0.0 && break

        # Sample time to next event
        t += randexp(rng) / λ_tot
        t >= tf && break

        # Sample event type
        ev_idx = wsampleindex(rng, λ)

        # Update state / apply event
        event = update_state!(rng, model, state, EVENTS[ev_idx])

        # Record event
        push!(events, event)

        # Record state and time
        push!(states, capture(state))
        push!(times, t)
    end
    return Simulation(model, times, states, events, seed)
end