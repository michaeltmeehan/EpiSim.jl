

struct Simulation
    states::Vector{State}
    events::Vector{Event}
    seed::Int
end


# Length/indices
Base.length(sim::Simulation) = length(sim.states)
Base.firstindex(sim::Simulation) = 1
Base.lastindex(sim::Simulation) = length(sim)
Base.IndexStyle(::Type{<:Simulation}) = IndexLinear()
Base.getindex(sim::Simulation, i::Int) = (sim.states[i], sim.events[i])

# Iteration: for step in sim ...
Base.iterate(sim::Simulation) = _iterate(sim, 1)
Base.iterate(sim::Simulation, i::Int) = _iterate(sim, i)

# internal helper
@inline function _iterate(sim::Simulation, i::Int)
    i > length(sim) && return nothing
    return (sim[i], i + 1)
end


# What a Simulation yields when iterated
Base.eltype(::Type{<:Simulation}) = Tuple{State,Event}

# Convenience views (no copying)
eachstate(sim::Simulation)  = sim.states
eachevent(sim::Simulation)  = sim.events

Base.summary(sim::Simulation) = begin
    n = length(sim)
    span = n == 0 ? "empty" : "t ∈ [$(sim.states[1].time), $(sim.states[end].time)]"
    "Simulation with $n steps ($span)"
end



struct Ensemble
    sims::Vector{Simulation}
end


Base.length(ens::Ensemble) = length(ens.sims)
Base.firstindex(::Ensemble) = 1
Base.lastindex(ens::Ensemble) = length(ens)
Base.IndexStyle(::Type{<:Ensemble}) = IndexLinear()
Base.eltype(::Type{<:Ensemble}) = Simulation
Base.getindex(ens::Ensemble, i::Int) = ens.sims[i]
Base.iterate(ens::Ensemble) = iterate(ens.sims)
Base.iterate(ens::Ensemble, st) = iterate(ens.sims, st)

# mutation passthroughs
Base.push!(ens::Ensemble, sim::Simulation) = (push!(ens.sims, sim); ens)
Base.append!(ens::Ensemble, sims::Vector{Simulation}) = (append!(ens.sims, sims); ens)

# convenience
eachsim(ens::Ensemble) = ens.sims
pairs(ens::Ensemble) = enumerate(ens)


Base.summary(ens::Ensemble) = "Ensemble with $(length(ens)) simulations"


function Base.show(io::IO, ::MIME"text/plain", ens::Ensemble)
    println(io, summary(ens))
    n = min(length(ens), 5)
    for i in 1:n
        sim = ens[i]
        span = isempty(sim.states) ? "empty" : "t∈[$(sim.states[1].time), $(sim.states[end].time)]"
        println(io, " [$i] ", typeof(sim), " (", length(sim), " steps, ", span, ", seed=", string(sim.seed), ")")
    end
    length(ens) > n && println(io, " …")
end