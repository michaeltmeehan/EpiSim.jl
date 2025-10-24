abstract type AbstractEvent end

time(event::AbstractEvent) = getfield(event, :time)

abstract type AbstractEpiEvent <: AbstractEvent end

host(event::AbstractEpiEvent) = getfield(event, :host)


struct Seed <: AbstractEpiEvent
    host::Union{Int, Nothing}
    time::Float64
end


struct Transmission <: AbstractEpiEvent
    infector::Union{Int, Nothing}
    infectee::Union{Int, Nothing}
    time::Float64
end


struct Sampling <: AbstractEpiEvent
    host::Union{Int, Nothing}
    time::Float64
end


struct Recovery <: AbstractEpiEvent
    host::Union{Int, Nothing}
    time::Float64
end


struct Activation <: AbstractEpiEvent
    host::Union{Int, Nothing}
    time::Float64
end


n_sampled(events::Vector{<:AbstractEvent}) = count(x -> x isa Sampling, events)
n_recovered(events::Vector{<:AbstractEvent}) = count(x -> x isa Recovery, events)
n_transmissions(events::Vector{<:AbstractEvent}) = count(x -> x isa Transmission, events)
n_activations(events::Vector{<:AbstractEvent}) = count(x -> x isa Activation, events)
n_seeds(events::Vector{<:AbstractEvent}) = count(x -> x isa Seed, events)
n_events(events::Vector{<:AbstractEvent}) = length(events)