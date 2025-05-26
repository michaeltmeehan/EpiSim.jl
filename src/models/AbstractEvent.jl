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


n_sampled(event_log::Vector{<:AbstractEvent}) = count(x -> x isa Sampling, event_log)
n_recovered(event_log::Vector{<:AbstractEvent}) = count(x -> x isa Recovery, event_log)
n_transmissions(event_log::Vector{<:AbstractEvent}) = count(x -> x isa Transmission, event_log)
n_activations(event_log::Vector{<:AbstractEvent}) = count(x -> x isa Activation, event_log)
n_seeds(event_log::Vector{<:AbstractEvent}) = count(x -> x isa Seed, event_log)
n_events(event_log::Vector{<:AbstractEvent}) = length(event_log)