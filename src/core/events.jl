abstract type AbstractEvent end


@inline host(event::AbstractEvent) = event.host
@inline time(event::AbstractEvent) = event.time


Base.isless(e1::AbstractEvent, e2::AbstractEvent) = e1.time < e2.time


struct Seeding <: AbstractEvent
    time::Float64
    host::Int
end


struct Transmission <: AbstractEvent
    time::Float64
    host::Int   # <-- infectee
    infector::Int
end


struct Sampling <: AbstractEvent
    time::Float64
    host::Int
end


struct Recovery <: AbstractEvent
    time::Float64
    host::Int
end



struct Activation <: AbstractEvent
    time::Float64
    host::Int
end


n_sampled(events::Vector{<:AbstractEvent}) = count(x -> x isa Sampling, events)
n_recovered(events::Vector{<:AbstractEvent}) = count(x -> x isa Recovery, events)
n_transmissions(events::Vector{<:AbstractEvent}) = count(x -> x isa Transmission, events)
n_activations(events::Vector{<:AbstractEvent}) = count(x -> x isa Activation, events)
n_seeds(events::Vector{<:AbstractEvent}) = count(x -> x isa Seeding, events)
n_events(events::Vector{<:AbstractEvent}) = length(events)