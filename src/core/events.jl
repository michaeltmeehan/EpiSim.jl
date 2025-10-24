abstract type AbstractEvent end


@inline host(event::AbstractEvent) = event.host


struct Seeding <: AbstractEvent
    host::Int
end


struct Transmission <: AbstractEvent
    host::Int   # <-- infectee
    infector::Int
end


struct Sampling <: AbstractEvent
    host::Int
end


struct Recovery <: AbstractEvent
    host::Int
end



struct Activation <: AbstractEvent
    host::Int
end


n_sampled(events::Vector{<:AbstractEvent}) = count(x -> x isa Sampling, events)
n_recovered(events::Vector{<:AbstractEvent}) = count(x -> x isa Recovery, events)
n_transmissions(events::Vector{<:AbstractEvent}) = count(x -> x isa Transmission, events)
n_activations(events::Vector{<:AbstractEvent}) = count(x -> x isa Activation, events)
n_seeds(events::Vector{<:AbstractEvent}) = count(x -> x isa Seeding, events)
n_events(events::Vector{<:AbstractEvent}) = length(events)