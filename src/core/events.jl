abstract type Event end


@inline host(event::Event) = event.host
@inline time(event::Event) = event.time


Base.isless(e1::Event, e2::Event) = e1.time < e2.time


struct Seeding <: Event
    time::Float64
    host::Int
end


struct Transmission <: Event
    time::Float64
    host::Int   # <-- infectee
    infector::Int
end


struct Sampling <: Event
    time::Float64
    host::Int
end

struct Recovery <: Event
    time::Float64
    host::Int
end


struct Activation <: Event
    time::Float64
    host::Int
end


n_sampled(events::Vector{<:Event}) = count(x -> x isa Sampling, events)
n_recovered(events::Vector{<:Event}) = count(x -> x isa Recovery, events)
n_transmissions(events::Vector{<:Event}) = count(x -> x isa Transmission, events)
n_activations(events::Vector{<:Event}) = count(x -> x isa Activation, events)
n_seeds(events::Vector{<:Event}) = count(x -> x isa Seeding, events)
n_events(events::Vector{<:Event}) = length(events)