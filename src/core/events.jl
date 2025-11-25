abstract type Event end
abstract type Sampling <: Event end


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

@inline hosts(event::Transmission) = [event.infector, event.host]

struct Recovery <: Event
    time::Float64
    host::Int
end


struct Activation <: Event
    time::Float64
    host::Int
end


struct FossilizedSampling <: Sampling
    time::Float64
    host::Int
end


struct SerialSampling <: Sampling
    time::Float64
    host::Int
end


n_sampled(events::Vector{<:Event}) = count(x -> x isa Sampling, events)
n_recovered(events::Vector{<:Event}) = count(x -> x isa Recovery, events)
n_transmissions(events::Vector{<:Event}) = count(x -> x isa Transmission, events)
n_activations(events::Vector{<:Event}) = count(x -> x isa Activation, events)
n_seeds(events::Vector{<:Event}) = count(x -> x isa Seeding, events)
n_events(events::Vector{<:Event}) = length(events)


@enum EventKind::UInt8 begin
    EK_None = 0
    EK_Seeding = 1
    EK_Transmission = 2
    EK_Recovery = 3
    EK_Activation = 4
    EK_FossilizedSampling = 5
    EK_SerialSampling = 6
end


struct EventLog
    time::Vector{Float64}
    host::Vector{Int}
    infector::Vector{Int}  # only for Transmission events; 0 otherwise
    kind::Vector{EventKind} # 0: None, 1: Seeding, 2: Transmission, 3: Recovery, 4: Activation, 5: FossilizedSampling, 6: SerialSampling
end


Base.length(el::EventLog) = length(el.time)

# Write iterators for time, host, infector, kind
eachtime(el::EventLog) = (el.time[i] for i in 1:length(el))
eachhost(el::EventLog) = (el.host[i] for i in 1:length(el))
eachinfector(el::EventLog) = (el.infector[i] for i in 1:length(el))
eachkind(el::EventLog) = (el.kind[i] for i in 1:length(el))