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