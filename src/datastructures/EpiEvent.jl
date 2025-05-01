module EpiEvent

import Base: show

using ColorSchemes
using Crayons
using StatsBase
using RecipesBase
using Plots

export Seed, Transmission, Activation, Sampling, Recovery, AbstractEpiEvent
export incidence
export transmission!, sampling!, recovery!, activation!, get_sampled_events

abstract type AbstractEpiEvent end

time(event::AbstractEpiEvent) = getfield(event, :time)
host(event::AbstractEpiEvent) = getfield(event, :host)

const SEED_COLOR         = Crayon(foreground = :green, bold = true)
const TRANSMISSION_COLOR = Crayon(foreground = :yellow)
const SAMPLING_COLOR     = Crayon(foreground = :blue)
const RECOVERY_COLOR     = Crayon(foreground = :red)
const ACTIVATION_COLOR   = Crayon(foreground = :cyan)


struct Seed <: AbstractEpiEvent
    host::Int
    time::Float64
end


struct Transmission <: AbstractEpiEvent
    infector::Int
    infectee::Int
    time::Float64
end


function transmission!(events::Vector{T}, infector::Int, infectee::Int, time::Float64) where T <: AbstractEpiEvent
    push!(events, Transmission(infector, infectee, time))
end


struct Sampling <: AbstractEpiEvent
    host::Int
    time::Float64
end


function sampling!(events::Vector{T}, sampled::Int, time::Float64) where T <: AbstractEpiEvent
    push!(events, Sampling(sampled, time))
end


struct Recovery <: AbstractEpiEvent
    host::Int
    time::Float64
end


function recovery!(events::Vector{T}, recovered::Int, time::Float64) where T <: AbstractEpiEvent
    push!(events, Recovery(recovered, time))
end


struct Activation <: AbstractEpiEvent
    host::Int
    time::Float64
end


function activation!(events::Vector{AbstractEpiEvent}, host::Int, time::Float64)
    push!(events, Activation(host, time))
end


function show(io::IO, e::Seed)
    msg = "Seed(host $(e.host), time = $(round(e.time, digits=2)))"
    print(io, SEED_COLOR(msg))
end

function show(io::IO, e::Transmission)
    msg = "Transmission($(e.infector) → $(e.infectee), time = $(round(e.time, digits=2)))"
    print(io, TRANSMISSION_COLOR(msg))
end

function show(io::IO, e::Sampling)
    msg = "Sampling(host $(e.host), time = $(round(e.time, digits=2)))"
    print(io, SAMPLING_COLOR(msg))
end

function show(io::IO, e::Recovery)
    msg = "Recovery(host $(e.host), time = $(round(e.time, digits=2)))"
    print(io, RECOVERY_COLOR(msg))
end

function show(io::IO, e::Activation)
    msg = "Activation(host $(e.host), time = $(round(e.time, digits=2)))"
    print(io, ACTIVATION_COLOR(msg))
end


"""
    incidence(events; bin_width=1.0, t_start=nothing, t_stop=nothing)

Compute a histogram of transmission times from a list of EpiEvents.
You can optionally specify:
  - `bin_width`: the width of each time bin (default = 1.0)
  - `t_start`: leftmost edge of the first bin (default = floor(min(t_times)))
  - `t_stop`: rightmost edge of the last bin (default = ceil(max(t_times)))
Returns a `Histogram` object.
"""
function incidence(
    events::Vector{<:AbstractEpiEvent};
    bin_width::Float64 = 1.0,
    t_start::Union{Nothing, Float64} = nothing,
    t_stop::Union{Nothing, Float64} = nothing
)
    t_times = [time(e) for e in events if e isa Transmission]
    isempty(t_times) && return fit(Histogram, Float64[], collect(0.0:bin_width:1.0))  # empty dummy histogram

    # Determine range
    t_min = isnothing(t_start) ? floor(minimum(t_times)) : t_start
    t_max = isnothing(t_stop)  ? ceil(maximum(t_times))  : t_stop

    edges = collect(t_min:bin_width:t_max)
    return fit(Histogram, t_times, edges)
end


function Base.show(io::IO, ::MIME"text/plain", events::Vector{<:AbstractEpiEvent})
    nevents = length(events)
    if nevents == 0
        println(io, "Event log is empty.")
        return
    end

    # Extract stats
    sorted = sort(events, by = time)
    times = time.(sorted)
    t_min, t_max = minimum(times), maximum(times)

    nseed = count(e -> e isa Seed, events)
    ntransmission = count(e -> e isa Transmission, events)
    nactivated = count(e -> e isa Activation, events)
    nsampled  = count(e -> e isa Sampling, events)
    nrecovered = count(e -> e isa Recovery, events)

    println(io, Crayon(bold=true)("Epidemic Event Log Summary"))
    println(io, "  Total events:    $nevents")
    nseed > 0 && println(io, SEED_COLOR("  Seeds:           $nseed"))
    ntransmission > 0 && println(io, TRANSMISSION_COLOR("  Transmissions:   $ntransmission"))
    nactivated > 0 && println(io, ACTIVATION_COLOR("  Activations:     $nactivated"))
    nsampled > 0 && println(io, SAMPLING_COLOR("  Samplings:       $nsampled"))
    nrecovered > 0 && println(io, RECOVERY_COLOR("  Recoveries:      $nrecovered"))
    print(io, "  Time range:      $(round(t_min, digits=2)) to $(round(t_max, digits=2))")
    # println()

    # Draw incidence
    # ninfected > 0 && draw_log_incidence(io, events; bin_width=1.0, t_start=t_min, t_stop=t_max, max_width=60)

    # Show tail
    # println(io, "\n  Last 5 events:")
    # for e in last(sorted, min(5, nevents))
    #     println(io, "    ", e)
    # end
end


function draw_log_incidence(io::IO, events::Vector{<:AbstractEpiEvent}; bin_width=1.0, t_start=nothing, t_stop=nothing, max_width=60)
    h = incidence(events; bin_width=bin_width)
    edges = h.edges[1]
    counts = h.weights

    n_bins = length(counts)
    symbols = collect(" ▁▂▃▄▅▆▇█")

    max_log = log10(maximum(counts) + 1)
    log_scaled = [log10(c + 1) / max_log for c in counts]
    heights = [round(Int, s * 8) for s in log_scaled]

    # Downsample to target width if needed
    if n_bins > max_width
        step = ceil(Int, n_bins / max_width)
        heights = [maximum(heights[i:min(i+step-1, end)]) for i in 1:step:n_bins]
    end

    n_bars = length(heights)
    line = join([symbols[clamp(h+1, 1, 9)] for h in heights])
    println(io, "  ", line)

    # Align left and right labels under the histogram
    label_left  = string(round(edges[1], digits=1))
    label_right = string(round(edges[end], digits=1))

    space_between = n_bars - length(label_left) - length(label_right)
    space_between = max(space_between, 1)  # avoid negative spacing

    println(io, "  ", label_left, repeat(" ", space_between), label_right)
end


@recipe function f(events::Vector{<:AbstractEpiEvent})

    transmission_times = [time(event) for event in events if (event isa Transmission || event isa Seed)]

    n_hosts = length(transmission_times)

    legend := false

    xlabel := "Time"
    ylabel := "Host ID"
    title := "Epidemic Event Log"

    linewidth_val = clamp(20.0 / sqrt(n_hosts), 1.0, 10.0)
    markersize_val = clamp(7.0 / sqrt(n_hosts), 2.0, 4.0)
    
    alpha := 0.7


    final_time = time(events[end])
    bar_drawn = BitVector(fill(false, n_hosts))
    colors = Vector{RGB{Float64}}(undef, n_hosts)

    for event in events
        if event isa Seed
            bar_drawn[event.host] = false
            colors[event.host] = ColorSchemes.tab20.colors[mod1(event.host, 20)]
        elseif event isa Transmission
            bar_drawn[event.infectee] = false
            colors[event.infectee] = ColorSchemes.tab20.colors[mod1(event.infectee, 20)]
            x = [event.time, event.time]
            y = [event.infector, event.infectee]
            @series begin
                linewidth := linewidth_val / 5.
                color := colors[event.infector]
                x, y
            end
        elseif event isa Sampling || event isa Recovery
            bar_drawn[event.host] = true
            x = [transmission_times[event.host], event.time]
            y = [event.host, event.host]
            @series begin
                linewidth := linewidth_val
                color := colors[event.host]
                x, y
            end

            @series begin
                seriestype := :scatter
                markercolor := colors[event.host]
                markerstrokecolor := :white
                markerstrokewidth := markersize_val / 4.
                markersize := markersize_val
                label := ""
                [transmission_times[event.host]], [y]
            end
        end
    end

    # Draw bars for remaining hosts
    for host in eachindex(transmission_times)
        if !bar_drawn[host]
            x = [transmission_times[host], final_time]
            y = [host, host]
            @series begin
                linewidth := linewidth_val
                color := colors[host]
                x, y
            end

            @series begin
                seriestype := :scatter
                markercolor := colors[host]
                markerstrokecolor := :white
                markerstrokewidth := markersize_val / 4.
                markersize := markersize_val
                label := ""
                [transmission_times[host]], [y]
            end
        end
    end
end



function get_sampled_events(events::Vector{AbstractEpiEvent})
    sampled_events = Vector{AbstractEpiEvent}()
    hosts = Vector{Int}()
    for event in reverse(events)
        if event isa Sampling
            push!(sampled_events, event)
            push!(hosts, event.host)
        elseif event isa Transmission && event.infectee in hosts
            push!(sampled_events, event)
            push!(hosts, event.infector)
        elseif event isa Seed && event.host in hosts
            push!(sampled_events, event)
            push!(hosts, event.host)
        end
    end
    return reverse(sampled_events)
end



end