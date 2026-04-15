"""
    TransmissionTreeView

Lightweight transmission-edge view derived from one [`EventLog`](@ref).

This is an event-log inspection layer, not a tree-native data structure. Each
entry records one `EK_Transmission` row as `(infector, infectee, time)` in event
log order. Seeding events, removals, activations, and sampling events are not
edges in this view.
"""
struct TransmissionTreeView
    infector::Vector{Int}
    infectee::Vector{Int}
    time::Vector{Float64}
end

"""
    TransmissionEdge

One transmission edge row from a [`TransmissionTreeView`](@ref).
"""
struct TransmissionEdge
    infector::Int
    infectee::Int
    time::Float64
end


Base.length(view::TransmissionTreeView) = length(view.time)
Base.isempty(view::TransmissionTreeView) = isempty(view.time)
Base.firstindex(::TransmissionTreeView) = 1
Base.lastindex(view::TransmissionTreeView) = length(view)
Base.eltype(::Type{TransmissionTreeView}) = TransmissionEdge


function Base.getindex(view::TransmissionTreeView, i::Integer)
    return TransmissionEdge(view.infector[i], view.infectee[i], view.time[i])
end


function Base.iterate(view::TransmissionTreeView, state::Int=firstindex(view))
    state > lastindex(view) && return nothing
    return (view[state], state + 1)
end


function Base.show(io::IO, view::TransmissionTreeView)
    host_count = length(union(view.infector, view.infectee))
    print(io, "TransmissionTreeView(", length(view), " edges")
    print(io, ", ", host_count, " linked hosts")
    if isempty(view)
        print(io, ", no transmission events")
    else
        print(io, ", time range ", view.time[1], " to ", view.time[end])
    end
    print(io, ")")
end


function Base.show(io::IO, edge::TransmissionEdge)
    print(io, "TransmissionEdge(", edge.infector, " -> ", edge.infectee,
          " at time ", edge.time, ")")
end


"""
    TransmissionChain

Ancestry chain for one host derived from a [`TransmissionTreeView`](@ref).

`host_id` is the requested focal host. `host_path` lists the known transmission
ancestry from the earliest host with no recorded infector in the view to
`host_id`.
`infection_time[i]` is `nothing` for the first host when that host has no known
infector in the transmission view; later entries are the transmission times that
created the corresponding host.
"""
struct TransmissionChain
    host_id::Int
    host_path::Vector{Int}
    infection_time::Vector{Union{Nothing,Float64}}
end

"""
    TransmissionChainStep

One host path step from a [`TransmissionChain`](@ref).
"""
struct TransmissionChainStep
    host_id::Int
    infection_time::Union{Nothing,Float64}
end


Base.length(chain::TransmissionChain) = length(chain.host_path)
Base.isempty(chain::TransmissionChain) = isempty(chain.host_path)
Base.firstindex(::TransmissionChain) = 1
Base.lastindex(chain::TransmissionChain) = length(chain)
Base.eltype(::Type{TransmissionChain}) = TransmissionChainStep


function Base.getindex(chain::TransmissionChain, i::Integer)
    return TransmissionChainStep(chain.host_path[i], chain.infection_time[i])
end


function Base.iterate(chain::TransmissionChain, state::Int=firstindex(chain))
    state > lastindex(chain) && return nothing
    return (chain[state], state + 1)
end


function Base.show(io::IO, chain::TransmissionChain)
    print(io, "TransmissionChain(host ", chain.host_id)
    print(io, ", ", length(chain), " hosts in ancestry path")
    if isempty(chain.host_path)
        print(io, ", host not observed in transmission view")
    elseif length(chain.host_path) == 1
        print(io, ", no known infector in transmission view")
    else
        print(io, ", source ", first(chain.host_path), " -> host ", last(chain.host_path))
    end
    print(io, ")")
end


function Base.show(io::IO, step::TransmissionChainStep)
    print(io, "TransmissionChainStep(host_id=", step.host_id,
          ", infection_time=", step.infection_time,
          ")")
end


"""
    transmission_tree(log::EventLog)

Return a [`TransmissionTreeView`](@ref) containing one edge per transmission
event in `log`.

The result preserves event-log order and stores only source host, target host,
and transmission time. It is intended for inspection of who-infected-whom
relationships while keeping the full event log as the canonical simulation
record.
"""
function transmission_tree(log::EventLog)
    n = event_count(log, EK_Transmission)
    infector = Vector{Int}(undef, n)
    infectee = Vector{Int}(undef, n)
    time = Vector{Float64}(undef, n)

    j = 0
    @inbounds for i in eachindex(log.kind)
        if log.kind[i] == EK_Transmission
            j += 1
            infector[j] = log.infector[i]
            infectee[j] = log.host[i]
            time[j] = log.time[i]
        end
    end

    return TransmissionTreeView(infector, infectee, time)
end


"""
    transmission_edges(view::TransmissionTreeView)
    transmission_edges(log::EventLog)

Return transmission edges as named tuples with fields `infector`, `infectee`,
and `time`.
"""
function transmission_edges(view::TransmissionTreeView)
    edges = Vector{NamedTuple{(:infector, :infectee, :time),Tuple{Int,Int,Float64}}}(undef, length(view))
    @inbounds for i in eachindex(view.time)
        edges[i] = (infector=view.infector[i], infectee=view.infectee[i], time=view.time[i])
    end
    return edges
end


transmission_edges(log::EventLog) = transmission_edges(transmission_tree(log))


"""
    transmission_chain(view::TransmissionTreeView, host_id)
    transmission_chain(log::EventLog, host_id)

Return the known transmission ancestry path for `host_id`.

The chain follows infectee-to-infector links backward through transmission
events, then returns the path in source-to-host order. If `host_id` has no
recorded infector, the returned path contains only `host_id`. This helper does
not infer seed status or non-transmission events; use the original
[`EventLog`](@ref) for the canonical event record.
"""
function transmission_chain(view::TransmissionTreeView, host_id::Integer)
    host = Int(host_id)
    host > 0 || throw(ArgumentError("host_id must be positive"))

    parent = Dict{Int,Tuple{Int,Float64}}()
    @inbounds for i in eachindex(view.time)
        parent[view.infectee[i]] = (view.infector[i], view.time[i])
    end

    reverse_path = Int[host]
    reverse_times = Union{Nothing,Float64}[nothing]
    seen = Set{Int}([host])
    current = host

    while haskey(parent, current)
        source, t = parent[current]
        source in seen && throw(ArgumentError("transmission edges contain a cycle involving host $source"))
        push!(reverse_path, source)
        push!(reverse_times, t)
        push!(seen, source)
        current = source
    end

    host_path = reverse(reverse_path)
    edge_times = reverse(reverse_times)
    infection_time = Vector{Union{Nothing,Float64}}(undef, length(edge_times))
    infection_time[1] = nothing
    @inbounds for i in 2:length(edge_times)
        infection_time[i] = edge_times[i - 1]
    end

    return TransmissionChain(host, host_path, infection_time)
end


transmission_chain(log::EventLog, host_id::Integer) = transmission_chain(transmission_tree(log), host_id)
