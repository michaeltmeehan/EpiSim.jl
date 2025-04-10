mutable struct HostTree
    infector::Int
    infectees::Vector{Int}
    infection_times::Vector{Float64}
    sampling_times::Vector{Float64}
end


HostTree(infector::Int) = HostTree(infector, Vector{Int}(), Vector{Float64}(), Vector{Float64}())


# Return a triple for the i-th infectee
function Base.getindex(tree::HostTree, idx::Int)
    return (tree.infectees[idx], tree.infection_times[idx], tree.sampling_times[idx])
end

# Total number of infectees in the tree
function Base.length(tree::HostTree)
    return length(tree.infectees)
end

# Iterator
function Base.iterate(tree::HostTree, state::Int = 1)
    state > length(tree) && return nothing
    return (tree[state], state + 1)
end


function Base.push!(tree::HostTree, infectee::Int, infection_time::Float64, sampling_time::Float64)
    push!(tree.infectees, infectee)
    push!(tree.infection_times, infection_time)
    push!(tree.sampling_times, sampling_time)
    return tree
end



function Base.show(io::IO, tree::HostTree)
    n = length(tree.infectees)
    println(io, "HostTree($n infectee$(n == 1 ? "" : "s")):")
    max_show = 2  # number of infectees to show before truncating
    for idx in 1:min(n, max_show)
        infectee = tree.infectees[idx]
        inf_time = round(tree.infection_times[idx], digits=2)
        samp_time = tree.sampling_times[idx]
        sampled_str = isnan(samp_time) ? "not sampled" : "sampled at $(round(samp_time, digits=2))"
        println(io, "  Infectee $infectee: infected at $inf_time, $sampled_str")
    end
    if n > max_show
        println(io, "  ...and $(n - max_show) more")
    end
end


struct TransmissionTree
    hosts::Dict{Int, HostTree}
    root::Int
end

Base.getindex(tree::TransmissionTree, idx::Int) = tree.hosts[idx]

function Base.iterate(tree::TransmissionTree)
    return iterate(tree.hosts)
end

function Base.iterate(tree::TransmissionTree, state)
    return iterate(tree.hosts, state)
end


Base.length(tree::TransmissionTree) = length(tree.hosts)
Base.keys(tree::TransmissionTree) = keys(tree.hosts)
Base.values(tree::TransmissionTree) = values(tree.hosts)
Base.pairs(tree::TransmissionTree) = pairs(tree.hosts)




function Base.show(io::IO, tree::TransmissionTree)
    println(io, "TransmissionTree rooted at individual ", tree.root)
    n = length(tree.hosts)
    max_show = 4  # number of infectors to show before truncating
    shown = 0
    for (infector, host_tree) in sort(collect(tree.hosts); by=first)
        println(io, "\nInfector $infector:")
        show(io, host_tree)
        shown += 1
        if shown >= max_show
            remaining = n - shown
            if remaining > 0
                println(io, "\n...and $remaining more infectors not shown")
            end
            break
        end
    end
end


function extract_sampled_tree(chain::TransmissionChain)
    return extract_sampled_tree(chain.infectors, chain.infection_times, chain.sampling_times)
end


function extract_sampled_tree(infectors::Vector{Int}, 
                              infection_times::Vector{Float64}, 
                              sampling_times::Vector{Float64})
    sampled_ancestors = fill(false, length(infectors))
    hosts = Dict{Int, HostTree}()
    for infectee in eachindex(infectors)
        if !isnan(sampling_times[infectee])
            ancestor = infectee
            while ancestor > 0 && !sampled_ancestors[ancestor]
                infector = infectors[ancestor]
                tree = get!(hosts, infector, HostTree(infector))
                push!(tree, ancestor, infection_times[ancestor], sampling_times[ancestor])
                sampled_ancestors[ancestor] = true
                ancestor = infector
            end
        end
    end
    sort_tree_by_infection_time!(hosts)
    return TransmissionTree(hosts, 0)
end


function sort_tree_by_infection_time!(tree::Dict{Int, HostTree})
    for i in values(tree)
        if !issorted(i.infection_times)
            perm = sortperm(i.infection_times)
            i.infectees = i.infectees[perm]
            i.infection_times = i.infection_times[perm]
            i.sampling_times = i.sampling_times[perm]
        end
    end
end
