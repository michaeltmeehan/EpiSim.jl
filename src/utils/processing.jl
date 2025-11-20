

function get_seeds(events::Vector{<:Event})
    return [event for event in events if event isa Seeding]
end

# TODO: Forward methods to Simulation type
# @forward get_seeds(events::Vector{<:Event}) Simulation


function get_sampled_events(events::Vector{<:Event}, max_id::Int)
    keep_id = falses(max_id)
    keep_ev = falses(length(events))
    seeds = Vector{Int}()
    for i in length(events):-1:1
        event = events[i]
        if event isa Sampling
            keep_ev[i] = true
            keep_id[event.host] = true

        elseif event isa Transmission
            if keep_id[event.host]
                keep_ev[i] = true
                keep_id[event.infector] = true
            end
        elseif event isa Activation
            if keep_id[event.host]
                keep_ev[i] = true
            end
        elseif event isa Seeding
            if keep_id[event.host]
                keep_ev[i] = true
                push!(seeds, event.host)
            end
        end
    end
    reverse!(seeds)
    return events[keep_ev]
end


function get_subtrees(events::Vector{<:Event}, max_id::Int)
    trees = Dict{Int,Vector{Event}}()
    host_tree = zeros(Int, max_id)

    for i in 1:length(events)
        event = events[i]
        if event isa Seeding
            host_tree[event.host] = event.host
        elseif event isa Transmission
            host_tree[event.host] = host_tree[event.infector]
        end
        if !haskey(trees, host_tree[event.host])
            trees[host_tree[event.host]] = Event[]
        end
        push!(trees[host_tree[event.host]], event)
    end
    return trees
end


function get_sampled_subtrees(events, max_id::Int)
    sampled_events = get_sampled_events(events, max_id)
    return get_subtrees(sampled_events, max_id)
end

