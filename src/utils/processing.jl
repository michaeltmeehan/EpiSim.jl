

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


function extract_sampled_trees(events::Vector{<:Event}, max_id::Int)
    # Initialize vector of tree assignments for each host
    host_tree = zeros(Int, max_id)

    # Initialize vector of nodes for each event
    nodes = Vector{Union{Nothing, Node}}(undef, length(events))
    fill!(nodes, nothing)

    # Initialize vector of tree assignments for each node
    node_tree = zeros(Int, length(events))

    # Initialize counter for tree IDs
    tree_count = 0

    # Keep a vector of leading nodes
    leading_node = Int[]
    for i in length(events):-1:1
        event = events[i]
        if event isa SerialSampling
            # Start new tree
            tree_count += 1
            # Assign host to tree
            host_tree[event.host] = tree_count
            # Mark node as leaf
            nodes[i] = SampledLeaf(i, event.time, event.host)
            # Add to leading nodes
            push!(leading_node, i)
        elseif event isa FossilizedSampling
            if host_tree[event.host] == 0   # Tree not yet assigned to host == host not subsequently sampled
                # Start new tree
                tree_count += 1
                # Assign host to tree
                host_tree[event.host] = tree_count
                # Mark node as leaf
                nodes[i] = SampledLeaf(i, event.time, event.host)
                # Add to leading nodes
                push!(leading_node, i)
            else    # Host has subsequent sampled event; already assigned to tree; add fossil as unary node
                child = leading_node[host_tree[event.host]]
                nodes[i] = SampledUnary(i, event.time, event.host, child)
                leading_node[host_tree[event.host]] = i
            end
        elseif event isa Activation
            if host_tree[event.host] != 0
                # Host has already been assigned to tree; add activation as unary node
                child = leading_node[host_tree[event.host]]
                nodes[i] = UnsampledUnary(i, event.time, event.host, child)
                # Update leading node
                leading_node[host_tree[event.host]] = i
            end
        elseif event isa Transmission
            if host_tree[event.host] != 0 && host_tree[event.infector] != 0
                # Both infector and infectee assigned to tree; add transmission as binary node
                left_child = leading_node[host_tree[event.infector]]
                right_child = leading_node[host_tree[event.host]]
                nodes[i] = Binary(i, event.time, event.host, left_child, right_child)
                host_tree[event.infector] = host_tree[event.host]
                leading_node[host_tree[event.host]] = i
            elseif host_tree[event.host] != 0
                # Only infectee assigned to tree; add transmission as unary node and assign infector to tree
                host_tree[event.infector] = host_tree[event.host]
                nodes[i] = UnsampledUnary(i, event.time, event.host, leading_node[host_tree[event.host]])
                leading_node[host_tree[event.host]] = i
            end
        elseif event isa Seeding
            if host_tree[event.host] != 0
                # Host already assigned to tree; add seed as root node
                nodes[i] = Root(i, event.time, event.host, leading_node[host_tree[event.host]])
            end
        end
    end
    # return nodes

    # Collect trees
    trees = Vector{Vector{Node}}()
    tree_count = 0
    node_tree = zeros(Int, length(nodes))
    for node in nodes
        if node isa Root
            tree_count += 1
            node_tree[node.id] = tree_count
            push!(trees, [node])
            node_tree[node.child] = node_tree[node.id]
        elseif node isa UnsampledUnary
            push!(trees[node_tree[node.id]], node)
            node_tree[node.child] = node_tree[node.id]
        elseif node isa Binary
            push!(trees[node_tree[node.id]], node)
            node_tree[node.left] = node_tree[node.id]
            node_tree[node.right] = node_tree[node.id]
        elseif node isa Leaf
            push!(trees[node_tree[node.id]], node)
        end
    end
    return trees
end