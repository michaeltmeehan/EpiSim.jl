

# function get_seeds(events::Vector{<:Event})
#     return [event for event in events if event isa Seeding]
# end

# # TODO: Forward methods to Simulation type
# # @forward get_seeds(events::Vector{<:Event}) Simulation


# function get_sampled_events(events::Vector{<:Event}, max_id::Int)
#     keep_id = falses(max_id)
#     keep_ev = falses(length(events))
#     seeds = Vector{Int}()
#     for i in length(events):-1:1
#         event = events[i]
#         if event isa Sampling
#             keep_ev[i] = true
#             keep_id[event.host] = true

#         elseif event isa Transmission
#             if keep_id[event.host]
#                 keep_ev[i] = true
#                 keep_id[event.infector] = true
#             end
#         elseif event isa Activation
#             if keep_id[event.host]
#                 keep_ev[i] = true
#             end
#         elseif event isa Seeding
#             if keep_id[event.host]
#                 keep_ev[i] = true
#                 push!(seeds, event.host)
#             end
#         end
#     end
#     reverse!(seeds)
#     return events[keep_ev]
# end


# function get_subtrees(events::Vector{<:Event}, max_id::Int)
#     trees = Dict{Int,Vector{Event}}()
#     host_tree = zeros(Int, max_id)

#     for i in 1:length(events)
#         event = events[i]
#         if event isa Seeding
#             host_tree[event.host] = event.host
#         elseif event isa Transmission
#             host_tree[event.host] = host_tree[event.infector]
#         end
#         if !haskey(trees, host_tree[event.host])
#             trees[host_tree[event.host]] = Event[]
#         end
#         push!(trees[host_tree[event.host]], event)
#     end
#     return trees
# end


# function get_sampled_subtrees(events, max_id::Int)
#     sampled_events = get_sampled_events(events, max_id)
#     return get_subtrees(sampled_events, max_id)
# end

# Forward pass
function _assign_nodes!(nodes::Vector{Union{Nothing, Node}}, events::Vector{<:Event}, max_id::Int)
    # Initialize vector of tree assignments for each host
    host_tree = zeros(Int, max_id)

    # Initialize counter for tree IDs
    tree_count = 0

    # Keep a vector of leading nodes (one for each tree)
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
end


# Backward pass
function _build_trees(nodes::Vector{Union{Nothing, Node}})
    trees = Vector{Vector{Node}}()
    tree_count = 0
    node_tree = zeros(Int, length(nodes))
    for node in nodes
        if node isa Root
            tree_count += 1
            node_tree[node.id] = tree_count
            push!(trees, [node])
            node_tree[node.child] = node_tree[node.id]
        elseif node isa Unary
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


function extract_sampled_trees(events::Vector{<:Event}, max_id::Int)

    # Initialize vector of nodes for each event
    nodes = Vector{Union{Nothing, Node}}(undef, length(events))
    fill!(nodes, nothing)

    # Assign nodes to events
    _assign_nodes!(nodes, events, max_id)

    return _build_trees(nodes)
end


# Generic tree tests
# 1. Each tree should always have one, and only one root node
# 2. Each node should have a unique ID
# 3. All nodes should be time-ordered (trees[i] == sort(trees[i]))
# 4. Each tree should be internally consistent (i.e., child nodes should have later times than parent nodes)
# 5. Each tree should be connected (i.e., no nodes are orphaned)
# 6. Children should appear after parents in the node list (i.e., nodes are listed in order of increasing time)

# Tests specific to trees constructed from epidemic events
# 1. The number of trees should ≤ the number of seeds
# 2. Events should not be duplicated within or across trees
# 3. All sampled individuals should be included in the trees
# 4. Every ancestor of a sampled node should be included in the tree
# 5. Each host should only appear in one tree
# 6. If sampling is complete then the number of trees should equal the number of seeds (set recovery to Inf)
# 7. If there are no sampling events then there should be no trees (set sampling rate to 0)
# 8. Root/Seeding correspondence: Each tree root should correspond to a seeding event for the same host at the same time
# 9. Leaf/Sampling correspondence: Each tree leaf should correspond to a sampling event for the same host at the same time
# 10. Activation/Unary correspondence: Every UnsampledUnary node whose events[id] isa Activation should correspond to an activation event for the same host at the same time
# 11. Transmission/Binary correspondence: Every Binary node should correspond to a transmission event for the same host at the same time, with left and right children corresponding to infector and infectee respectively
# 12. Recoviery events should not appear in the trees
# 13. Host field consistency: For every node, the host field should match the host of the corresponding event
# 14. Every tree has at least one sampled node
# 15. Every tree terminates with a sampled leaf node
# 16. Seeds and hosts without subsequent sampling events do not appear in the trees

