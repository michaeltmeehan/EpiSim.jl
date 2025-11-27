

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
# function _assign_nodes!(nodes::Vector{Union{Nothing, Node}}, events::Vector{<:Event}, max_id::Int)
#     # Initialize vector of tree assignments for each host
#     host_tree = zeros(Int, max_id)

#     # Initialize counter for tree IDs
#     tree_count = 0

#     # Keep a vector of leading nodes (one for each tree)
#     leading_node = Int[]
#     for i in length(events):-1:1
#         event = events[i]
#         if event isa SerialSampling
#             # Start new tree
#             tree_count += 1
#             # Assign host to tree
#             host_tree[event.host] = tree_count
#             # Mark node as leaf
#             nodes[i] = SampledLeaf(i, event.time, event.host)
#             # Add to leading nodes
#             push!(leading_node, i)
#         elseif event isa FossilizedSampling
#             if host_tree[event.host] == 0   # Tree not yet assigned to host == host not subsequently sampled
#                 # Start new tree
#                 tree_count += 1
#                 # Assign host to tree
#                 host_tree[event.host] = tree_count
#                 # Mark node as leaf
#                 nodes[i] = SampledLeaf(i, event.time, event.host)
#                 # Add to leading nodes
#                 push!(leading_node, i)
#             else    # Host has subsequent sampled event; already assigned to tree; add fossil as unary node
#                 child = leading_node[host_tree[event.host]]
#                 nodes[i] = SampledUnary(i, event.time, event.host, child)
#                 leading_node[host_tree[event.host]] = i
#             end
#         elseif event isa Activation
#             if host_tree[event.host] != 0
#                 # Host has already been assigned to tree; add activation as unary node
#                 child = leading_node[host_tree[event.host]]
#                 nodes[i] = UnsampledUnary(i, event.time, event.host, child)
#                 # Update leading node
#                 leading_node[host_tree[event.host]] = i
#             end
#         elseif event isa Transmission
#             if host_tree[event.host] != 0 && host_tree[event.infector] != 0
#                 # Both infector and infectee assigned to tree; add transmission as binary node
#                 left_child = leading_node[host_tree[event.infector]]
#                 right_child = leading_node[host_tree[event.host]]
#                 nodes[i] = Binary(i, event.time, event.host, left_child, right_child)
#                 host_tree[event.infector] = host_tree[event.host]
#                 leading_node[host_tree[event.host]] = i
#             elseif host_tree[event.host] != 0
#                 # Only infectee assigned to tree; add transmission as unary node and assign infector to tree
#                 host_tree[event.infector] = host_tree[event.host]
#                 nodes[i] = UnsampledUnary(i, event.time, event.host, leading_node[host_tree[event.host]])
#                 leading_node[host_tree[event.host]] = i
#             end
#         elseif event isa Seeding
#             if host_tree[event.host] != 0
#                 # Host already assigned to tree; add seed as root node
#                 nodes[i] = Root(i, event.time, event.host, leading_node[host_tree[event.host]])
#             end
#         end
#     end
# end


# function _assign_nodes(events::EventLog, max_id::Int)
#     # Calculate number of events
#     n_events = length(events.time)

#     # Initialize vector of tree properties
#     all_ids = collect(1:n_events)
#     all_lefts = zeros(Int, n_events)
#     all_rights = zeros(Int, n_events)
#     all_kinds = fill(K_None, n_events)

#     # Initialize vector of tree assignments for each host
#     host_tree = zeros(Int, max_id)

#     # Initialize counter for tree IDs
#     leaf_count = 0
#     tree_ids = Int[]

#     # Keep a vector of leading nodes (one for each tree)
#     leading_node = Int[]
#     for i in n_events:-1:1
#         event_kind = events.kind[i]
#         host = events.host[i]
#         left, right = 0, 0
#         node_kind = K_None
#         if event_kind == EK_SerialSampling || (event_kind == EK_FossilizedSampling && host_tree[events.host[i]] == 0)
#             # Start new tree
#             leaf_count += 1; push!(tree_ids, leaf_count)
#             # Assign host to tree
#             host_tree[host] = leaf_count
#             host_tree[host] = @view leaf_count
#             # Mark node as leaf
#             left, right = 0, 0
#             node_kind = K_SampledLeaf
#             # Add to leading nodes
#             push!(leading_node, i)
#         elseif event_kind == EK_FossilizedSampling
#             # Host has subsequent sampled event; already assigned to tree; add fossil as unary node
#             left, right = leading_node[host_tree[host]], 0
#             node_kind = K_SampledUnary
#             leading_node[host_tree[host]] = i
#         elseif event_kind == EK_Activation && host_tree[host] != 0
#             # Host has already been assigned to tree; add activation as unsampled unary node
#             left, right = leading_node[host_tree[host]], 0
#             node_kind = K_UnsampledUnary
#             # Update leading node
#             leading_node[host_tree[host]] = i
#         elseif event_kind == EK_Transmission
#             if host_tree[host] != 0 && host_tree[events.infector[i]] != 0
#                 # Both infector and infectee assigned to tree; add transmission as binary node
#                 left, right = leading_node[host_tree[host]], leading_node[host_tree[events.infector[i]]]
#                 node_kind = K_Binary
#                 host_tree[events.infector[i]] = host_tree[events.host[i]]
#                 leading_node[host_tree[events.host[i]]] = i
#             elseif host_tree[events.host[i]] != 0
#                 # Only infectee assigned to tree; add transmission as unary node and assign infector to tree
#                 left, right = leading_node[host_tree[events.host[i]]], 0
#                 node_kind = K_UnsampledUnary
#                 leading_node[host_tree[events.host[i]]] = i
#                 host_tree[events.infector[i]] = host_tree[events.host[i]]
#             end
#         elseif event_kind == EK_Seeding && host_tree[events.host[i]] != 0
#             # Host already assigned to tree; add seed as root node
#             left, right = leading_node[host_tree[events.host[i]]], 0
#             node_kind = K_Root
#         end
#         all_lefts[i] = left
#         all_rights[i] = right
#         all_kinds[i] = node_kind
#     end
#     println("Host tree: ", host_tree)
#     return all_ids, all_lefts, all_rights, all_kinds
# end

function _assign_nodes(events::EventLog, max_id::Int)
    n_events = length(events.time)
    left_children  = zeros(Int, n_events)
    right_children = zeros(Int, n_events)
    node_kinds     = fill(NK_None, n_events)

    host_leading = zeros(Int, max_id)  # host → current leading node

    for i in n_events:-1:1
        kind = events.kind[i]
        host = events.host[i]

        if kind == EK_SerialSampling
            @assert host_leading[host] == 0  # if your model guarantees one sample per host
            node_kinds[i] = NK_SampledLeaf
            host_leading[host] = i

        elseif kind == EK_FossilizedSampling
            if host_leading[host] == 0
                node_kinds[i] = NK_SampledLeaf
                host_leading[host] = i
            else
                node_kinds[i]  = NK_SampledUnary
                left_children[i] = host_leading[host]
                host_leading[host] = i
            end

        elseif kind == EK_Activation && host_leading[host] != 0
            node_kinds[i]  = NK_UnsampledUnary
            left_children[i] = host_leading[host]
            host_leading[host] = i

        elseif kind == EK_Transmission
            infector = events.infector[i]
            hL = host_leading[host]
            iL = host_leading[infector]

            if hL != 0 && iL != 0
                # both lineages present → binary merge
                node_kinds[i]  = NK_Binary
                left_children[i]  = hL
                right_children[i] = iL
                # before the transmission (further back in time) only the infector lineage exists
                host_leading[infector] = i
                host_leading[host]     = 0

            elseif hL != 0
                # only host lineage present; propagate it back to infector
                node_kinds[i]  = NK_UnsampledUnary
                left_children[i] = hL
                host_leading[infector] = i
                host_leading[host]     = 0
            end

        elseif kind == EK_Seeding && host_leading[host] != 0
            node_kinds[i]  = NK_Root
            left_children[i] = host_leading[host]
            # break
            host_leading[host] = 0
        end
    end

    return left_children, right_children, node_kinds
end


function _label_trees(left_children::Vector{Int},
                      right_children::Vector{Int},
                      node_kinds::Vector{NodeKind})
    n = length(node_kinds)
    node_tree_map = zeros(Int, n)
    tree_id = 0

    @inbounds for i in 1:n
        # start a new component at each root not yet visited
        if node_kinds[i] == NK_Root && node_tree_map[i] == 0
            tree_id += 1
            stack = [i]
            while !isempty(stack)
                v = pop!(stack)
                node_tree_map[v] != 0 && continue
                node_tree_map[v] = tree_id

                lc = left_children[v]
                rc = right_children[v]
                lc != 0 && node_kinds[lc] != NK_None && push!(stack, lc)
                rc != 0 && node_kinds[rc] != NK_None && push!(stack, rc)
            end
        end
    end

    return node_tree_map
end


# function _compress_single_tree(times::Vector{Float64},
#                                hosts::Vector{Int},
#                                left_children::Vector{Int},
#                                right_children::Vector{Int},
#                                node_kinds::Vector{NodeKind})
#     n_nodes = length(left_children)
#     @assert n_nodes == length(right_children) == length(node_kinds) == length(times) == length(hosts)
#     node_map = zeros(Int, n_nodes)  # old index → new index
#     new_index = 0

#     for i in 1:n_nodes
#         if node_kinds[i] != NK_None
#             new_index += 1
#             node_map[i] = new_index
#         end
#     end

#     # println("Node map: ", node_map)

#     n_new = new_index
#     new_time = zeros(Float64, n_new)
#     new_host = zeros(Int, n_new)
#     new_left  = zeros(Int, n_new)
#     new_right = zeros(Int, n_new)
#     new_kinds = Vector{NodeKind}(undef, n_new)

#     for i in 1:n_nodes
#         ni = node_map[i]
#         if ni != 0
#             new_time[ni] = times[i]
#             new_host[ni] = hosts[i]
#             new_kinds[ni] = node_kinds[i]
#             lc = left_children[i]
#             rc = right_children[i]
#             new_left[ni]  = lc == 0 ? 0 : node_map[lc]
#             new_right[ni] = rc == 0 ? 0 : node_map[rc]
#         end
#     end

#     return Tree(new_time, new_host, node_map, new_left, new_right, new_kinds)
# end



# function _assign_nodes(events::EventLog, max_id::Int)
#     # Calculate number of events
#     n_events = length(events.time)

#     # Initialize vector of tree properties
#     left_children = zeros(Int, n_events)
#     right_children = zeros(Int, n_events)
#     node_kinds = fill(NK_None, n_events)

#     # Initialize vector of tree provisional assignments for each host
#     host_tree_map = zeros(Int, max_id)

#     # Initialize counters for leaves (which start as provisional trees)
#     leaf_count = 0

#     # Provisional tree IDs
#     tree_ids = Int[]

#     # Number of nodes per tree
#     # node_count = Int[]  # length(tree_ids) == length(node_counts)

#     # Provisional tree assignment for each node
#     node_tree_map = zeros(Int, n_events)

#     # Keep a vector of leading nodes (one for each tree)
#     leading_node = Int[]

#     # Iterate backwards through events (i.e., from latest to earliest)
#     for i in n_events:-1:1  # Use i to index nodes and events
#         event_kind = events.kind[i]
#         host = events.host[i]
#         host_tree = host_tree_map[host]
#         left, right = 0, 0  # Pre-allocate child indices
#         node_kind = NK_None  # Pre-allocate node kind

#         # Check for a leaf event
#         if event_kind == EK_SerialSampling || (event_kind == EK_FossilizedSampling && host_tree == 0)
#             # Start new tree
#             leaf_count += 1; push!(tree_ids, leaf_count)
#             # push!(node_count, 1)
#             # Assign host to tree
#             host_tree_map[host] = leaf_count
#             # Assign node to tree
#             node_tree_map[i] = leaf_count
#             # Mark node as leaf
#             node_kind = NK_SampledLeaf
#             # Add to leading nodes
#             push!(leading_node, i)
#         elseif event_kind == EK_FossilizedSampling
#             # Host has subsequent sampled event; already assigned to tree; add fossil as unary node
#             left = leading_node[host_tree]
#             node_kind = NK_SampledUnary
#             node_tree_map[i] = host_tree
#             leading_node[host_tree] = i
#             # node_count[host_tree] += 1
#         elseif event_kind == EK_Activation && host_tree != 0
#             # Host has already been assigned to tree; add activation as unsampled unary node
#             left = leading_node[host_tree]
#             node_kind = NK_UnsampledUnary
#             node_tree_map[i] = host_tree
#             # Update leading node
#             leading_node[host_tree] = i
#             # node_count[host_tree] += 1
#         elseif event_kind == EK_Transmission
#             infector = events.infector[i]
#             infector_tree = host_tree_map[infector]
#             if host_tree != 0 && infector_tree != 0
#                 # Both infector and infectee assigned to tree; add transmission as binary node
#                 left, right = leading_node[host_tree], leading_node[infector_tree]
#                 node_kind = NK_Binary
#                 node_tree_map[i] = host_tree
#                 # Merge trees
#                 tree_ids[infector_tree] = tree_ids[host_tree]
#                 # node_count[host_tree] += node_count[infector_tree] + 1
#                 # node_count[infector_tree] = 0
#                 host_tree_map[infector] = host_tree
#                 leading_node[host_tree] = i
#             elseif host_tree != 0
#                 # Only infectee assigned to tree; add transmission as unary node and assign infector to tree
#                 left = leading_node[host_tree]
#                 node_kind = NK_UnsampledUnary
#                 node_tree_map[i] = host_tree
#                 leading_node[host_tree] = i
#                 # node_count[host_tree] += 1
#                 host_tree_map[infector] = host_tree
#             end
#         elseif event_kind == EK_Seeding && host_tree != 0
#             # Host already assigned to tree; add seed as root node
#             left = leading_node[host_tree]
#             node_kind = NK_Root
#             node_tree_map[i] = host_tree
#             # node_count[host_tree] += 1
#         end
#         left_children[i] = left
#         right_children[i] = right
#         node_kinds[i] = node_kind
#         println("Event $i: e_kind=$(event_kind), n_kind=$(node_kind), left=$(left), right=$(right), tree=$(host_tree_map[host]), node_tree=$(node_tree_map)")
#     end

#     # Clean up node_tree to reflect tree IDs
#     for i in eachindex(node_tree_map)
#         t = node_tree_map[i]
#         if t != 0
#             node_tree_map[i] = tree_ids[t]
#         end
#     end

#     return left_children, right_children, node_kinds, node_tree_map #, node_count
# end


function _compute_node_counts(node_tree_map::Vector{Int})
    max_tree_id = maximum(node_tree_map)
    max_tree_id == 0 && return Int[]
    node_count = zeros(Int, max_tree_id)
    @inbounds for t in node_tree_map
        if t != 0
            node_count[t] += 1
        end
    end
    return node_count
end


function extract_sampled_trees(events::EventLog, max_id::Int)
    left, right, node_kinds = _assign_nodes(events, max_id)

    node_tree_map = _label_trees(left, right, node_kinds)
    node_counts = _compute_node_counts(node_tree_map)

    trees  = [Tree(n) for n in node_counts]
    nextpos = ones(Int, length(trees))

    @inbounds for i in eachindex(node_tree_map)
        t = node_tree_map[i]
        t == 0 && continue

        k = nextpos[t]
        nextpos[t] += 1

        tr = trees[t]
        tr.time[k]  = events.time[i]
        tr.host[k]  = events.host[i]
        tr.id[k]    = i
        tr.left[k]  = left[i]
        tr.right[k] = right[i]
        tr.kind[k]  = node_kinds[i]
    end

    return trees
end


# Backward pass
# function _build_trees(nodes::Vector{Union{Nothing, Node}})
#     trees = Vector{Vector{Node}}()
#     tree_count = 0
#     node_tree = zeros(Int, length(nodes))
#     for node in nodes
#         if node isa Root
#             tree_count += 1
#             node_tree[node.id] = tree_count
#             push!(trees, [node])
#             node_tree[node.child] = node_tree[node.id]
#         elseif node isa Unary
#             push!(trees[node_tree[node.id]], node)
#             node_tree[node.child] = node_tree[node.id]
#         elseif node isa Binary
#             push!(trees[node_tree[node.id]], node)
#             node_tree[node.left] = node_tree[node.id]
#             node_tree[node.right] = node_tree[node.id]
#         elseif node isa Leaf
#             push!(trees[node_tree[node.id]], node)
#         end
#     end
#     return trees
# end


# function _build_trees(all_times::Vector{Float64}, 
#                       all_hosts::Vector{Int}, 
#                       all_ids::Vector{Int},
#                       all_lefts::Vector{Int}, 
#                       all_rights::Vector{Int}, 
#                       all_kinds::Vector{NodeKind})
#     times = Vector{Vector{Float64}}()
#     hosts = Vector{Vector{Int}}()
#     ids = Vector{Vector{Int}}()
#     lefts = Vector{Vector{Int}}()
#     rights = Vector{Vector{Int}}()
#     kinds = Vector{Vector{NodeKind}}()
#     tree_count = 0
#     node_tree = zeros(Int, length(all_ids))
#     for i in eachindex(all_times)
#         kind = all_kinds[i]
#         if kind == K_None
#             continue
#         end
#         time = all_times[i]
#         host = all_hosts[i]
#         id = all_ids[i]
#         left = all_lefts[i]
#         right = all_rights[i]
#         if kind == K_Root
#             tree_count += 1
#             node_tree[i] = tree_count
#             push!(times, [time])
#             push!(hosts, [host])
#             push!(ids, [id])
#             push!(lefts, [left])
#             push!(rights, [right])
#             push!(kinds, [kind])
#             node_tree[left] = node_tree[i]
#         else
#             push!(times[node_tree[i]], time)
#             push!(hosts[node_tree[i]], host)
#             push!(ids[node_tree[i]], id)
#             push!(lefts[node_tree[i]], left)
#             push!(rights[node_tree[i]], right)
#             push!(kinds[node_tree[i]], kind)
#             if kind == K_UnsampledUnary || kind == K_SampledUnary
#                 node_tree[left] = node_tree[i]
#             elseif kind == K_Binary
#                 node_tree[left] = node_tree[i]
#                 node_tree[right] = node_tree[i]
#             end
#         end
#     end
#     trees = Vector{Tree}()
#     for i in 1:tree_count
#         push!(trees, Tree(times[i], hosts[i], ids[i], lefts[i], rights[i], kinds[i]))
#     end
#     return trees # times, hosts, ids, lefts, rights, kinds
# end


# function extract_sampled_trees(events::Vector{<:Event}, max_id::Int)

#     # Initialize vector of nodes for each event
#     nodes = Vector{Union{Nothing, Node}}(undef, length(events))
#     fill!(nodes, nothing)

#     # Assign nodes to events
#     _assign_nodes!(nodes, events, max_id)

#     return _build_trees(nodes)
# end


# function extract_sampled_trees(events::EventLog, max_id::Int)

#     all_times, all_hosts, all_ids, all_lefts, all_rights, all_kinds = _assign_nodes(events, max_id)

#     return _build_trees(all_times, all_hosts, all_ids, all_lefts, all_rights, all_kinds)
# end


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

