
function extract_sampled_tree(log::EventLog)

    # Determine maximum host id for active lineage tracking
    max_host = maximum(log.host)

    tree = Tree()
    active = zeros(Int, max_host)  # active[h] = node index or 0

    # Node creation primitive
    function new_node!(t::Time, h::HostID, l::Int, r::Int, k::NodeKind)
        push!(tree.time, t)
        push!(tree.host, h)
        push!(tree.left, l)
        push!(tree.right, r)
        push!(tree.parent, 0)
        push!(tree.kind, k)
        return length(tree.time)
    end

    for ev in reverse(log)

        ########################################
        # 1. Serial Sampling → SampledLeaf
        ########################################
        if ev.kind == SerialSampling

            # Should not already be active
            active[ev.host] != 0 && error("SerialSampling encountered but host already active.")

            i = new_node!(ev.t, ev.host, 0, 0, SampledLeaf)
            active[ev.host] = i

        ########################################
        # 2. Fossilised Sampling → SampledUnary
        ########################################
        elseif ev.kind == FossilisedSampling

            child = active[ev.host]

            # If lineage not active yet, treat like leaf
            if child == 0
                i = new_node!(ev.t, ev.host, 0, 0, SampledLeaf)
                active[ev.host] = i
            else
                i = new_node!(ev.t, ev.host, child, 0, SampledUnary)
                tree.parent[child] = i
                active[ev.host] = i
            end

        ########################################
        # 3. Transmission
        ########################################
        elseif ev.kind == Transmission

            child_host = ev.host
            parent_host = ev.src

            child_node = active[child_host]

            # If no sampled descendants → irrelevant
            child_node == 0 && continue

            parent_node = active[parent_host]

            if parent_node != 0
                # TRUE COALESCENCE → Binary
                i = new_node!(ev.t, parent_host, parent_node, child_node, Binary)

                tree.parent[parent_node] = i
                tree.parent[child_node]  = i

                active[parent_host] = i
                active[child_host]  = 0

            else
                # Only child active → Unary transfer upward
                i = new_node!(ev.t, parent_host, child_node, 0, UnsampledUnary)

                tree.parent[child_node] = i

                active[parent_host] = i
                active[child_host]  = 0
            end

        ########################################
        # 4. Seeding → Root
        ########################################
        elseif ev.kind == Seeding

            child_node = active[ev.host]

            child_node == 0 && continue  # No sampled descendants — ignore


            # Create explicit root node
            i = new_node!(ev.t, ev.host, child_node, 0, Root)

            tree.parent[child_node] = i

            active[ev.host] = 0
        end
    end

    any(active .!= 0) &&  error("Active lineages remain after processing log. Missing seeding events?")

    return tree
end