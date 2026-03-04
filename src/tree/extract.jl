
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


function validate_tree_against_log(log::EventLog, tree::Tree)

    ############################################
    # Build event indices
    ############################################

    event_index = Dict{Tuple{EventKind,HostID,Time},Bool}()
    tx_index    = Dict{Tuple{HostID,HostID,Time},Bool}()

    sampling_count = 0

    for ev in log
        event_index[(ev.kind, ev.host, ev.t)] = true

        if ev.kind == Transmission
            # key = (infectee, infector, time)
            tx_index[(ev.host, ev.src, ev.t)] = true
        end

        if ev.kind == SerialSampling || ev.kind == FossilisedSampling
            sampling_count += 1
        end
    end

    ############################################
    # 4. No sampling → no tree
    ############################################

    if sampling_count == 0 && length(tree) != 0
        error("Tree exists but no sampling events in log.")
    end

    ############################################
    # 1. Root ↔ Seeding correspondence
    ############################################

    for (_, node) in pairs(tree)
        if node.kind == Root
            if !haskey(event_index, (Seeding, node.host, node.time))
                error("Root does not correspond to Seeding event (host=$(node.host), time=$(node.time)).")
            end
        end
    end

    ############################################
    # 2. Leaf ↔ Sampling correspondence
    ############################################

    for (_, node) in pairs(tree)
        if node.kind == SampledLeaf
            ok = haskey(event_index, (SerialSampling, node.host, node.time)) ||
                 haskey(event_index, (FossilisedSampling, node.host, node.time))
            if !ok
                error("SampledLeaf does not correspond to sampling event (host=$(node.host), time=$(node.time)).")
            end
        end
    end

    ############################################
    # Transmission ↔ Binary / UnsampledUnary
    ############################################

    for (i, node) in pairs(tree)

        if node.kind == Binary

            infectee_host = tree.host[node.right]
            infector_host = node.host

            if !haskey(tx_index, (infectee_host, infector_host, node.time))
                error("Binary node at index $i does not match Transmission event.")
            end

        elseif node.kind == UnsampledUnary

            infectee_host = tree.host[node.left]
            infector_host = node.host

            if !haskey(tx_index, (infectee_host, infector_host, node.time))
                error("UnsampledUnary node at index $i does not match Transmission event.")
            end
        end
    end

    ############################################
    # SampledUnary ↔ FossilisedSampling
    ############################################

    for (_, node) in pairs(tree)
        if node.kind == SampledUnary
            if !haskey(event_index, (FossilisedSampling, node.host, node.time))
                error("SampledUnary does not correspond to FossilisedSampling event (host=$(node.host), time=$(node.time)).")
            end
        end
    end

    ############################################
    # 6. Removal events must not appear
    ############################################

    for (_, node) in pairs(tree)
        if haskey(event_index, (Removal, node.host, node.time))
            error("Removal event appears in tree (host=$(node.host), time=$(node.time)).")
        end
    end

    ############################################
    # 3. Every sampled leaf traces to a root
    ############################################

    for (i, node) in pairs(tree)
        if node.kind == SampledLeaf
            current = i
            while tree.parent[current] != 0
                current = tree.parent[current]
            end
            if tree.kind[current] != Root
                error("Sampled leaf at index $i does not trace to a Root.")
            end
        end
    end

    ############################################
    # 8. Every root has ≥1 sampled descendant
    ############################################

    for (i, node) in pairs(tree)
        if node.kind == Root

            found_sample = false
            stack = [i]

            while !isempty(stack)
                j = pop!(stack)

                if tree.kind[j] == SampledLeaf
                    found_sample = true
                    break
                end

                l = tree.left[j]
                r = tree.right[j]

                l != 0 && push!(stack, l)
                r != 0 && push!(stack, r)
            end

            if !found_sample
                error("Root at index $i has no sampled descendant.")
            end
        end
    end

    ############################################
    # 9. No extraneous nodes (every node has sampled descendant)
    ############################################

    for (i, _) in pairs(tree)

        found_sample = false
        stack = [i]

        while !isempty(stack)
            j = pop!(stack)

            if tree.kind[j] == SampledLeaf
                found_sample = true
                break
            end

            l = tree.left[j]
            r = tree.right[j]

            l != 0 && push!(stack, l)
            r != 0 && push!(stack, r)
        end

        if !found_sample
            error("Node at index $i has no sampled descendant.")
        end
    end

    ############################################
    # No duplicate event-node mapping
    ############################################

    seen = Set{Tuple{NodeKind,HostID,Time}}()

    for (_, node) in pairs(tree)
        key = (node.kind, node.host, node.time)
        if key in seen
            error("Duplicate node-event mapping detected for host=$(node.host), time=$(node.time).")
        end
        push!(seen, key)
    end

    return true
end