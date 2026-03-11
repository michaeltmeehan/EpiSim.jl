############################################################
# Tree statistics
############################################################

struct TreeStats
    ntips::Int
    height::Float64
    total_branch_length::Float64

    mean_internal_branch::Float64
    mean_terminal_branch::Float64
    internal_terminal_ratio::Float64

    mean_node_depth::Float64
    var_node_depth::Float64

    sackin::Int
    colless::Int
    cherries::Int

    branch_length_cv::Float64

    largest_clade_fraction::Float64
end


############################################################
# Main statistics routine
############################################################
# TODO: Revisit largest_clade_fraction definition (currently largest_clade / ntips) - always 1.
function tree_statistics(tree::Tree)

    n = length(tree)

    subtree_size = zeros(Int, n)
    node_depth   = zeros(Int, n)

    ################################################
    # accumulators
    ################################################

    total_branch_length = 0.0

    internal_sum = 0.0
    terminal_sum = 0.0
    internal_count = 0
    terminal_count = 0

    branch_sum = 0.0
    branch_sqsum = 0.0
    branch_count = 0

    ntips = 0

    root_time = Inf
    max_time  = -Inf

    colless = 0
    cherries = 0
    largest_clade = 0

    ################################################
    # PASS 1  (tips → root)
    ################################################

    for (i,node) in pairs(tree)

        t = node.time
        max_time = max(max_time, t)

        if node.parent == 0
            root_time = min(root_time, t)
        end

        ################################################
        # initialize leaf subtree sizes
        ################################################

        if node.kind == SampledLeaf
            subtree_size[i] = 1
            ntips += 1
        end

        ################################################
        # propagate subtree sizes upward
        ################################################

        l = node.left
        r = node.right

        if l != 0
            subtree_size[i] += subtree_size[l]
        end

        if r != 0
            subtree_size[i] += subtree_size[r]
        end

        ################################################
        # balance statistics (binary nodes only)
        ################################################

        if node.kind == Binary

            sl = subtree_size[l]
            sr = subtree_size[r]

            colless += abs(sl - sr)

            if sl == 1 && sr == 1
                cherries += 1
            end

            largest_clade = max(largest_clade, subtree_size[i])
        end

        ################################################
        # branch-length statistics
        ################################################

        if node.parent == 0
            continue
        end

        # do not start branches at unary nodes
        if node.kind == UnsampledUnary || node.kind == SampledUnary
            continue
        end

        p = node.parent
        bl = node.time - tree.time[p]

        # collapse unary chains upward
        while p != 0 &&
              (tree.kind[p] == UnsampledUnary ||
               tree.kind[p] == SampledUnary)

            gp = tree.parent[p]
            bl += tree.time[p] - tree.time[gp]
            p = gp
        end

        total_branch_length += bl

        branch_sum   += bl
        branch_sqsum += bl*bl
        branch_count += 1

        if node.kind == SampledLeaf
            terminal_sum += bl
            terminal_count += 1
        elseif node.kind == Binary
            internal_sum += bl
            internal_count += 1
        end

    end

    ################################################
    # PASS 2  (root → tips)
    ################################################

    sackin = 0

    for (i,node) in pairs(reverse(tree))

        l = node.left
        r = node.right

        if l != 0
            node_depth[l] = node_depth[i] + 1
        end

        if r != 0
            node_depth[r] = node_depth[i] + 1
        end

        if node.kind == SampledLeaf
            sackin += node_depth[i]
        end
    end

    ################################################
    # derived statistics
    ################################################

    height = max_time - root_time

    mean_internal = internal_sum / internal_count
    mean_terminal = terminal_sum / terminal_count

    internal_terminal_ratio = mean_internal / mean_terminal

    branch_mean = branch_sum / branch_count
    branch_var  = branch_sqsum / branch_count - branch_mean^2
    branch_length_cv = sqrt(branch_var) / branch_mean

    ################################################
    # node-time depth statistics
    ################################################

    depth_sum = 0.0
    depth_sqsum = 0.0
    depth_count = 0

    for node in tree
        if node.kind == Binary
            d = max_time - node.time
            depth_sum += d
            depth_sqsum += d*d
            depth_count += 1
        end
    end

    mean_node_depth = depth_sum / depth_count
    var_node_depth  = depth_sqsum / depth_count - mean_node_depth^2

    largest_clade_fraction = largest_clade / ntips

    ################################################

    TreeStats(
        ntips,
        height,
        total_branch_length,
        mean_internal,
        mean_terminal,
        internal_terminal_ratio,
        mean_node_depth,
        var_node_depth,
        sackin,
        colless,
        cherries,
        branch_length_cv,
        largest_clade_fraction
    )

end