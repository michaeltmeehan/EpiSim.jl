@enum NodeKind::UInt8 begin
    Root = 0
    Binary = 1
    UnsampledUnary = 2
    SampledUnary = 3
    SampledLeaf = 4
end


mutable struct Tree
    time::Vector{Time}
    host::Vector{HostID}
    left::Vector{Int}
    right::Vector{Int}
    parent::Vector{Int}
    kind::Vector{NodeKind}
end


Tree() = Tree(Time[], HostID[], Int[], Int[], Int[], NodeKind[])


struct Node
    time::Time
    host::HostID
    left::Int
    right::Int
    parent::Int
    kind::NodeKind
end


Base.getindex(tree::Tree, i::Int) =
    Node(tree.time[i],
         tree.host[i],
         tree.left[i],
         tree.right[i],
         tree.parent[i],
         tree.kind[i])


Base.length(tree::Tree) = length(tree.time)

function Base.iterate(tree::Tree, i::Int=1)
    i > length(tree) && return nothing
    return (tree[i], i + 1)
end

Base.size(tree::Tree) = (length(tree),)


Base.iterate(r::Base.Iterators.Reverse{Tree},
             i::Int=length(r.itr)) =
    i < 1 ? nothing :
    (r.itr[i], i - 1)

Base.reverse(tree::Tree) = Base.Iterators.Reverse(tree)


Base.pairs(tree::Tree) = zip(1:length(tree), tree)


Base.pairs(r::Base.Iterators.Reverse{Tree}) =
    zip(length(r.itr):-1:1, r)
    

function validate_tree(tree::Tree)

    n = length(tree.time)

    ########################################
    # 1. Length consistency
    ########################################
    if !(length(tree.host)   == n &&
         length(tree.left)   == n &&
         length(tree.right)  == n &&
         length(tree.parent) == n &&
         length(tree.kind)   == n)
        error("Tree arrays have inconsistent lengths.")
    end

    ########################################
    # 2. Index validity
    ########################################
    for i in 1:n

        l = tree.left[i]
        r = tree.right[i]
        p = tree.parent[i]

        if l != 0 && (l < 1 || l > n || l == i)
            error("Invalid left index at node $i")
        end

        if r != 0 && (r < 1 || r > n || r == i)
            error("Invalid right index at node $i")
        end

        if p != 0 && (p < 1 || p > n || p == i)
            error("Invalid parent index at node $i")
        end
    end

    ########################################
    # 3. Parent–child consistency
    ########################################
    for i in 1:n
        l = tree.left[i]
        r = tree.right[i]

        if l != 0 && tree.parent[l] != i
            error("Parent-child mismatch (left) at node $i")
        end

        if r != 0 && tree.parent[r] != i
            error("Parent-child mismatch (right) at node $i")
        end
    end

    ########################################
    # 4. Degree constraints by NodeKind
    ########################################
    for i in 1:n
        l = tree.left[i]
        r = tree.right[i]
        p = tree.parent[i]
        k = tree.kind[i]

        if k == SampledLeaf
            if l != 0 || r != 0 || p == 0
                error("Invalid SampledLeaf structure at node $i")
            end

        elseif k == SampledUnary || k == UnsampledUnary
            if l == 0 || r != 0 || p == 0
                error("Invalid Unary structure at node $i")
            end

        elseif k == Binary
            if l == 0 || r == 0 || p == 0
                error("Invalid Binary structure at node $i")
            end

        elseif k == Root
            if l == 0 || r != 0 || p != 0
                error("Invalid Root structure at node $i")
            end

        else
            error("Unknown NodeKind at node $i")
        end
    end

    ########################################
    # 5. Time ordering
    ########################################
    for i in 1:n
        l = tree.left[i]
        r = tree.right[i]

        if l != 0 && !(tree.time[i] < tree.time[l])
            error("Time ordering violated (left child) at node $i")
        end

        if r != 0 && !(tree.time[i] < tree.time[r])
            error("Time ordering violated (right child) at node $i")
        end
    end

    ########################################
    # 6. Acyclicity check
    ########################################
    for i in 1:n
        visited = Set{Int}()
        current = i

        while current != 0
            if current in visited
                error("Cycle detected starting at node $i")
            end
            push!(visited, current)
            current = tree.parent[current]
        end
    end

    return true
end