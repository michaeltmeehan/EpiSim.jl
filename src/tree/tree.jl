abstract type Node end
abstract type Internal <: Node end
abstract type Unary <: Internal end
abstract type Leaf <: Node end


struct SampledLeaf <: Leaf
    id::Int
    time::Float64
    host::Int
end


struct UnsampledLeaf <: Leaf
    id::Int
    time::Float64
    host::Int
end


struct Binary <: Internal
    id::Int
    time::Float64
    host::Int
    left::Int
    right::Int
end


struct SampledUnary <: Unary
    id::Int
    time::Float64
    host::Int
    child::Int
end


struct UnsampledUnary <: Unary
    id::Int
    time::Float64
    host::Int
    child::Int
end


struct Root <: Node
    id::Int
    time::Float64
    host::Int
    child::Int
end


Base.isless(n1::Node, n2::Node) = n1.time < n2.time


@enum NodeKind::UInt8 begin
    NK_None = 0
    NK_Root = 1
    NK_Binary = 2
    NK_UnsampledUnary = 3
    NK_SampledUnary = 4
    NK_SampledLeaf = 5
    NK_UnsampledLeaf = 6
end


const NODEKIND_LABELS = Dict(
    NK_None           => "None",
    NK_Root           => "Root",
    NK_Binary         => "Binary",
    NK_UnsampledUnary => "UnsampledUnary",
    NK_SampledUnary   => "SampledUnary",
    NK_SampledLeaf    => "SampledLeaf",
    NK_UnsampledLeaf  => "UnsampledLeaf"
)


function Base.show(io::IO, x::NodeKind)
    print(io, NODEKIND_LABELS[x])
end

mutable struct Tree
    time::Vector{Float64}   # NaN: none
    host::Vector{Int}   # 0 : none
    id::Vector{Int}   # 0 : none
    left::Vector{Int}   # 0 : none
    right::Vector{Int}   # 0 : none
    kind::Vector{NodeKind}  # 0: None, 1: Root, 2: Binary, 3: UnsampledUnary, 4: SampledUnary, 5: SampledLeaf, 6: UnsampledLeaf
end


function Tree(n::Int)
    return Tree(Vector{Float64}(undef, n), Vector{Int}(undef, n), Vector{Int}(undef, n), Vector{Int}(undef, n), Vector{Int}(undef, n), Vector{NodeKind}(undef, n))
end


Base.length(t::Tree) = length(t.time)

