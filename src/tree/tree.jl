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