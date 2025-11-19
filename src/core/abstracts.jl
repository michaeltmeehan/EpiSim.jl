abstract type AbstractModel end

abstract type AbstractState end

abstract type AbstractCount <: AbstractState end

abstract type AbstractAgent <: AbstractState end


struct State
    time::Float64
    S::Int64
    E::Int64
    I::Int64
    R::Int64
end