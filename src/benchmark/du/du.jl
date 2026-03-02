module DU

using LinearAlgebra
using SparseArrays

# ------------------------------------------------------------
# Abstract generator interface
# ------------------------------------------------------------

abstract type AbstractGenerator end

# Must be implemented by concrete generators
function Q_mul!(out::Vector{Float64},
                v::Vector{Float64},
                gen::AbstractGenerator)
    throw(MethodError(Q_mul!, (out, v, gen)))
end

function uniformization_rate(gen::AbstractGenerator)::Float64
    throw(MethodError(uniformization_rate, (gen,)))
end

# ------------------------------------------------------------
# Include components
# ------------------------------------------------------------

include("uniformization.jl")
include("sir_generator.jl")
include("marginal.jl")
include("api.jl")


export du_infected_distribution, SIRGenerator

end # module