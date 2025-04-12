module Models

# using ..Transmission
using ..EpiEvent
import LinearAlgebra: ⋅
import StatsBase: sample, wsample
using ..EpiSim: pop_random!

include("AbstractEpiModel.jl")
# include("SIRModel.jl")
# include("SEIRModel.jl")
include("BirthDeathModel.jl")
# include("MultiTypeBirthDeathModel.jl")

export BirthDeathModel, simulate_outbreak
# export SIRModel, SEIRModel, BirthDeathModel, MultiTypeBirthDeathModel, simulate_chain

end