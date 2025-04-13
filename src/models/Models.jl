module Models

using Random
using ..EpiEvent
import LinearAlgebra: ⋅
import StatsBase: sample, wsample
using ..EpiSim: pop_random!

include("AbstractEpiModel.jl")
include("SIRModel.jl")
include("SEIRModel.jl")
include("BirthDeathModel.jl")
include("MultiTypeBirthDeathModel.jl")

export BirthDeathModel, SIRModel, SEIRModel, MultiTypeBirthDeathModel
export simulate_outbreak

end