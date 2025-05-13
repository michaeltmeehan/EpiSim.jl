module EpiSim

# using CairoMakie # Experimental
# using ColorSchemes
import Base.Threads
using Crayons
using DataFrames
using DataFramesMeta
using Lazy
using Parameters
using Plots
# using Printf
using Random
using RecipesBase
using RecipesPipeline
using RecursiveArrayTools
using StatsBase
using UnPack


include("DataStructures/EpiEvent.jl")
using .EpiEvent

export incidence, get_sampled_events

include("utils/random_ops.jl")

include("models/Models.jl")
using .Models

export Models
export BirthDeathModel, SIRModel, SEIRModel, MultiTypeBirthDeathModel, SuperSpreaderModel
export simulate_events

include("simulation/simulate.jl")

export simulate, Outbreak, Ensemble

end # module EpiSim
