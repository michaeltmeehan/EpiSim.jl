module EpiSim

# using CairoMakie # Experimental
# using ColorSchemes
using DataFrames
using DataFramesMeta
using Lazy
using Parameters
using Plots
using RecipesBase
using RecipesPipeline
using RecursiveArrayTools
using StatsBase
using UnPack

include("parameters.jl")
include("outbreak.jl")
include("models/bdutils.jl")
include("models/crbd.jl")
include("models/mtbd.jl")
include("process.jl")
include("plot.jl")

export EpiParameters, BirthDeathParameters, CRBDParameters, MTBDParameters, SIRParameters, Outbreak
export initialize_linelist, initialize_active, check_parameters
export simulate_outbreak, summarize, type_dist, offspring_dist, n_sampled, n_deceased
export plot_outbreak, create_transmission_graph



end # module EpiSim
