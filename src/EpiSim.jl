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

include("utils/random_ops.jl")

include("models/Models.jl")
using .Models

export Model
export AbstractModel, Model
export MTBDModel, BDModel, SIRModel, SEIRModel, SuperSpreaderModel

include("simulation/simulate.jl")

export simulate, Outbreak, Ensemble

include("simulation/processing.jl")

export get_prevalence_timeseries, get_prevalence, get_prevalence_quantiles, plot_prevalence, isextinct, get_extinction_probability

end # module EpiSim
