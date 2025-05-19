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
export calc_R0, calc_infectious_period, calc_sampling_fraction, calc_extinction_probability
export BirthDeathModel, SIRModel, SEIRModel, MultiTypeBirthDeathModel, SuperSpreaderModel
export simulate_events

include("simulation/simulate.jl")

export simulate, Outbreak, Ensemble

include("simulation/processing.jl")

export get_prevalence_timeseries, get_prevalence, get_prevalence_quantiles, plot_prevalence, isextinct, get_extinction_probability

end # module EpiSim
