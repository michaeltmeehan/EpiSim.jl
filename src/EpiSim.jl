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
export MTBDParameters, BDParameters, SIRParameters, SEIRParameters, SuperSpreaderParameters
export AgenticMTBDState, AgenticBDState, AgenticSIRState, AgenticSEIRState, AgenticSuperSpreaderState
export AggregateMTBDState, AggregateBDState, AggregateSIRState, AggregateSEIRState, AggregateSuperSpreaderState

include("simulation/simulate.jl")

export simulate, Simulation, Ensemble
export eachstate, eachevent, eachsim

include("simulation/processing.jl")

export event_counts, get_state

export plot_prevalence, isextinct, get_extinction_probability

end # module EpiSim
