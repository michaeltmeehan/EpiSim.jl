module EpiSim

# using CairoMakie # Experimental
# using ColorSchemes
import Base.Threads
using Crayons
using DataFrames
using DataFramesMeta
using LinearAlgebra
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


include("core/abstracts.jl")

include("utils/random.jl")

export popr!, wsample, wsampleindex, wsampleindex_cols


include("core/events.jl")

export Seeding, Transmission, Sampling, Recovery, Activation


include("core/simulate.jl")

export simulate, eacht, eachsim, eachstate, eachevent, Simulation, Ensemble


include("models/sir.jl")

export SIR, SIRAgent, SIRCount


include("models/bd.jl")

export BD, BDAgent, BDCount


include("models/mtbd.jl")

export MTBD, MTBDAgent, MTBDCount, SSBD, SSBDAgent, SSBDCount


end # module EpiSim
