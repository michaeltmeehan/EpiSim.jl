module EpiSim

# using CairoMakie # Experimental
# using ColorSchemes
import Base.Threads
using Crayons
using DataFrames
using DataFramesMeta
using DataStructures: BinaryMinHeap, push!, pop!, top
using Distributions
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

export State

include("utils/random.jl")

export popr!, wsample, wsampleindex, wsampleindex_cols


include("core/events.jl")

export Seeding, Transmission, Sampling, FossilizedSampling, SerialSampling, Recovery, Activation, Event

include("core/simulation.jl")

export Simulation, Ensemble


# include("core/simulate.jl")

# export simulate, eacht, eachsim, eachstate, eachevent, Simulation, Ensemble


# include("models/sir.jl")

# export SIR, SIRAgent, SIRCount


# include("models/bd.jl")

# export BD, BDAgent, BDCount


# include("models/mtbd.jl")

# export MTBD, MTBDAgent, MTBDCount, SSBD, SSBDAgent, SSBDCount


include("core/sellke.jl")

export sellke, Susceptible, Infected, Agent, TraitDists, Traits
export isdone, make_infected, next_event, next_time


include("core/gillespie.jl")

export gillespie


include("utils/processing.jl")

export get_seeds, get_subtree, get_sampled_events, get_subtrees, get_sampled_subtrees


end # module EpiSim
