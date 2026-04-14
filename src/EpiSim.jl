module EpiSim

# using CairoMakie # Experimental
# using ColorSchemes
using DataStructures: BinaryMinHeap, push!, pop!
using Distributions
using Random
using UnPack


# include("core/abstracts.jl")

# export State

include("utils/random.jl")

export popr!, wsample, wsampleindex, wsampleindex_cols


include("core/events.jl")

export EventLog, EventKind, EK_None, EK_Seeding, EK_Transmission, EK_FossilizedSampling, EK_SerialSampling, EK_Recovery, EK_Activation

# include("core/simulation.jl")

# export Simulation, Ensemble


# include("core/simulate.jl")

# export simulate, eacht, eachsim, eachstate, eachevent, Simulation, Ensemble


# include("models/sir.jl")

# export SIR, SIRAgent, SIRCount


# include("models/bd.jl")

# export BD, BDAgent, BDCount


# include("models/mtbd.jl")

# export MTBD, MTBDAgent, MTBDCount, SSBD, SSBDAgent, SSBDCount


include("core/sellke.jl")

export sellke, TraitDists


include("core/gillespie.jl")

export gillespie


include("tree/tree.jl")

export NodeKind, NK_None, NK_Root, NK_Binary, NK_SampledLeaf, NK_UnsampledLeaf, NK_SampledUnary, NK_UnsampledUnary
export Tree

include("utils/processing.jl")

export extract_sampled_trees


include("birthdeath/likelihood.jl")

export likelihood, p₀, γ, γ₀, γ₁, β

end # module EpiSim
