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

export EventLog, EventKind, EK_None, EK_Seeding, EK_Transmission, EK_FossilizedSampling, EK_SerialSampling, EK_Recovery, EK_Activation, validate_event_log

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

include("utils/processing.jl")

# include("birthdeath/likelihood.jl")

end # module EpiSim
