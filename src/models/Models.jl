module Models

using Random
using ..EpiEvent
import LinearAlgebra: ⋅
import StatsBase: sample, wsample
using ..EpiSim: pop_random!

include("AbstractEpiModel.jl")
include("AbstractEpiState.jl")
include("SIRModel.jl")
include("SEIRModel.jl")
include("BirthDeathModel.jl")
include("MultiTypeBirthDeathModel.jl")

export BirthDeathModel, SIRModel, SEIRModel, MultiTypeBirthDeathModel
export simulate_outbreak

# include("Gillespie.jl")

export AbstractEpiModel, AbstractEpiState

export SIRState
export EpiState, get_event_types, get_default_stop_condition, initialize_event_log, update_event_rates!, update_state!


end