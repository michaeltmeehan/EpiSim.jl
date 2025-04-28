module Models

using Random
using ..EpiEvent
import LinearAlgebra: ⋅
import StatsBase: sample, wsample
using ..EpiSim: pop_random!

VALIDATE_STATE = false

include("AbstractEpiModel.jl")
include("AbstractEpiState.jl")
include("SIRModel.jl")
include("SEIRModel.jl")
include("BirthDeathModel.jl")
include("MultiTypeBirthDeathModel.jl")
include("SuperSpreaderModel.jl")

export BirthDeathModel, SIRModel, SEIRModel, MultiTypeBirthDeathModel, SuperSpreaderModel
export simulate_outbreak


export AbstractEpiModel, AbstractEpiState

export SIRState
export EpiState, get_event_types, get_default_stop_condition, initialize_event_log, update_event_rates!, update_state!


end