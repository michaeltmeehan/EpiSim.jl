module Models

using DataFrames
using Lazy
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

export AbstractEpiModel, AbstractEpiState, AbstractEpiStateSlice
export BirthDeathModel, SIRModel, SEIRModel, MultiTypeBirthDeathModel, SuperSpreaderModel
export simulate_events
export capture, isagentic



export SIRState
export EpiState, get_event_types, get_default_stop_condition, initialize_event_log, update_event_rates!, update_state!

end