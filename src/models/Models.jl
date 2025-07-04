module Models

using DataFrames
using Lazy
using Random
import LinearAlgebra: ⋅, diagm
import StatsBase: sample, wsample, quantile
using ..EpiSim: pop_random!

VALIDATE_STATE = false

include("AbstractEvent.jl")

export AbstractEvent, AbstractEpiEvent, Seed, Transmission, Activation, Sampling, Recovery
export time, host

include("AbstractParameters.jl")
include("AbstractState.jl")
include("AbstractModel.jl")

export AbstractModel, Model, isagentic, capture

include("SIRModel.jl")
include("SEIRModel.jl")
include("BirthDeathModel.jl")
include("MultiTypeBirthDeathModel.jl")
include("SuperSpreaderModel.jl")

export MTBDModel, BDModel, SIRModel, SEIRModel, SuperSpreaderModel
export MTBDParameters, BDParameters, SIRParameters, SEIRParameters, SuperSpreaderParameters
export AgenticMTBDState, AgenticBDState, AgenticSIRState, AgenticSEIRState, AgenticSuperSpreaderState
export AggregateMTBDState, AggregateBDState, AggregateSIRState, AggregateSEIRState, AggregateSuperSpreaderState

export get_default_stop_condition, initialize_event_log, update_event_rates!, update_state!
export n_sampled, n_recovered, n_transmissions, n_activations, n_seeds, n_events

end