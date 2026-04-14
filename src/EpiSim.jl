module EpiSim

using DataStructures: BinaryMinHeap, push!, pop!
using Distributions
using Random
using UnPack

include("utils/random.jl")

export popr!, wsample, wsampleindex, wsampleindex_cols

include("core/events.jl")

export EventLog, EventKind, EK_None, EK_Seeding, EK_Transmission, EK_FossilizedSampling, EK_SerialSampling, EK_Recovery, EK_Activation, validate_event_log

include("core/sellke.jl")

export sellke, TraitDists

include("core/gillespie.jl")

export gillespie

end # module EpiSim
