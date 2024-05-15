module EpiSim

using DataFrames
using DataFramesMeta
using Lazy
using Parameters
using RecursiveArrayTools
using StatsBase
using UnPack

include("parameters.jl")
include("outbreak.jl")
include("models/bdutils.jl")
include("models/crbd.jl")
include("models/mtbd.jl")
include("process.jl")


export EpiParameters, BirthDeathParameters, CRBDParameters, MTBDParameters, Outbreak
export initialize_linelist, initialize_active, check_parameters
export simulate, summarize, type_dist, offspring_dist, n_sampled



end # module EpiSim
