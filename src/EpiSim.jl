module EpiSim

using DataStructures
using Distributions
using Random
using UnPack

# Utilities
include("utils/random.jl")

# Public API
include("core/types.jl")
include("core/eventlog.jl")

# Models
include("models/abstract_model.jl")
include("models/seir.jl")
include("models/sir.jl")

# Engines
include("engines/abstract_engine.jl")
include("engines/sellke.jl")
include("engines/gillespie.jl")

# Post-processing
include("postprocess/trajectories.jl")


# =========================================
# Exports
# =========================================

export SEIRModel,
       SIRModel,
       SellkeEngine,
       GillespieEngine,
       simulate,
       reconstruct_trajectory,
       StateTrajectory

end