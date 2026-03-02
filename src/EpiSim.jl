module EpiSim

using DataStructures
using Distributions
using Random
using UnPack

# Utilities
# include("utils/random.jl")

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
include("core/ensemble.jl")

# Post-processing
include("postprocess/trajectories.jl")
include("postprocess/summaries.jl")
include("postprocess/ensemble_stats.jl")


# =========================================
# Exports
# =========================================

export SEIRModel,
       SIRModel,
       SellkeEngine,
       GillespieEngine,
       simulate,
       simulate_ensemble,
       reconstruct_trajectory,
       reconstruct_on_grid,
       infected_distribution_matrix,
       StateTrajectory

end