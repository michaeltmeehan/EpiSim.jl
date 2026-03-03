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
include("core/stopping.jl")

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
include("postprocess/linelist.jl")

# Differential uniformization
include("benchmark/du/du.jl")
using .DU: du_infected_distribution, SIRGenerator, SEIRGenerator

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
       StateTrajectory,
       du_infected_distribution,
       SIRGenerator,
       SEIRGenerator,
       LineList,
       secondary_cases,
       infectious_periods

end