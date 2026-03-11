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
include("models/birthdeath.jl")
include("models/reparameterize.jl")

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

# Trees
include("tree/tree.jl")
include("tree/extract.jl")
include("tree/tree_stats.jl")


# Probability
include("birthdeath/pgf.jl")
include("birthdeath/likelihood.jl")

# Inference
include("birthdeath/utils.jl")
include("birthdeath/inference/mle_constant.jl")
include("birthdeath/inference/ensemble_fit.jl")


# =========================================
# Exports
# =========================================

export SEIRModel,
       SIRModel,
       BirthDeathModel,
       SellkeEngine,
       GillespieEngine,
       StopWhenCumulativeInfected,
       StopWhenTimeReached,
       StopWhenCumulativeSampled,
       NoStopping,
       simulate,
       simulate_ensemble,
       nsampled,
       reconstruct_trajectory,
       reconstruct_on_grid,
       infected_distribution_matrix,
       StateTrajectory,
       du_infected_distribution,
       SIRGenerator,
       SEIRGenerator,
       LineList,
       secondary_cases,
       infectious_periods,
       generation_intervals,
       empirical_R,
       empirical_R_completed,
       empirical_sampling_proportion,
       susceptible_fraction,
       normalized_offspring,
       pₙ,
       γ, α, β,
       extract_sampled_tree,
       bd_loglikelihood_constant,
       fit_bd_ensemble_mle,
       bd_pairs_plot,
       epidemic_parameter_mapping,
       Tree,
       TreeStats,
       tree_statistics,
       R0DeltaParameterization,
       RateParameterization,
       BDFixedSpec,
       fit_bd_full,
       fit_bd_pars

end