module EpiSim

using DataStructures: BinaryMinHeap, push!, pop!
using Distributions
using Random
using UnPack

include("utils/random.jl")

export popr!, wsample, wsampleindex, wsampleindex_cols

# Core event-log and simulation API.
include("core/events.jl")

export EventLog, EventKind, EK_None, EK_Seeding, EK_Transmission, EK_FossilizedSampling, EK_SerialSampling, EK_Removal, EK_Activation, validate_event_log

include("core/sellke.jl")

export sellke, TraitDists

include("core/gillespie.jl")

export gillespie

# Derived analysis and summary API.
include("analysis/ensemble.jl")

export EnsembleSummary, run_ensemble, final_size, event_count, final_time, mean_final_size, mean_final_time, attack_rate

include("analysis/state_trajectory.jl")

export StateCountTrajectory, event_time_state_counts

include("analysis/event_log_helpers.jl")

export event_indices, event_times, first_event_time, last_event_time, has_event_kind, total_events, observed_hosts, distinct_host_count, HostEventSummary, host_event_summary

include("analysis/ensemble_derived.jl")

export ensemble_state_trajectories, ensemble_host_event_summaries

include("analysis/visualization_support.jl")

export trajectory_series, trajectory_final_sizes, trajectory_final_times, peak_infectious, peak_infectious_time, host_series

include("analysis/aggregate_summaries.jl")

export ScalarSummary, EnsembleAggregateSummary, TrajectoryAggregateSummary, HostAggregateSummary, ensemble_aggregate_summary, trajectory_aggregate_summary, host_aggregate_summary

end # module EpiSim
