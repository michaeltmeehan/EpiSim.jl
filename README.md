# EpiSim.jl

EpiSim.jl is a Julia package for stochastic epidemic outbreak simulation and event-log generation.

The current recovered package scope is intentionally narrow:

- stochastic epidemic simulation through the active `sellke` and `gillespie` engines
- semantically documented event logs
- validation utilities for event-log consistency
- compact ensemble summaries with optional raw event-log retention
- event-time state-count recovery for a single event log
- compact event-log processing helpers for single logs

The package does not provide tree extraction, tree likelihoods, birth-death analytical likelihoods, or orchestration logic for the broader modelling ecosystem.

## Event-time state counts

Recover a compact trajectory from one `EventLog` by supplying the initial state.
The first trajectory entry is the initial state at time zero; each later entry is
the state immediately after one event-log row. Simultaneous events remain as
separate entries with repeated times.

```julia
using EpiSim

log = gillespie(100, 0, 1, 0.8, 1.0, 1.0, 0.0, 0.0)
trajectory = event_time_state_counts(log; S0=100, E0=0, I0=1)
time, infectious = trajectory_series(trajectory, :I)
```

## Event-log helpers

Interrogate one `EventLog` without changing its compact columnar storage:

```julia
transmission_times = event_times(log, EK_Transmission)
hosts = host_event_summary(log)
host_ids, transmissions = host_series(hosts, :transmissions)
```

## Ensemble summaries

Use `run_ensemble` with any RNG-explicit simulator that returns an `EventLog`.
By default, only compact per-replicate statistics are retained:

```julia
using EpiSim
using Random

summary = run_ensemble(
    rng -> gillespie(rng, 100, 0, 1, 0.8, 1.0, 1.0, 0.0, 0.0),
    1_000;
    rng=Random.MersenneTwister(42),
)

mean_final_size(summary)
mean_final_time(summary)
attack_rate(summary, 101)
ensemble_stats = ensemble_aggregate_summary(summary)
```

Set `retain_logs=true` when later derived layers need the full event logs:

```julia
summary = run_ensemble(
    rng -> sellke(rng, 100, 0, 1, 0.8, 1.0, 1.0, Inf, 0.0),
    100;
    rng=Random.MersenneTwister(42),
    retain_logs=true,
)

trajectories = ensemble_state_trajectories(summary; S0=100, E0=0, I0=1)
host_summaries = ensemble_host_event_summaries(summary)
peaks = peak_infectious(trajectories)
trajectory_stats = trajectory_aggregate_summary(trajectories)
host_stats = host_aggregate_summary(host_summaries)
```

## Installation

To install EpiSim.jl, use the following command in the Julia REPL:
```julia
using Pkg
Pkg.add(url="https://github.com/michaeltmeehan/EpiSim.jl.git")
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.
