# EpiSim.jl

EpiSim.jl is a Julia package for stochastic epidemic outbreak simulation and event-log generation.

The current recovered package scope is intentionally narrow:

- stochastic epidemic simulation through the active `sellke` and `gillespie` engines
- semantically documented event logs
- validation utilities for event-log consistency
- compact ensemble summaries with optional raw event-log retention
- event-time state-count recovery for a single event log
- compact event-log processing helpers for single logs
- lightweight transmission-edge and ancestry-chain views derived from event logs

The package does not provide tree-native algorithms, tree likelihoods,
birth-death analytical likelihoods, plotting, fitting, or orchestration logic
for the broader modelling ecosystem.

## Installation

To install EpiSim.jl, use the following command in the Julia REPL:

```julia
using Pkg
Pkg.add(url="https://github.com/michaeltmeehan/EpiSim.jl.git")
```

## Output Model

Simulation functions return an `EventLog`. It is the canonical record of what
happened in one simulated outbreak: event time, affected host, transmission
source when relevant, and event kind.

The other analysis objects are derived views over that record:

- `HostEventSummary`: per-observed-host participation counts for one log
- `StateCountTrajectory`: recovered SEIR counts at raw event times for one log
- `TransmissionTreeView`: one who-infected-whom edge per transmission event
- `TransmissionChain`: known transmission ancestry for one host
- aggregate summaries: lightweight mean/min/max summaries across replicates or
  across already-derived per-replicate views

These views are complementary. They make common questions easier to ask without
replacing the event log or changing its compact representation.

## Single-Log Workflow

Start with one simulation, then derive only the views needed for the question at
hand:

```julia
using EpiSim

log = gillespie(100, 0, 1, 0.8, 1.0, 1.0, 0.0, 0.0)
```

Inspect raw event rows when timing or event types matter:

```julia
total_events(log)
transmission_times = event_times(log, :transmission)
first_removal = first_event_time(log, :removal)
```

Use host summaries when the question is about host participation:

```julia
hosts = host_event_summary(log)
host_ids, transmissions_caused = host_series(hosts, :transmissions)
```

Use a state trajectory when the question is about compartment counts over the
event sequence:

```julia
trajectory = event_time_state_counts(log; S0=100, E0=0, I0=1)
time, infectious = trajectory_series(trajectory, :I)
```

Use transmission views when the question is who infected whom:

```julia
tree = transmission_tree(log)
edges = transmission_edges(tree)
chain = transmission_chain(tree, 10)
```

`tree` is not a tree-native object. It is a compact transmission-edge view over
the same event log, and `chain` follows those edges for one host.

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

Interrogate one `EventLog` without changing its compact columnar storage.
Event-kind helpers preserve log order:

```julia
transmission_times = event_times(log, :transmission)
activation_rows = event_indices(log, :activation)
first_removal = first_event_time(log, :removal)
```

The enum values, such as `EK_Transmission`, remain available for low-level code,
but common event-log helpers also accept readable symbols. Accepted names are
`:seeding`, `:transmission`, `:activation`, `:removal`,
`:fossilized_sampling`, `:fossilised_sampling`, and `:serial_sampling`; plural
forms and matching strings are also accepted.

## Transmission views

Extract a compact who-infected-whom view when transmission relationships are
the focus:

```julia
tree = transmission_tree(log)
edges = transmission_edges(tree)
chain = transmission_chain(tree, 10)
```

`tree` is a `TransmissionTreeView`: one edge per transmission event, preserving
event-log order with `infector`, `infectee`, and `time` vectors. `edges` is a
plain named-tuple list for tabular inspection. `chain` is a `TransmissionChain`
for one host, listing the known ancestry path from the earliest known source to
that host.

These views deliberately omit activation, removal, and sampling events. Use the
original `EventLog` for the canonical event record, `HostEventSummary` for
per-host participation counts, and `StateCountTrajectory` for compartment
counts over event time.

## Ensemble summaries

Use `run_ensemble` when the question is about many independent replicates. The
simulator must return an `EventLog`. By default, only compact per-replicate
statistics are retained:

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

`summary` is an `EnsembleSummary`: it stores one value per replicate for final
size, final time, event counts, and sampling counts. `ensemble_stats` summarizes
those per-replicate vectors with mean, minimum, and maximum values.

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
peak_times = peak_infectious_time(trajectories)
trajectory_stats = trajectory_aggregate_summary(trajectories)
host_stats = host_aggregate_summary(host_summaries)
```

`trajectory_stats` summarizes scalar values derived independently from each
trajectory; it does not interpolate onto a shared time axis. Peak time is the
first time the peak is reached within each trajectory, then summarized across
trajectories. `host_stats` computes per-host quantities within each replicate
first, then summarizes those replicate-level values across replicates; host IDs
are not pooled across the ensemble.


## License

This project is licensed under the MIT License - see the LICENSE file for details.
