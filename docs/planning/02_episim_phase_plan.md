## `EpiSim.jl` next-phase plan

### Role

`EpiSim.jl` is the epidemic simulation core responsible for stochastic outbreak simulation and event-log generation.

### Current position

`EpiSim.jl` has been recovered, validated, benchmarked, and cleaned of legacy or out-of-scope code.

### Immediate priorities

#### 1. Ensemble infrastructure

Goals:
- support repeated simulations cleanly
- provide compact ensemble outputs
- support configurable retention of raw trajectories / logs / summaries
- preserve reproducibility and efficient RNG handling

#### 2. State-time-series recovery

Goals:
- reconstruct epidemic state trajectories over time from simulation outputs
- support state trajectories on event times and optionally user-defined grids
- support recovery of counts such as susceptible / exposed / infectious / removed where model-appropriate

#### 3. Host summaries

Candidate summaries:
- realized secondary case counts
- realized infection source counts
- infectious periods
- exposure-to-infectious delays where applicable
- sampling histories
- removal times

#### 4. Processing utilities

- filter event logs by host
- filter event logs by event type
- derive host histories
- aggregate outbreak-level statistics
- convert raw outputs into analysis-friendly tables or lightweight summary structs

#### 5. Visualization support

- epidemic curves
- state trajectories
- host-summary histograms
- outbreak summary distributions across ensembles

### Explicit deferrals

- direct dependency on `TreeSim.jl`
- direct dependency on `BDUtils.jl`
- embedding tree-native logic inside `EpiSim.jl`
- full calibration / optimization framework

### Architectural guardrails

- event logs remain a primary internal truth
- efficient simulation kernels remain central
- derived state recovery should be layered on top
- host summaries should be computed from compact internals
- ensemble support should be designed for scale

### Suggested implementation sequence

| Order | Work packet |
|---|---|
| 1 | Ensemble API design |
| 2 | Compact ensemble summary layer |
| 3 | State-time-series recovery utilities |
| 4 | Host summary functions |
| 5 | Event-log processing helpers |
| 6 | Visualization helpers |
