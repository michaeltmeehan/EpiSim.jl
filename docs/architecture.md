# EpiSim.jl – Architecture Overview (Sellke-First Design)

This document defines the intended clean architecture of EpiSim.jl after refactoring toward a Sellke-first simulation framework with systematic benchmarking support.

Primary goals:

1. Efficient stochastic simulation of epidemic models using a Sellke construction.
2. Full agent-level transmission histories.
3. Clean derivation of aggregate trajectories from event logs.
4. Compatibility with aggregate CTMC formulations for differentiated uniformization (DU).
5. Straightforward future extension to multi-type models.

The Sellke engine is the primary simulation engine. Aggregate and Gillespie engines exist for validation and benchmarking.

---

# 1. Canonical Event Representation

The canonical simulation output is an `EventLog` stored in **struct-of-arrays** form for performance and post-processing efficiency.

Each event must contain:

- `t::Float64` — event time
- `kind::EventKind` — event type
- `host::Int` — host ID experiencing the event
- `src::Int` — infector ID (meaningful only for transmission; otherwise 0)

No additional metadata field is required. Distinct event kinds are used instead of a generic sampling event with subtype flags.

## Event Kinds

```julia
@enum EventKind::UInt8 begin
    Seeding
    Transmission
    Activation
    Removal
    SerialSampling
    FossilisedSampling
end
```

### Semantics

| Event              | Meaning                                                 |
| ------------------ | ------------------------------------------------------- |
| Seeding            | Initial infected host (no parent)                       |
| Transmission       | Infection caused by another host                        |
| Activation         | Latent → infectious transition (if latent stage exists) |
| Removal            | Infectious → removed transition                         |
| SerialSampling     | Sampling while infectious                               |
| FossilisedSampling | Sampling associated with removal (tip event)            |

This structure ensures clarity, eliminates subtype flags, and simplifies post-processing.

---

# 2. Event Log Invariants

The EventLog must satisfy the following invariants:

- Event times are non-decreasing.
- Each infected host has exactly one Seeding or Transmission event.
- Every Transmission event has a valid `src`.
- A host cannot activate before infection.
- A host cannot be removed before activation (if latent stage exists).
- Each host has at most one Removal event.
- FossilisedSampling must occur at or immediately before Removal.
- Derived aggregate counts must never become negative.

The event log must be fully sufficient to reconstruct:

- Transmission trees
- Individual infection timelines
- Aggregate S/E/I/R trajectories
- Sampling histories

No hidden state outside the event log may be required.

---

# 3. Sellke Engine (Primary Agentic Engine)

The Sellke engine is the core simulation engine.

It implements:

- Pre-drawn susceptibility resistances
- Heap-based scheduling of infection pressure
- Explicit tracking of active infectors
- Explicit scheduling of activation, removal, and sampling events

The engine assumes an SEIR-like structure but does not hard-code specific distributions.

## Model Inputs

The model must provide:

- Initial counts: `S0`, `E0`, `I0`
- Infection pressure formulation (mass-action)
- Incubation time draw (optional; zero for SIR)
- Infectious duration draw
- Sampling time draw(s) (serial and/or fossilised)
- Transmission scaling parameters (e.g., β)

The engine calls model-provided draw functions and does not assume exponential distributions.

## Output

The Sellke engine returns:

- `EventLog`

All aggregate trajectories must be derived from the event log in post-processing.

This ensures:

- Single source of truth
- No duplicated state tracking
- Clean benchmarking interface
- Efficient memory usage

---

# 4. Aggregate CTMC Interface (For DU Benchmarking)

To enable differentiated uniformization benchmarking, the epidemic process must also be expressible as a finite-state CTMC.

This layer defines:

- State encoding: `(S, E, I, R)`
- Transition types:
  - Infection
  - Activation
  - Removal
- Rate functions for each transition

This interface enables:

- Construction of the infinitesimal generator
- Exact transient distribution computation via DU
- Parameter sensitivity analysis

The aggregate CTMC must match the induced aggregate process of the Sellke model under equivalent parameterisation.

---

# 5. Gillespie Reference Engine (Minimal)

A minimal Gillespie implementation is retained for:

- Aggregate SIR/SEIR simulations
- Cross-validation against Sellke-derived aggregate trajectories
- Validation against DU results

This engine:

- Operates only on counts
- Uses the same CTMC specification as the DU layer
- Does not record agent-level histories

It serves strictly as a validation tool.

---

# 6. Benchmarking Strategy

Benchmarking has two layers.

## 6.1 Deterministic Invariant Checks

Applied to every simulation:

- Monotone event times
- Valid host state transitions
- Correct infection source accounting
- Valid aggregate counts derived from the event log
- No negative compartment sizes

These checks must live in the test suite.

## 6.2 Statistical Validation

For selected parameter sets:

- Compare Sellke-derived aggregate distributions to:
  - Gillespie simulations
  - Differentiated uniformization results
- Compare:
  - Final size distribution
  - Extinction probability
  - Prevalence distributions at fixed time points
- Validate convergence under increasing Monte Carlo replicates

---

# 7. Multi-Type Future Extension

To support future multi-type models:

- Host IDs remain contiguous integers.
- Event log format remains type-agnostic.
- Internal aggregate representation must be easily generalisable from scalars to vectors per type.
- Avoid hard-coding scalar `S, E, I, R` fields in core abstractions.

No multi-type implementation is required yet, but architectural choices must not prevent it.

---

# 8. Design Principles

- Single canonical event representation
- Sellke engine as the primary simulation engine
- Aggregate trajectories derived only from the event log
- Clear separation between:
  - Agent simulation
  - Aggregate CTMC specification
  - Benchmarking utilities
- Minimal duplication of logic
- Refactor before extending

This document defines the target architecture for EpiSim.jl moving forward.
