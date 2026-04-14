# Validation Plan

## Purpose

This document defines the validation work required for `EpiSim.jl` to be considered a trustworthy simulation package for stochastic epidemic outbreak simulation and event-log generation.

`EpiSim.jl` is intended to remain independent of `TreeSim.jl` and `BDUtils.jl`, with its core responsibility being simulation and event-log generation rather than downstream analytical interpretation.

---

## Current scope

`EpiSim.jl` currently supports:

- stochastic epidemic outbreak simulation
- event-log generation
- the recovered and stabilised active simulation core
- a cleaned public API focused on simulator and event-log functionality

---

## Out of scope

The following are explicitly out of scope for the current phase:

- tree-package functionality
- birth-death analytical likelihoods
- tight integration with `TreeSim.jl`
- tight integration with `BDUtils.jl`
- orchestration logic for the full ecosystem

---

## Core validation questions

Validation work must establish that:

1. event logs are scientifically interpretable and internally consistent
2. different simulation engines agree in controlled settings where they should
3. seeded simulations are reproducible where expected
4. boundary and edge-case regimes do not silently corrupt outputs
5. the package’s current simulation scope is explicit and trustworthy

---

## Required validation areas

### 1. Event-log semantics

Validation must begin by documenting:

- the meaning of each event kind
- the meaning of host and source fields
- event ordering expectations
- what invariants downstream consumers may assume
- what is not guaranteed by the log format

### 2. Event-log invariants

Tests must cover:

- monotone time ordering where required
- valid host references
- valid source references when applicable
- event-specific source constraints
- internal consistency between events and state transitions
- consistency of cumulative counts derived from logs

### 3. Cross-engine validation

Where multiple engines support comparable scenarios, tests must compare:

- extinction probability
- prevalence/incidence summaries
- event-count summaries
- distributional behaviour at selected times where feasible
- qualitative agreement in limiting/simple regimes

### 4. Reproducibility validation

Tests must verify:

- deterministic behaviour under fixed seeds where expected
- correct RNG isolation
- reproducibility of replicate workflows where documented
- absence of unintended global RNG dependence

### 5. Boundary and stress validation

Tests should probe:

- small populations
- high and low rate regimes
- early extinction
- long simulations
- rare-event paths
- stopping-rule edge cases

---

## Canonical benchmark suite

A benchmark suite should include:

- small simple epidemic scenarios with expected qualitative behaviour
- shared-model cases for engine comparison
- fixed-seed regression cases
- edge cases expected to fail or terminate early
- examples with hand-checkable event-log properties

Each benchmark should document:

- model
- parameters
- initial conditions
- RNG setup
- expected outcome or invariant

---

## Downstream validation relevance

Although `EpiSim.jl` remains independent of `TreeSim.jl` and `BDUtils.jl`, validation should still specify what downstream users may safely assume from an event log.

This should be documented as a semantic contract, not as a dependency.

---

## Exit criteria for Phase 2

`EpiSim.jl` should not be considered fully validated for this phase until:

- event-log semantics are documented
- event-log invariants are tested thoroughly
- cross-engine validation exists for simple comparable settings
- fixed-seed reproducibility is tested
- stress tests exist for major boundary regimes
- the minimal stable simulation API is clearly defined

---

## Evidence of successful validation

Evidence should include:

- event-log invariant tests
- cross-engine comparison tests
- seeded reproducibility tests
- benchmark scenarios with documented expectations
- code-level documentation aligning event semantics with implementation

---

## Phase-2 first-pass hardening added

The first hardening pass establishes a minimal event-log semantic contract in
code and tests:

- completed logs have equal-length columns
- event times are finite, non-negative, and monotonically non-decreasing
- host ids are positive and optionally bounded by known population size
- `infector` is a transmission source only for `EK_Transmission`
- non-transmission events use `infector == 0` as the no-source sentinel
- `EK_None` is treated as an internal sentinel, not a completed log event
- seeded and transmitted hosts are not duplicated as newly infected hosts
- non-seeding/non-transmission events must refer to a previously seeded or
  transmitted host

Both active engines validate their generated logs before returning. The test
suite now includes canonical valid and malformed logs, fixed-seed
reproducibility checks, generated-log validation, and a simple finite-population
boundary case with zero transmission.

This remains intentionally lightweight. Future scientific validation can add
ensemble summaries for small SI/SIS/SIR-style finite-state models and compare
those summaries against exact or high-precision reference probabilities from
external tools such as `DifferentiatedUniformization.jl`. That comparator should
remain a validation asset rather than an architectural dependency of `EpiSim.jl`.

---

## Phase-2 second-pass scientific validation added

The second focused validation pass adds compact ensemble-level checks without
expanding the package API:

- test-local trajectory summaries derive final size and transmission counts from
  `EventLog`
- `sellke` and `gillespie` are compared in a shared finite-population setting
  using ensemble mean final size and mean transmission counts
- a shared zero-transmission limiting case checks that both active engines
  produce no transmissions and only resolve the initially infected/exposed hosts
- a monotonicity benchmark checks that higher transmission produces larger
  ensemble mean final size and more transmissions in a simple finite-population
  setting

These tests are stochastic but deliberately summary-based, not exact-path
comparisons. They are intended to catch large scientific regressions while
remaining fast enough for the normal test suite. They do not replace future
exact finite-state validation against an external CTMC probability propagator.

---

## Phase-2 exact finite-state validation added

The next validation passes add exact finite-state transient comparisons for a
tiny SEIR CTMC matching the active engines with no sampling:

- finite population size `N = 4`
- initial state `(S, E, I) = (3, 0, 1)`
- infection rate `β S I / N`
- activation rate `α E`
- recovery rate `γ I`
- fixed observation time `t = 1.2`

The exact reference distribution over `(S, E, I)` states is computed in the test
suite by finite-state uniformization, the same class of reference calculation
provided by packages such as `DifferentiatedUniformization.jl`. This remains a
validation-only calculation and does not add a runtime dependency or public API.

The first exact test compared `gillespie` ensemble estimates against exact
transient values for:

- expected infectious count `E[I(t)]`
- extinction probability `P(I(t) = 0)`

The strengthened exact test now compares the full infectious-count marginal
distribution:

- `P(I(t) = k)` for every `k = 0, ..., N`
- total variation distance between empirical and exact marginals
- maximum elementwise probability error

The same exact reference validates both `gillespie` and `sellke`. For `sellke`,
the Markovian overlap uses constant `β`, exponential exposed and infectious
durations with means `1 / α` and `1 / γ`, and a degenerate infinite sampling time
so that recovery is the only infectious removal event.

This is the first distribution-level exact reference comparison for `EpiSim.jl`.
It complements, rather than replaces, the stochastic cross-engine and
monotonicity checks above. Future work can cross-check the in-test
uniformization calculation directly against `DifferentiatedUniformization.jl` if
that package is available in the validation environment.
