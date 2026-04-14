# Trust Criteria

## Role of this package

`EpiSim.jl` is the standalone stochastic epidemic simulation and event-log package in the recovered outbreak-modelling ecosystem.

Its current role is to provide trustworthy simulation outputs and semantically interpretable event logs without depending on downstream tree or analytical packages.

---

## Trust goal

`EpiSim.jl` is trustworthy when its current simulation core can be used to generate reproducible and scientifically interpretable outbreak simulations within the documented scope.

---

## What trust does mean here

Trust in this phase means:

- the package’s current simulation scope is explicit
- the event-log format has documented semantics
- event logs satisfy tested invariants
- engines behave consistently in controlled comparable settings
- seeded reproducibility works as documented
- edge cases fail clearly or behave in explicitly documented ways

---

## What trust does not mean here

Trust in this phase does **not** imply:

- correctness for every conceivable epidemic model extension
- built-in tree extraction as part of the stable public API
- birth-death likelihood functionality
- ecosystem orchestration
- full agreement across engines in settings where they are not expected to coincide

---

## Stable trust boundary

The current trust boundary includes:

- the stabilised active simulation core
- public simulator functionality intentionally retained
- event-log generation functionality intentionally retained

Anything outside this boundary should be considered provisional unless explicitly promoted into the stable API.

---

## Conditions required for trust

### 1. Semantic clarity

The event log must have explicit documented meaning.

### 2. Internal consistency

Generated logs and simulation outputs must satisfy documented invariants.

### 3. Controlled agreement

Comparable engines/settings must agree to the extent scientifically expected.

### 4. Reproducibility

Seeded workflows must behave deterministically where promised.

### 5. Boundary safety

Extreme or invalid scenarios must not silently produce misleading outputs.

---

## Known trust risks

Current or likely risks include:

- under-specified event semantics
- subtle divergence across engines
- edge-case failures in stopping rules or rare-event paths
- hidden reproducibility issues
- mismatch between simulated state transitions and logged events

---

## Required evidence before calling this package trustworthy

The following evidence is required:

- a written event-log semantics note
- invariant tests for logs and state consistency
- cross-engine benchmark tests in simple settings
- reproducibility tests under fixed seeds
- stress tests across important parameter regimes

---

## Phase-2 completion standard

For the purposes of this project phase, `EpiSim.jl` is trustworthy when:

1. its current simulation scope is clearly stated
2. its event logs are semantically documented
3. its core invariants are enforced by tests
4. comparable engines have been validated in controlled settings
5. it is reliable enough to serve as the standalone outbreak simulation layer for research workflows
