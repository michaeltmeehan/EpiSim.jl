# Load and test summary

## Imports changed

- Kept only active top-level imports in `src/EpiSim.jl`: `DataStructures`, `Distributions`, `Random`, and `UnPack`.
- Removed the unused `Base.Threads` import.
- Replaced deprecated `DataStructures.top` usage in the Sellke engine with `first`, matching the current heap API without changing simulator behavior.

## Non-core dependencies

- Removed direct `Project.toml` dependencies that are not used by the active included code path: plotting, dataframe, benchmarking, sparse/linear-algebra helper, and legacy utility packages.
- Kept active runtime dependencies: `DataStructures`, `Distributions`, `Random`, and `UnPack`.
- Moved `Test` to test-only metadata via `[extras]` and `[targets]`.

## Tests added

- `EventLog` constructor and field-length invariants.
- `update_event_log!` append behavior.
- Sellke smoke/invariant test.
- Gillespie smoke/invariant test.
- Small sampled-tree extraction smoke test.

## Remaining confidence limits

- Tests are intentionally smoke-level and do not validate statistical correctness.
- Tree extraction is covered only by a tiny deterministic event log.
- Birth-death likelihood has no numerical reference-value tests yet.
