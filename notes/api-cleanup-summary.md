# API cleanup summary

> Superseded note: this file records an earlier API cleanup stage. A later
> legacy-code cleanup removed the remaining internal tree extraction,
> birth-death likelihood, and generic model/simulation remnants from
> `EpiSim.jl`. The retained package surface is now the random helpers,
> event-log API, `sellke`, and `gillespie`.

## Exports removed or corrected

- Removed stale event export `EK_Sampling`; the active enum defines `EK_FossilizedSampling` and `EK_SerialSampling` instead.
- Removed stale Sellke-related exports `Traits`, `isdone`, `make_infected`, `next_event`, and `next_time`; these names are not defined by the active included files.
- Removed stale tree node exports `Node`, `Root`, `Binary`, `SampledLeaf`, `SampledUnary`, and `UnsampledUnary`; later cleanup removed the remaining internal tree representation as out of scope.
- Removed stale processing exports `get_seeds`, `get_subtree`, `get_sampled_events`, `get_subtrees`, and `get_sampled_subtrees`; their active implementations are not present.
- Earlier cleanup corrected the birth-death likelihood method, but later cleanup removed birth-death likelihood code from `EpiSim.jl` as out of scope.
- Removed unused top-level imports from `src/EpiSim.jl` so the active source no longer directly imports non-core plotting/dataframe dependencies.

## Exports retained

- Random helpers: `popr!`, `wsample`, `wsampleindex`, `wsampleindex_cols`.
- Event-log API: `EventLog`, `EventKind`, `EK_None`, `EK_Seeding`, `EK_Transmission`, `EK_FossilizedSampling`, `EK_SerialSampling`, `EK_Recovery`, `EK_Activation`.
- Simulator API: `sellke`, `TraitDists`, `gillespie`.

## Ambiguous names left unchanged

- `TraitDists` remains exported because it is part of the active Sellke input normalization path and was already exported.
- Internal Sellke implementation details such as `StateKind`, `Population`, `ActiveList`, and scheduling helpers remain unexported.
