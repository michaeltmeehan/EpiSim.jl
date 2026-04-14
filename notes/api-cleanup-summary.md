# API cleanup summary

## Exports removed or corrected

- Removed stale event export `EK_Sampling`; the active enum defines `EK_FossilizedSampling` and `EK_SerialSampling` instead.
- Removed stale Sellke-related exports `Traits`, `isdone`, `make_infected`, `next_event`, and `next_time`; these names are not defined by the active included files.
- Removed stale tree node exports `Node`, `Root`, `Binary`, `SampledLeaf`, `SampledUnary`, and `UnsampledUnary`; the active tree representation is the columnar `Tree` plus `NodeKind`.
- Removed stale processing exports `get_seeds`, `get_subtree`, `get_sampled_events`, `get_subtrees`, and `get_sampled_subtrees`; their active implementations are not present.
- Corrected the active birth-death likelihood method to use the defined `NK_*` node-kind constants rather than stale `K_*` names.
- Removed unused top-level imports from `src/EpiSim.jl` so the active source no longer directly imports non-core plotting/dataframe dependencies.

## Exports retained

- Random helpers: `popr!`, `wsample`, `wsampleindex`, `wsampleindex_cols`.
- Event-log API: `EventLog`, `EventKind`, `EK_None`, `EK_Seeding`, `EK_Transmission`, `EK_FossilizedSampling`, `EK_SerialSampling`, `EK_Recovery`, `EK_Activation`.
- Simulator API: `sellke`, `TraitDists`, `gillespie`.
- Tree API: `Tree`, `NodeKind`, `NK_None`, `NK_Root`, `NK_Binary`, `NK_SampledLeaf`, `NK_UnsampledLeaf`, `NK_SampledUnary`, `NK_UnsampledUnary`.
- Tree extraction: `extract_sampled_trees`.
- Birth-death helpers: `likelihood`, `p₀`, `γ`, `γ₀`, `γ₁`, `β`.

## Ambiguous names left unchanged

- `TraitDists` remains exported because it is part of the active Sellke input normalization path and was already exported.
- `NodeKind` and `Tree` remain exported because they are consumed by active tree extraction and birth-death likelihood code.
- Internal Sellke implementation details such as `StateKind`, `Population`, `ActiveList`, and scheduling helpers remain unexported.
