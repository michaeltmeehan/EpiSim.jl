# EpiSim.jl

EpiSim.jl is a Julia package for stochastic epidemic outbreak simulation and event-log generation.

The current recovered package scope is intentionally narrow:

- stochastic epidemic simulation through the active `sellke` and `gillespie` engines
- semantically documented event logs
- validation utilities for event-log consistency

The package does not provide tree extraction, tree likelihoods, birth-death analytical likelihoods, or orchestration logic for the broader modelling ecosystem.

## Installation

To install EpiSim.jl, use the following command in the Julia REPL:
```julia
using Pkg
Pkg.add(url="https://github.com/michaeltmeehan/EpiSim.jl.git")
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.
