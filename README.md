# EpiSim.jl

EpiSim.jl is a Julia package designed for stochastic simulation of epidemic outbreaks. The package allows for the simulation of different epidemic models, including compartmental models, renewal models, and birth-death models. It returns detailed event diaries (linelists) for individual members of the population, which can be used to generate phylogenetic trees.

## Installation

To install EpiSim.jl, use the following command in the Julia REPL:
```julia
using Pkg
Pkg.add(url="https://github.com/michaeltmeehan/EpiSim.jl.git")
```
## Usage

### Basic Usage

To simulate an epidemic outbreak, you can use the provided model functions. Here is an example using the constant-rate birth-death (CRBD) model:
```julia
using EpiSim

# Define parameters for the CRBD model
params = CRBDParameters(
    N₀=10,
    λ=0.5,
    μ=0.2,
    ψ=0.1,
    ρ₀=0.1,
    r=0.1,
    t_max=100.0
)

# Simulate the outbreak
outbreak = simulate_outbreak(params, N_max=1000, S_max=100)

# Summarize the outbreak
summary = summarize(outbreak)
```
### Multi-Type Birth-Death (MTBD) Model

To simulate an outbreak using the multi-type birth-death (MTBD) model:
```julia
using EpiSim

# Define parameters for the MTBD model
params = MTBDParameters(
    n_types=2,
    N₀=[10, 5],
    λ=[0.5 0.3; 0.2 0.4],
    μ=[0.1, 0.2],
    γ=[0.0 0.1; 0.1 0.0],
    ψ=[0.1, 0.1],
    ρ₀=[0.1, 0.1],
    r=[0.1, 0.1],
    t_max=50.0
)

# Simulate the outbreak
outbreak = simulate_outbreak(params, N_max=1000, S_max=100)

# Summarize the outbreak
summary = summarize(outbreak)
```
### Helper Functions

EpiSim.jl provides several utility functions to aid in analyzing and processing simulation results. Some examples include:

- n_sampled(linelist::DataFrame): Returns the number of sampled individuals.
- n_deceased(linelist::DataFrame): Returns the number of deceased individuals.
- offspring_dist(linelist::DataFrame; height=Inf, deceased_only=true): Returns the distribution of offspring for each individual.
- type_dist(out::Outbreak): Returns the distribution of types in the outbreak.

### Testing

Unit tests are provided to ensure the correctness of the implementation. To run the tests, use the following command in the Julia REPL:
```julia
using Pkg
Pkg.test("EpiSim")
```
## Contributing

Contributions to EpiSim.jl are welcome. Please feel free to open issues or submit pull requests on the GitHub repository: https://github.com/michaeltmeehan/EpiSim.jl

## License

This project is licensed under the MIT License - see the LICENSE file for details.

