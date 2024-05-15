
"""
# Outbreak

A struct to represent the outcome of an epidemic simulation.

## Fields
- `parms::EpiParameters`: The parameters used for the epidemic simulation.
- `traj::Matrix{Float64}`: The trajectory of the epidemic over time, where each column represents a time point and rows represent the prevalence of sub-populations.
- `linelist::DataFrame`: A detailed event diary for each individual in the simulation, including times of infection/birth, recovery/death, and sampling events.

## Example
```julia
using EpiSim

# Example parameters (for demonstration purposes)
params = CRBDParameters(N₀=10, λ=0.5, μ=0.2, ψ=0.1, ρ₀=0.1, r=0.1, t_max=100.0)

# Simulate an outbreak
outbreak = simulate(params, N_max=1000, S_max=100)
```
"""
struct Outbreak
    parms::EpiParameters
    traj::Matrix{Float64}
    linelist::DataFrame
end


"""
    show(io::IO, out::Outbreak)

Custom show method for the Outbreak struct.

Prints a summary of the outbreak, including:
- Initial population size
- Final population size and cumulative count
- Number of sampled individuals
- Time span of the outbreak

# Example
```julia
params = CRBDParameters(N₀=10, λ=0.5, μ=0.2, ψ=0.1, ρ₀=0.1, r=0.1, t_max=100.0)
outbreak = simulate(params, N_max=1000, S_max=100)
println(outbreak)
```
"""
function Base.show(io::IO, out::Outbreak)
    println(io, "Outbreak Summary")
    println(io, "================")
    println(io, "Initial pop.     : ", out.parms.N₀)
    println(io, "Final pop. (cum) : ", out.traj[2:end, end], " (", nrow(out.linelist), ")")
    println(io, "Sampled          : ", n_sampled(out))
    print(io, "Time span        : ", round(out.traj[1, end] - out.traj[1,1], digits=2), " time units")
end