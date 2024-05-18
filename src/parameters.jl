"""
    EpiParameters

Abstract type for all epidemic parameters. This serves as the base type for all parameter sets used in epidemic simulations.
"""
abstract type EpiParameters end

"""
    BirthDeathParameters <: EpiParameters

Abstract type for birth-death model parameters. Inherits from `EpiParameters`.
"""
abstract type BirthDeathParameters <: EpiParameters end


"""
    CRBDParameters

A type for storing parameters of the Constant Rate Birth-Death (CRBD) model.

## Fields
- `N₀::Int64`: Initial population size.
- `λ::Float64`: Birth rate.
- `μ::Float64`: Death rate.
- `ψ::Float64`: Extinct/ancestral sampling rate.
- `ρ₀::Float64`: Extant sampling rate.
- `r::Float64`: Removal probability (upon sampling).
- `t_max::Float64`: Maximum simulation time.
"""
@with_kw mutable struct CRBDParameters <: BirthDeathParameters
    N₀::Int64 = 1   # initial population size
    λ::Float64      # birth rate
    μ::Float64      # death rate
    ψ::Float64      # extinct / ancestral sampling rate
    ρ₀::Float64     # extant sampling rate
    r::Float64      # removal probability (upon sampling)
    t_max::Float64  # maximum simulation time
end


"""
    MTBDParameters

Parameters for the Multi-Type Birth-Death (MTBD) model.

# Fields
- `n_types::Int64`: Number of distinct subgroups/types.
- `N₀::Vector{Int64}`: Initial population distribution of each type.
- `λ::Matrix{Float64}`: Birth rates where `λ[a, b]` represents the rate of type `a` giving birth to type `b`.
- `μ::Vector{Float64}`: Death rates for each type.
- `γ::Matrix{Float64}`: Mutation rates where `γ[a, b]` represents the rate of type `a` mutating to type `b`.
- `ψ::Vector{Float64}`: Extinct/ancestral sampling rates for each type.
- `ρ₀::Vector{Float64}`: Extant sampling rates for each type.
- `r::Vector{Float64}`: Removal probabilities upon sampling for each type.
- `t_max::Float64`: Maximum simulation time.

# Example
```julia
params = MTBDParameters(
    n_types=2,
    N₀=[1, 0],
    λ=[0.5 0.6; 1. 2.],
    μ=[0.2, 0.3],
    γ=[0. 0.1; 0.2 0.],
    ψ=[0.1, 0.1],
    ρ₀=[0.1, 0.2],
    r=[0.1, 0.1],
    t_max=100.0
)
```
"""
@with_kw mutable struct MTBDParameters <: BirthDeathParameters
    n_types::Int64          # number of distinct subgroups / types
    N₀::Vector{Int64}       # initial population distribution of each type a
    λ::Matrix{Float64}      # birth rate of type a -> b
    μ::Vector{Float64}      # death rate
    γ::Matrix{Float64}      # mutation rate of type a -> b
    ψ::Vector{Float64}      # extinct / ancestral sampling rate
    ρ₀::Vector{Float64}     # extant sample rate
    r::Vector{Float64}      # removal probability (upon sampling)
    t_max::Float64  # maximum simulation time
end


"""
    CompartmentalParameters <: EpiParameters

Abstract type for compartmental model parameters. Inherits from `EpiParameters`.
"""
abstract type CompartmentalParameters <: EpiParameters end




"""
    RenewalParameters <: EpiParameters

Abstract type for renewal model parameters. Inherits from `EpiParameters`.
"""
abstract type RenewalParameters <: EpiParameters end

"""
    summary(params::EpiParameters)

Prints a summary of the parameters for the given parameter set.
"""
function summary(params::EpiParameters)
    println("Parameter Summary:")
    for field in fieldnames(typeof(params))
        println("$(field): $(getfield(params, field))")
    end
end

"""
    update_parameters!(params::EpiParameters; kwargs...)

Updates the fields of the given parameter set with new values provided as keyword arguments.
"""
function update_parameters!(params::EpiParameters; kwargs...)
    for (key, value) in kwargs
        if hasfield(typeof(params), key)
            setfield!(params, key, value)
        else
            println("Warning: $(key) is not a valid field for $(typeof(params))")
        end
    end
end


struct OutbreakSummary
    δbar::Vector{Float64}
    λbar::Matrix{Float64}
    ψbar::Vector{Float64}
    Ribar::Matrix{Float64}
    Rbar::Float64
    frequency::Vector{Float64}
    sampled::Vector{Int64}
end


function Base.show(io::IO, summ::OutbreakSummary)
    println(io, "Sampled          : ", summ.sampled)
    println(io, "Removal rate     : ", round.(summ.δbar, digits=3))
    println(io, "Birth rate       : ", round.(summ.λbar, digits=3))
    println(io, "Sample rate      : ", round.(summ.ψbar, digits=3))
    println(io, "Sample prop.     : ", round.(summ.ψbar ./ summ.δbar, digits=3))
    println(io, "R (type)         : ", round.(summ.Ribar, digits=3))
    println(io, "R (total)        : ", round.(summ.Rbar, digits=3))
    print(io, "Frequency        : ", round.(summ.frequency, digits=3))
end