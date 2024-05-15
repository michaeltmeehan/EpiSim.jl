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
