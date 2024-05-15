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
