"""
    initialize_linelist(N::Int64; type::Union{Nothing, Vector{Int64}}=nothing)::DataFrame

Initializes a linelist DataFrame to store event details for each individual in the simulation.

# Arguments
- `N::Int64`: The initial population size.
- `type::Union{Nothing, Vector{Int64}}`: An optional vector specifying the types of individuals. If `nothing`, all individuals are assigned a type of 1.

# Returns
- `DataFrame`: A DataFrame with the following columns:
    - `event`: Event type (initialized to 1 for all individuals).
    - `child_id`: Unique ID for each individual.
    - `child_type`: Type of the individual.
    - `parent_id`: Parent ID (initialized to 0 for all individuals).
    - `parent_type`: Type of the parent (initialized to 0 for all individuals).
    - `t_birth`: Time of birth/infection (initialized to 0.0 for all individuals).
    - `t_death`: Time of death/recovery (initialized to Inf for all individuals).
    - `t_sam`: Time of sampling (initialized to -1.0 for all individuals).

# Example
```julia
linelist = initialize_linelist(10)
```
"""
function initialize_linelist(N::Int64; type::Union{Nothing, Vector{Int64}}=nothing)::DataFrame
    type = isnothing(type) ? fill(1,N) : type
    return DataFrame(
        event = fill(1, N),
        child_id = collect(1:N),
        child_type = type,
        parent_id = fill(0, N),
        parent_type = fill(0, N),
        t_birth = fill(0., N),
        t_death = fill(Inf, N),
        t_sam = fill(-1., N))
end


"""
    initialize_active(N::Vector{Int64})::Tuple{Vector{Vector{Int64}}, Vector{Int64}}

Initializes the active individuals list and their types for the simulation.

# Arguments
- `N::Vector{Int64}`: A vector where each element represents the initial population size of a specific type.

# Returns
- `Tuple{Vector{Vector{Int64}}, Vector{Int64}}`: A tuple containing:
    - `active::Vector{Vector{Int64}}`: A vector of vectors, where each inner vector contains the IDs of active individuals of a specific type.
    - `type::Vector{Int64}`: A vector containing the type of each individual.

# Example
```julia
active, type = initialize_active([10, 5])
```
"""
function initialize_active(N::Vector{Int64})::Tuple{Vector{Vector{Int64}}, Vector{Int64}}
    N_tot = 0
    active = [Vector{Int64}() for _ in 1:length(N)]
    type = Vector{Int64}()
    for t in eachindex(N)
        for _ in 1:N[t]
            N_tot += 1
            push!(active[t], N_tot)
            push!(type, t)
        end
    end
    return active, type
end


"""
    check_parameters(params::CRBDParameters)

Checks the validity of the parameters for the CRBD (Constant Rate Birth-Death) model.

# Arguments
- `params::CRBDParameters`: A `CRBDParameters` struct containing the parameters for the CRBD model.

# Ensures
- `N₀` is non-negative.
- `λ` (birth rate) is non-negative.
- `μ` (death rate) is non-negative.
- `ψ` (extinct sampling rate) is non-negative.
- `ρ₀` (extant sampling rate) is between 0 and 1.
- `r` (removal probability) is between 0 and 1.
- `t_max` (maximum simulation time) is positive.

# Errors
Throws an `AssertionError` if any parameter is invalid, with a message indicating which parameter is out of range and its value.
"""
function check_parameters(params::CRBDParameters)
    @unpack N₀, λ, μ, ψ, ρ₀, r, t_max = params
    # Ensure input parameters are valid
    @assert N₀ ≥ 0 "Initial population size (N₀) must be non-negative"
    @assert λ ≥ 0. "Birth rate (λ) must be non-negative"
    @assert μ ≥ 0. "Death rate (μ) must be non-negative"
    @assert ψ ≥ 0. "Extinct sampling rate (ψ) must be non-negative"
    @assert 0. ≤ ρ₀ ≤ 1. "Extant sampling rate (ρ₀) must be between 0 and 1"
    @assert 0. ≤ r ≤ 1. "Removal probability (r) must be between 0 and 1"
    @assert t_max > 0 "Maximum simulation time (t_max) must be positive"
end


"""
    check_parameters(params::MTBDParameters)

Checks the validity of the parameters for the MTBD (Multi-Type Birth-Death) model.

# Arguments
- `params::MTBDParameters`: A `MTBDParameters` struct containing the parameters for the MTBD model.

# Ensures
- Parameters have consistent dimensions.
- `N₀` is non-negative.
- `λ` (birth rates) are non-negative.
- `μ` (death rates) are non-negative.
- `γ` (transition rates) have consistent dimensions.
- `ψ` (extinct sampling rates) are non-negative.
- `ρ₀` (extant sampling rates) are between 0 and 1.
- `r` (removal probabilities) are between 0 and 1.
- `t_max` (maximum simulation time) is positive.

# Errors
Throws an `AssertionError` if any parameter is invalid, with a message indicating which parameter is out of range and its value.
"""
function check_parameters(params::MTBDParameters)
    @unpack n_types, N₀, λ, μ, γ, ψ, ρ₀, r, t_max = params
    @assert size(λ)[1] == size(λ)[2] == length(μ) == size(γ)[1] == size(γ)[2] == length(ψ) == length(ρ₀) == length(r) "Parameters must have consistent dimensions (=n_types)"
    @assert all(N₀ .≥ 0) "Initial population size (N₀) must be non-negative"
    @assert all(λ .≥ 0.) "Birth rate (λ) must be non-negative"
    @assert all(μ .≥ 0.) "Death rate (μ) must be non-negative"
    @assert all(ψ .≥ 0.) "Extinct sampling rate (ψ) must be non-negative"
    @assert all(0. .≤ ρ₀ .≤ 1.) "Extant sampling rate (ρ₀) must be between 0 and 1"
    @assert all(0. .≤ r .≤ 1.) "Removal probability (r) must be between 0 and 1"
    @assert t_max > 0 "Maximum simulation time (t_max) must be positive"
end