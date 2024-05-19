
"""
    summarize(parms::CRBDParameters)

Prints a summary of key parameters for the Constant Rate Birth-Death (CRBD) model.

# Arguments
- `parms::CRBDParameters`: A `CRBDParameters` struct containing the parameters for the CRBD model.

# Output
Prints the following summaries to the console:
- Lifespan: Expected lifespan of an individual, calculated as the reciprocal of the removal rate.
- Growth rate: Net growth rate, calculated as the difference between the birth rate and the removal rate.
- R₀: Basic reproduction number, calculated as the birth rate divided by the removal rate.
- Sample rate: Rate at which individuals are sampled, calculated as the extinct/ancestral sampling rate divided by the removal rate.

# Example
```julia
params = CRBDParameters(
    N₀=10,
    λ=0.5,
    μ=0.2,
    ψ=0.1,
    ρ₀=0.1,
    r=0.1,
    t_max=100.0
)
summarize(params)
```
"""
function summarize(parms::CRBDParameters)
    removal_rate = parms.μ + parms.r * parms.ψ
    println("Lifespan   : ", 1. / removal_rate)
    println("Growth rate: ", parms.λ - removal_rate)
    println("R₀         : ", parms.λ / removal_rate)
    print("Sample rate: ", parms.ψ / removal_rate)
end


# TODO: Update to incorporate mutation (currently not included)
"""
    summarize(parms::MTBDParameters)

Prints a summary of key parameters for the Multi-Type Birth-Death (MTBD) model.

# Arguments
- `parms::MTBDParameters`: An `MTBDParameters` struct containing the parameters for the MTBD model.

# Output
Prints the following summaries to the console:
- Lifespan: Expected lifespan of individuals of each type, calculated as the reciprocal of the removal rate.
- Growth rate: Net growth rate for each type, calculated as the difference between the birth rate and the removal rate.
- Rᵢ: Basic reproduction number for each type, calculated as the birth rate divided by the removal rate.
- R₀: Overall basic reproduction number considering the contributions of all types.
- Sample rate: Rate at which individuals are sampled for each type, calculated as the extinct/ancestral sampling rate divided by the removal rate.
- Frequencies: The equilibrium frequencies of the types.

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
summarize(params)
```
"""
function summarize(parms::MTBDParameters)
    removal_rate = parms.μ .+ parms.r .* parms.ψ
    println("Lifespan   : ", 1. ./ removal_rate)
    println("Growth rate: ", parms.λ .- removal_rate)
    Rᵢ = parms.λ ./ removal_rate
    println("Rᵢ         : ", Rᵢ)
    Λ = parms.λ[1, 1] - parms.λ[2,2] - removal_rate[1] + removal_rate[2]
    c = sqrt(Λ^2 + 4 * parms.λ[1,2] * parms.λ[2,1])
    f₁ = (c + Λ) / (c + Λ + 2. * parms.λ[1,2])
    println("R₀         : ", f₁ * (Rᵢ[1,1] + Rᵢ[1, 2]) + (1. - f₁) * (Rᵢ[2,1] + Rᵢ[2, 2]))
    println("Sample rate: ", parms.ψ ./ removal_rate)
    print("Frequencies :", [f₁, 1. - f₁])
end


"""
    n_sampled(linelist::DataFrame) -> Int

Calculates the number of sampled individuals in the linelist.

# Arguments
- `linelist::DataFrame`: A DataFrame containing the linelist of individuals, with a column `:t_sam` indicating the sampling time.

# Returns
- `Int`: The number of sampled individuals, where sampling time (`:t_sam`) is greater than 0.

# Example
```julia
using DataFrames

linelist = DataFrame(
    t_sam = [0.0, -1.0, 1.0, 2.0, -1.0]
)

n = n_sampled(linelist)  # n will be 2
```
"""
function n_sampled(linelist::DataFrame)
    sampled = 0
    @eachrow! linelist begin
        if :t_sam > 0.
            sampled += 1
        end
    end
    return sampled
end

@forward Outbreak.linelist n_sampled


"""
    n_deceased(linelist::DataFrame) -> Int

Calculates the number of deceased individuals in the linelist.

# Arguments
- `linelist::DataFrame`: A DataFrame containing the linelist of individuals, with a column `:t_death` indicating the time of death.

# Returns
- `Int`: The number of deceased individuals, where time of death (`:t_death`) is greater than 0 and finite.

# Example
```julia
using DataFrames

linelist = DataFrame(
    t_death = [0.0, -1.0, 1.0, 2.0, Inf]
)

n = n_deceased(linelist)  # n will be 2
```
"""
function n_deceased(linelist::DataFrame)
    deceased = 0
    @eachrow! linelist begin
        if :t_death > 0. && isfinite(:t_death)
            deceased += 1
        end
    end
    return deceased
end


@forward Outbreak.linelist n_deceased


"""
    offspring_dist(linelist::DataFrame; height=Inf, deceased_only=true) -> Vector{Int}

Calculates the offspring distribution for individuals in the linelist.

# Arguments
- `linelist::DataFrame`: A DataFrame containing the linelist of individuals, with columns `:parent_id`, `:child_id`, and `:t_death`.
- `height::Float64`: A threshold time. Only individuals with `:t_death` less than this value will be considered deceased. Defaults to `Inf`.
- `deceased_only::Bool`: If `true`, only includes individuals who are deceased (determined by `:t_death < height`). If `false`, includes all individuals. Defaults to `true`.

# Returns
- `Vector{Int}`: A vector where each element represents the number of offspring for each individual, filtered by the `deceased_only` and `height` criteria.

# Example
```julia
using DataFrames

linelist = DataFrame(
    parent_id = [0, 1, 1, 2, 3],
    child_id = [1, 2, 3, 4, 5],
    t_death = [Inf, 0.0, 1.0, 2.0, Inf]
)

dist = offspring_dist(linelist, height=1.0, deceased_only=true)  # dist will be [2, 1]
```
"""
function offspring_dist(linelist::DataFrame; height=Inf, deceased_only=true)
    dist = fill(0, nrow(linelist))
    deceased = fill(false, nrow(linelist))
    @eachrow! linelist begin
        if :parent_id > 0
            dist[:parent_id] += 1
        end
        deceased[:child_id] = (:t_death < height || !deceased_only) ? true : false
    end
    return dist[deceased]
end


@forward Outbreak.linelist offspring_dist


function offspring_dist(out; deceased_only)
    return offspring_dist(out.linelist, height=out.height, deceased_only=deceased_only)
end


function offspring_dist(out::Outbreak; deceased_only=true)
    @unpack linelist, height = out
    dist = fill(0, nrow(linelist))
    deceased = fill(false, nrow(linelist))
    @eachrow! linelist begin
        if :parent_id > 0
            dist[:parent_id] += 1
        end
        deceased[:child_id] = (:t_death < height || !deceased_only) ? true : false
    end
    return dist[deceased]
end


"""
    type_dist(out::Outbreak) -> Vector{Int}

Calculates the distribution of types in the outbreak.

# Arguments
- `out::Outbreak`: An `Outbreak` struct containing the outbreak data, including parameters and the linelist.

# Returns
- `Vector{Int}`: A vector where each element represents the number of individuals of each type.

# Example
```julia
using DataFrames

# Example Outbreak struct
struct Outbreak
    parms
    traj
    linelist::DataFrame
end

# Example linelist DataFrame
linelist = DataFrame(
    child_type = [1, 2, 1, 2, 1]
)

# Example parms with n_types
parms = (n_types=2, )

# Create an outbreak
outbreak = Outbreak(parms, nothing, linelist)

# Get the type distribution
dist = type_dist(outbreak)  # dist will be [3, 2]
```
"""
function type_dist(out::Outbreak)
    n_types = hasproperty(out.parms, :n_types) ? out.parms.n_types : 1
    dist = fill(0, n_types)
    @eachrow! out.linelist begin
        dist[:child_type] += 1
    end
    return dist
end


"""
    summarize(linelist::DataFrame, t::Float64) -> Tuple{Vector{Int}, Vector{Bool}, Vector{Bool}, Vector{Float64}, Vector{Int}}

Summarizes the linelist data, providing information about the type, lifespan, offspring count, and status of individuals.

# Arguments
- `linelist::DataFrame`: A DataFrame containing the linelist of individuals, with columns `:child_id`, `:child_type`, `:t_death`, `:t_sam`, and `:t_birth`.
- `t::Float64`: The current time for the simulation, used to calculate lifespan for individuals still alive.

# Returns
- `Tuple{Vector{Int}, Vector{Bool}, Vector{Bool}, Vector{Float64}, Vector{Int}}`:
  - `type::Vector{Int}`: Vector indicating the type of each individual.
  - `deceased::Vector{Bool}`: Vector indicating whether each individual is deceased.
  - `sampled::Vector{Bool}`: Vector indicating whether each individual is sampled.
  - `lifespan::Vector{Float64}`: Vector indicating the lifespan of each individual.
  - `offspring::Vector{Int}`: Vector indicating the number of offspring for each individual.

# Example
```julia
using DataFrames

linelist = DataFrame(
    child_id = [1, 2, 3, 4, 5],
    child_type = [1, 2, 1, 2, 1],
    t_death = [Inf, 0.0, 1.0, 2.0, Inf],
    t_sam = [0.0, -1.0, 1.0, 2.0, -1.0],
    t_birth = [0.0, 0.0, 0.0, 0.0, 0.0]
)
t = 3.0

type, deceased, sampled, lifespan, offspring = summarize(linelist, t)
```
"""
function summarize(linelist::DataFrame, t::Float64)
    n = nrow(linelist)
    type = fill(-1, n)
    lifespan = zeros(n)
    offspring = fill(0, n)
    deceased = fill(false, n)
    sampled = fill(false, n)

    type[1] = linelist[1, :child_type]
    deceased[1] = !isinf(linelist[1, :t_death])
    sampled[1] = linelist[1, :t_sam] > 0.
    lifespan[1] = deceased[1] ? linelist[1, :t_death] - linelist[1, :t_birth] : t - linelist[1, :t_birth]

    @eachrow! linelist[2:end, :] begin
        type[:child_id] = :child_type
        deceased[:child_id] = !isinf(:t_death)
        sampled[:child_id] = :t_sam > 0.
        lifespan[:child_id] = deceased[:child_id] ? :t_death - :t_birth : t - :t_birth
        offspring[:parent_id] += 1
    end
    return type, deceased, sampled, lifespan, offspring
end


"""
    summarize(out::Outbreak)


    OutbreakSummary

A struct to represent the summary of an epidemic outbreak, including various rates and metrics.

# Fields
- `δbar::Vector{Float64}`: Average death rates for each type.
- `λbar::Matrix{Float64}`: Average birth rates between types.
- `ψbar::Vector{Float64}`: Average sampling rates for each type.
- `Ribar::Matrix{Float64}`: Basic reproduction number matrix for each type pair.
- `Rbar::Float64`: Overall basic reproduction number.
- `frequency::Vector{Float64}`: Equilibrium frequencies of the types.
- `sampled::Vector{Int64}`: Number of sampled individuals for each type.

# Example
```julia
using EpiSim
using DataFrames

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
outbreak = simulate(params, N_max=1000, S_max=100)

# Summarize the outbreak
summary = summarize(outbreak)
```
"""
function summarize(out::Outbreak)
    n_types = typeof(out.parms) == MTBDParameters ? out.parms.n_types : 1
    n = nrow(out.linelist)
    end_time = out.traj[1, end]

    sampled = zeros(n_types)
    deceased = zeros(n_types)
    frequency = zeros(n_types)


    # Track cumulative lifespans and offspring
    lifespan = zeros(n_types) .+ eps()
    offspring = zeros(n_types, n_types)
    mutations = zeros(n_types, n_types)

    @eachrow! out.linelist[2:end, :] begin
        if :parent_id > 0
            frequency[:child_type] += 1.
            if :event == 1
                offspring[:parent_type, :child_type] += 1
            elseif :event == 2
                mutations[:parent_type, :child_type] += 1
            end
        end

        if !isinf(:t_death)
            deceased[:child_type] += 1
            lifespan[:child_type] += :t_death - :t_birth
        else
            lifespan[:child_type] += end_time - :t_birth
        end

        if :t_sam > 0.
            sampled[:child_type] += 1
        end
    end

    # Calculate expected lifespan and offspring
    δbar = deceased ./ lifespan .+ eps()
    λbar = offspring ./ lifespan
    ψbar = sampled ./ lifespan
    Ribar = λbar ./ δbar
    frequency ./= sum(frequency)
    Rbar = frequency' * Ribar * ones(n_types)

    return  OutbreakSummary(δbar, λbar, ψbar, Ribar, Rbar, frequency, sampled)
end


@forward Outbreak.out summarize