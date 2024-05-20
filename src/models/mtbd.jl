"""
    simulate_outbreak(N₀::Vector{Int64}, λ::Matrix{Float64}, μ::Vector{Float64}, γ::Matrix{Float64}, ψ::Vector{Float64}, ρ₀::Vector{Float64}, r::Vector{Float64};
             N_max::Union{Float64, Int64}=Inf, S_max::Union{Int64, Float64}=Inf, t_max::Union{Int64, Float64}=Inf)::Outbreak

Simulates an epidemic outbreak using the Multi-Type Birth-Death (MTBD) model and returns the resulting outbreak.

# Arguments
- `N₀::Vector{Int64}`: Initial population distribution across types.
- `λ::Matrix{Float64}`: Birth rates where `λ[a, b]` is the rate of type `a` giving birth to type `b`.
- `μ::Vector{Float64}`: Death rates for each type.
- `γ::Matrix{Float64}`: Mutation rates where `γ[a, b]` is the rate of type `a` mutating to type `b`.
- `ψ::Vector{Float64}`: Sampling rates for each type.
- `ρ₀::Vector{Float64}`: Extant sampling probabilities for each type.
- `r::Vector{Float64}`: Removal probabilities for each type.
- `N_max::Union{Float64, Int64}=Inf`: Maximum population size.
- `S_max::Union{Int64, Float64}=Inf`: Maximum sample size.
- `t_max::Union{Int64, Float64}=Inf`: Maximum simulation time.

# Returns
- `Outbreak`: An `Outbreak` object containing the parameters, trajectory, and linelist of the simulated outbreak.

# Example
```julia
N₀ = [10, 5]
λ = [0.5 0.3; 0.2 0.4]
μ = [0.2, 0.1]
γ = [0.1 0.2; 0.3 0.4]
ψ = [0.1, 0.05]
ρ₀ = [0.1, 0.2]
r = [0.1, 0.1]
outbreak = simulate_outbreak(N₀, λ, μ, γ, ψ, ρ₀, r, N_max=1000, S_max=100, t_max=100.0)
```
"""
function simulate_outbreak(N₀::Vector{Int64},                   # Initial population distribution (across types)
                  λ::Matrix{Float64},                  # Birth rate: a -> b
                  μ::Vector{Float64},                  # Death rate
                  γ::Matrix{Float64},                  # Mutation rate: a -> b
                  ψ::Vector{Float64},                  # Sampling rate
                  ρ₀::Vector{Float64},                 # Extant sampling probability
                  r::Vector{Float64};                  # Removal probability
                  N_max::Union{Float64, Int64} = Inf,  # Maximum population size
                  S_max::Union{Int64, Float64} = Inf,  # Maximum sample size
                  t_max::Union{Int64, Float64} = Inf)::Outbreak

    n_types = size(λ)[1]

    params = MTBDParameters(n_types, N₀, λ, μ, γ, ψ, ρ₀, r, t_max)

    return simulate_outbreak(params, N_max=N_max, S_max=S_max)
end




"""
    simulate_outbreak(params::MTBDParameters; N_max::Union{Float64, Int64}=Inf, S_max::Union{Float64, Int64}=Inf)::Outbreak

Simulates an epidemic outbreak using the Multi-Type Birth-Death (MTBD) model with the given parameters and returns the resulting outbreak.

# Arguments
- `params::MTBDParameters`: A `MTBDParameters` struct containing the parameters for the MTBD model.
- `N_max::Union{Float64, Int64}=Inf`: Maximum population size.
- `S_max::Union{Float64, Int64}=Inf`: Maximum sample size.
- `t_max::Union{Float64, Int64}=Inf`: Maximum simulation time.

# Returns
- `Outbreak`: An `Outbreak` object containing the parameters, trajectory, and linelist of the simulated outbreak.

# Example
```julia
params = MTBDParameters(n_types=2, N₀=[10, 5], λ=[0.5 0.3; 0.2 0.4], μ=[0.2, 0.1], γ=[0.1 0.2; 0.3 0.4], ψ=[0.1, 0.05], ρ₀=[0.1, 0.2], r=[0.1, 0.1], t_max=100.0)
outbreak = simulate_outbreak(params, N_max=1000, S_max=100)
```
"""
function simulate_outbreak(params::MTBDParameters; N_max::Union{Float64, Int64}=Inf, S_max::Union{Float64, Int64}=Inf)::Outbreak

    @unpack n_types, N₀, λ, μ, γ, ψ, ρ₀, r, t_max = params
    check_parameters(params)

    # Determine the number of types
    n_types = length(μ)

    # Specify a maximum population size
    N_max = isinf(min(N_max, S_max, t_max)) ? 1e3 : N_max + 0.

    # Calculate total birth and mutation rates for each type
    Λ = sum(λ, dims=2)
    Γ = sum(γ, dims=2)

    # Throughout the simulation keep track of two key vectors
    # 1. [N]ᵢ - number of individuals of type i
    # 2. [active]ⱼ - ids of active individuals j ∈ 1, ..., ∑_(i=1)^(n_type) Nᵢ

    # Generate the initial distribution of N and active individuals given inputs
    active, type = initialize_active(N₀)
    N_tot = sum(N₀)

    # Initialize linelist of all events
    # - each birth and mutation event appears in a new row
    # - sampling and death update existing entries
    linelist = initialize_linelist(N_tot, type=type)

    # Initialize the simulation
    N = reshape(N₀ .+ 0., (1, n_types))
    t = 0.
    t_out = [t]
    N_out = VectorOfArray([N[:]])

    # Record the cumulative event count
    row_num = N_tot

    # Track cumulative sampled
    S = 0.

    while (N_tot > 0) && (row_num < N_max) && (S < S_max)

        birth_rate = (N * Λ)[1]     # Cumulate birth rate (across all individuals and types)
        death_rate = (N * μ)[1]     # Cumulate death rate (across all individuals and types)
        mutation_rate = (N * Γ)[1]  # Cumulate mutation rate (across all individuals and types)
        sample_rate = (N * ψ)[1]    # Cumulate sample rate (across all individuals and types)

        # Compute the total event rate
        total_event_rate = birth_rate + death_rate + mutation_rate + sample_rate

        # Sample a random value
        rnd = rand()
        t -= log(rnd) / total_event_rate            # Update time
        t > t_max && break                          # Compare time with stopping criteria

        # Birth event
        if rnd <= birth_rate / total_event_rate
            row_num += 1
            parent_type = wsample(1:n_types, (N' .* Λ)[:])              # Sample type of parent
            parent_id = sample(active[parent_type])                     # Sample id of parent
            child_type = wsample(1:n_types, @view λ[parent_type, :])    # Sample type of new active individual
            push!(linelist, (1, row_num, child_type, 
                             parent_id, parent_type, t, Inf, -1.))      # Update linelist
            push!(active[child_type], row_num)                          # Update active
            N[child_type] += 1.                                         # Update population sizes
            N_tot += 1                                                  # Update total popuation size

        # Death event
        elseif rnd <= (birth_rate + death_rate) / total_event_rate
            deceased_type = wsample(1:n_types, (N' .* μ)[:])            # Sample type of deceased
            deceased_idx = sample(1:length(active[deceased_type]))      # Sample index of deceased individual
            deceased = active[deceased_type][deceased_idx]              # Determine id of deceased individual
            linelist[deceased, :t_death] = t                            # Update time of death in linelist
            deleteat!(active[deceased_type], deceased_idx)              # Update list of active individuals
            N[deceased_type] -= 1.                                      # Update population sizes
            N_tot -= 1

        # Mutation event
        elseif rnd <= (birth_rate + death_rate + mutation_rate) / total_event_rate
            row_num += 1                                        
            old_type = wsample(1:n_types, (N' .* Γ)[:])                 # Sample type of mutating individual
            mutant_idx = sample(1:length(active[old_type]))             # Sample index of mutating individual 
            mutant_id = active[old_type][mutant_idx]                    # Retrieve id of mutating individual
            new_type = wsample(1:n_types, @view γ[old_type, :])         # Sample new mutated type
            linelist[mutant_id, :t_death] = t                           # Record death of old individual from linelist
            push!(linelist, (0, row_num, new_type, 
                             mutant_id, old_type, t, Inf, -1.))         # Update linelist with new mutant
            deleteat!(active[old_type], mutant_idx)                     # Update list of active individuals
            push!(active[new_type], row_num)                            # Update list of active individuals
            N[old_type] -= 1.                                           # Update population sizes
            N[new_type] += 1.                                           

        # Sampling event
        else
            sampled_type = wsample(1:n_types, (N' .* ψ)[:])             # Sample type of sampled individual
            sampled_idx = sample(1:length(active[sampled_type]))        # Sample index of sampled individual
            sampled = active[sampled_type][sampled_idx]                 # Retrieve id of sampled individual
            if linelist[sampled, :t_sam] < 0.                           # Check if sampled individual has been sampled previously
                S += 1.                                                 # If not, increment cumulative sample count
            end
            linelist[sampled, :t_sam] = t                               # Update linelist with new sample time

            # Removal
            if rand() < r[sampled_type]                                 # Check if sampled individal is removed
                linelist[sampled, :t_death] = t                         # Update removal time of sampled individual
                deleteat!(active[sampled_type], sampled_idx)            # Update list of active individuals
                N[sampled_type] -= 1.                                   # Update population sizes
                N_tot -= 1                                              # Update total population size
            end
        end

        # Update outputs
        push!(t_out, t)
        push!(N_out, N[:])
    end

    # Present day sampling
    if t >= t_max && sum(N) > 0.                 # If any active individuals remain
        for t in 1:n_types
            for i in active[t]
                linelist[i, :t_death] = t_max    # Truncate their time of death
                if rand() <  ρ₀[t]                  # If they are sampled
                    linelist[i, :t_sam] = t_max  # Update their sampling time
                end
            end
        end
    end


    return Outbreak(
                params,
                vcat(t_out', convert(Array, N_out)),
                linelist
    )
end
