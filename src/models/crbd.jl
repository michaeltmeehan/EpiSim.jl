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
    simulate(N₀::Int64, λ::Float64, μ::Float64, ψ::Float64, ρ₀::Float64, r::Float64;
             N_max::Union{Float64, Int64}=Inf, S_max::Union{Int64, Float64}=Inf, t_max::Union{Int64, Float64}=Inf)::Outbreak

Simulates an outbreak using the specified parameters for the CRBD model.

## Arguments
- `N₀::Int64`: Initial population size.
- `λ::Float64`: Birth rate.
- `μ::Float64`: Death rate.
- `ψ::Float64`: Extinct/ancestral sampling rate.
- `ρ₀::Float64`: Extant sampling rate.
- `r::Float64`: Removal probability (upon sampling).
- `N_max::Union{Float64, Int64}`: The maximum number of individuals to simulate (default is `Inf`).
- `S_max::Union{Int64, Float64}`: The maximum number of sampled individuals (default is `Inf`).
- `t_max::Union{Int64, Float64}`: Maximum simulation time (default is `Inf`).

## Returns
- `Outbreak`: The resulting outbreak object containing the linelist and trajectory.

This function creates a `CRBDParameters` instance using the provided arguments and calls the `simulate` function that takes `CRBDParameters` as input.
"""
function simulate(N₀::Int64, 
                  λ::Float64, 
                  μ::Float64, 
                  ψ::Float64,
                  ρ₀::Float64,
                  r::Float64;
                  N_max::Union{Float64, Int64} = Inf,
                  S_max::Union{Int64, Float64} = Inf,
                  t_max::Union{Int64, Float64} = Inf)::Outbreak

    params = CRBDParameters(N₀, λ, μ, ψ, ρ₀, r, t_max)

    return simulate(params, N_max=N_max, S_max=S_max)
end


"""
    simulate(params::CRBDParameters; N_max::Union{Float64, Int64}=Inf, S_max::Union{Float64, Int64}=Inf)

Simulates an outbreak using the CRBD model.

# Arguments
- `params::CRBDParameters`: The parameters for the CRBD model.
- `N_max::Union{Float64, Int64}`: The maximum number of individuals to simulate.
- `S_max::Union{Float64, Int64}`: The maximum number of sampled individuals.

# Returns
- `Outbreak`: The resulting outbreak object containing the linelist and trajectory.
"""
function simulate(params::CRBDParameters; N_max::Union{Float64, Int64}=Inf, S_max::Union{Float64, Int64}=Inf)::Outbreak

    @unpack N₀, λ, μ, ψ, ρ₀, r, t_max = params
    check_parameters(params)
    
    # Make sure some stopping criteria is specified
    N_max = isinf(min(N_max, S_max, t_max)) ? 1e3 : N_max + 0.
    
    # Initialize the linelist DataFrame to store event details for each individual
    linelist = initialize_linelist(N₀)
    
    # Convert population size to Float64
    N = N₀ + 0.
    
    # Initialize trajectory
    t = 0.
    t_out = [t]
    N_out = [N]
    
    # Pre-calculate event rates
    total_event_rate = λ + μ + ψ
    birth_threshold = λ / total_event_rate
    death_threshold = (λ + μ) / total_event_rate
    
    # Create a cohort of active individuals
    active = collect(1:N₀)
    cum_inc = N₀
    
    # Track cumulative sampled individuals
    S = 0.
    
    # Main simulation loop
    while (N > 0.) && (cum_inc < N_max) && (S < S_max)
        rnd = rand()
        t -= log(rnd) / (total_event_rate * N)
        t > t_max && break
        if rnd ≤ birth_threshold                                       # Birth event
            cum_inc += 1                                                # Increment cumulative incidence
            parent = sample(active)                                     # Sample random parent from current cohort
            push!(linelist, (1, cum_inc, 1, parent, 1, t, Inf, -1.))    # Add event to linelist
            push!(active, cum_inc)                                      # Update list of active individuals
            N += 1.                                                     # Update total population size
        elseif rnd ≤ death_threshold                                   # Death event
            deceased_idx = sample(1:length(active))                     # Sample index of deceased individual
            deceased = active[deceased_idx]                             # Retrieve id of deceased individual
            linelist[deceased, :t_death] = t                            # Set the time of death for deceased individual
            deleteat!(active, deceased_idx)                             # Update list of active individuals
            N -= 1.                                                     # Update total population size
        else                                                           # Sampling event
            sampled_idx = sample(1:length(active))                      # Sample index of sampled individual
            sampled = active[sampled_idx]                               # Retrieve id of sampled individual
            if linelist[sampled, :t_sam] < 0.                           # Increment cumulative sampled count (if not already sampled)
                S += 1.
            end
            linelist[sampled, :t_sam] = t                               # Update time of sampling for sampled individual
            if rand() < r                                               # Removal upon sampling
                linelist[sampled, :t_death] = t                         # Update time of death / removal time
                deleteat!(active, sampled_idx)                          # Update list of active individuals
                N -= 1.                                                 # Update total population size
            end
        end
        push!(t_out, t)                                                 # Update trajectory
        push!(N_out, N)
    end

    # Present day sampling
    if t ≥ t_max && N > 0.
        for i in active                                        
            linelist[i, :t_death] = t_max                               # Truncate survival time for remaining active individuals
            if rand() <  ρ₀                                             # Sample with probability ρ₀
                linelist[i, :t_sam] = t_max
            end
        end
    end

    return Outbreak(
                params,
                vcat(t_out', N_out'),
                linelist
    )
end