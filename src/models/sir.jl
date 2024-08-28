function initialize_linelist(params::SIRParameters)::DataFrame
    N = params.I₀
    return DataFrame(
        event = fill(1, N),
        child_id = collect(1:N),
        parent_id = fill(0, N),
        t_infection = fill(0., N),
        t_recovery = fill(Inf, N),
        t_sam = fill(-1., N))
end


function check_params(params::SIRParameters)
    @assert params.β > 0.0 "Transmission rate must be positive."
    @assert params.γ > 0.0 "Recovery rate must be positive."
    @assert params.ψ >= 0.0 "Sampling rate must be non-negative."
    @assert params.N₀ > 0 "Initial population size must be positive."
    @assert 0 < params.I₀ ≤ params.N₀ "Initial number of infected individuals must be positive and less than N₀."
    @assert params.t_max > 0.0 "Maximum simulation time must be positive."
end


# Add docstring for simulate_outbreak
"""
    simulate_outbreak(params::SIRParameters; S_max::Union{Float64, Int64}=Inf)

Simulates an outbreak using the SIR model.

# Arguments
- `params::SIRParameters`: The parameters for the SIR model.
- `S_max::Union{Float64, Int64}`: The maximum number of sampled individuals.

# Returns
- `Outbreak`: The resulting outbreak object containing the linelist and trajectory.

# Example
```julia
params = SIRParameters(β=2.0, γ=1.0, ψ=0.5, N₀=1000, I₀=1, t_max=100.0)
outbreak = simulate_outbreak(params, S_max=100)
```
"""
function simulate_outbreak(params::SIRParameters; S_max::Union{Float64, Int64}=Inf)::Outbreak
    
        @unpack N₀, β, γ, ψ, r, I₀, t_max = params
        check_params(params)
        
        # Make sure some stopping criteria is specified
        S_max = isinf(min(S_max, t_max)) ? 1e3 : S_max + 0.
        
        # Initialize the linelist DataFrame to store event details for each individual
        linelist = initialize_linelist(params)
        
        # Convert population size to Float64
        N = N₀ + 0.

        # Initialize susceptible individuals and recovered individuals
        Susceptible = N - I₀
        Infected = I₀ + 0.
        Recovered = 0.
        
        # Initialize trajectory
        t = 0.
        t_out = [t]
        S_out = [Susceptible]
        I_out = [Infected]
        R_out = [Recovered]
        
        # Create a cohort of active individuals
        active = collect(1:I₀)
        cum_inc = I₀

        # Track cumulative sampled individuals
        Sampled = 0.
        
        # Main simulation loop
        while (Susceptible > 0) && (t < t_max) && (Sampled < S_max)

            # Calculate event rates
            total_event_rate = (β * Susceptible / N + γ + ψ) * Infected
            infection_threshold = (β * Susceptible * Infected / N) / total_event_rate
            recovery_threshold = (β * Susceptible * Infected / N + γ * Infected) / total_event_rate
            
            rnd = rand()
            t -= log(rnd) / (total_event_rate)
            t > t_max && break
            
            if rnd < infection_threshold        # Infection event
                cum_inc += 1                                                # Increment cumulative incidence
                parent = sample(active)                                     # Sample random parent from current cohort
                push!(linelist, (1, cum_inc, parent, t, Inf, -1.))          # Add event to linelist
                push!(active, cum_inc)                                      # Update list of active individuals
                Susceptible -= 1                                            # Update susceptible population size
                Infected += 1.                                              # Update infected population size
            elseif rnd < recovery_threshold     # Recovery event
                recovered_idx = sample(1:length(active))                    # Sample index of recovered individual
                recovered = active[recovered_idx]                           # Retrieve id of recovered individual
                linelist[recovered, :t_recovery] = t                        # Set the time of recovery for recovered individual
                deleteat!(active, recovered_idx)                            # Update list of active individuals
                Infected -= 1.                                              # Update infected population size
                Recovered += 1.                                             # Update recovered population size
            else                                                           # Sampling event
                sampled_idx = sample(1:length(active))                      # Sample index of sampled individual
                sampled = active[sampled_idx]                               # Retrieve id of sampled individual
                if linelist[sampled, :t_sam] < 0.                           # Increment cumulative sampled count (if not already sampled)
                    Sampled += 1.
                end
                linelist[sampled, :t_sam] = t                               # Update time of sampling for sampled individual
                if rand() < r                                               # Removal upon sampling
                    linelist[sampled, :t_recovery] = t                      # Update time of removal
                    deleteat!(active, sampled_idx)                          # Update list of active individuals
                    Infected -= 1.                                          # Update total population size
                    Recovered += 1.                                         # Update recovered population size
                end
        end
        push!(t_out, t)                                                 # Update trajectory
        push!(S_out, Susceptible)
        push!(I_out, Infected)
        push!(R_out, Recovered)

    end

    return Outbreak(
        params,
        vcat(t_out', S_out', I_out', R_out'),
        linelist
        )
end