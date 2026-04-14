function gillespie(S₀::Real,
                   E₀::Int,
                   I₀::Int,
                   β::Float64,
                   α::Float64,
                   γ::Float64,
                   ψ::Float64,
                   r::Float64)

    # Initialize state
    S, E, I, R = S₀, E₀, I₀, 0

    # Initialize time
    t = 0.0

    # Initialize total population size
    N = S + E + I + R

    # Initialize population of exposeds
    exposeds = collect(1:E₀)

    # Initialize population of infectives
    infectives = collect((E₀+1):(E₀+I₀))

    # Initialize event log
    times = Float64[]
    hosts = Int[]
    infectors = Int[]
    kinds = EventKind[]
    for id in 1:(E₀ + I₀)
        push!(times, 0.0)
        push!(hosts, id)
        push!(infectors, 0)
        push!(kinds, EK_Seeding)
    end

    # Initialize rates
    λ = Vector{Float64}(undef, 4)    # Transmission, Activation, Recovery, Sampling

    # Initial rates
    λ[1] = isfinite(S) ? β * S * I / N : β * I      # Transmission
    λ[2] = iszero(E) ? 0.0 : α * E              # Activation
    λ[3] = γ * I              # Recovery
    λ[4] = ψ * I              # Sampling

    while sum(λ) > 0.

        # Total rate
        λ_total = sum(λ)

        # Time to next event
        t += randexp() / λ_total

        # Determine which event occurs
        r_event = rand() * λ_total
        cumulative_rate = 0.0
        event_type = 0

        for i in eachindex(λ)
            cumulative_rate += λ[i]
            if r_event ≤ cumulative_rate
                event_type = i
                break
            end
        end

        infector = 0

        # Update state based on event type
        if event_type == 1  # Transmission
            S -= 1; E += 1
            new_exposed = E + I + R
            kind = EK_Transmission
            host = new_exposed
            infector = rand(infectives)
            push!(exposeds, new_exposed)  # Add new exposed to population

        elseif event_type == 2  # Activation
            E -= 1; I += 1
            new_infected = popr!(exposeds)
            kind = EK_Activation
            host = new_infected
            push!(infectives, new_infected)  # Add new infected to population

        elseif event_type == 3  # Recovery
            I -= 1; R += 1
            new_recovered = popr!(infectives)
            kind = EK_Recovery
            host = new_recovered

        elseif event_type == 4  # Sampling
            new_sampled = popr!(infectives)

            # Determine if sampled individual recovers/removes
            if rand() < r
                I -= 1; R += 1
                kind = EK_SerialSampling
                host = new_sampled
            else
                # If not recovered, put back into infectives
                kind = EK_FossilizedSampling
                host = new_sampled
                push!(infectives, new_sampled)
            end

        end

        push!(times, t)
        push!(hosts, host)
        push!(infectors, infector)
        push!(kinds, kind)

        # Update rates
        λ[1] = isfinite(S) ? β * S * I / N : β * I      # Transmission
        λ[2] = iszero(E) ? 0.0 : α * E              # Activation
        λ[3] = γ * I              # Recovery
        λ[4] = ψ * I              # Sampling

    end
    el = EventLog(times, hosts, infectors, kinds)
    population_size = isfinite(N) && isinteger(N) ? Int(N) : nothing
    validate_event_log(el; population_size=population_size)
    return el
end
