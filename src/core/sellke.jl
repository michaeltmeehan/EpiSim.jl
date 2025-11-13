

function sellke(I₀::Int, 
                S₀::Int, 
                β::Union{Float64, Distribution}, 
                τ::Union{Float64, Distribution})

    @assert I₀ ≥ 0 "Number of initial infected I₀ must be at least 0"
    @assert S₀ ≥ 0 "Number of susceptibles S₀ must be at least 0"

    # Initialize number of infected, susceptible and recovered
    I = I₀; S = S₀; R = 0; N = I₀ + S₀

    Λ = 0.0     # Cumulative infection pressure
    t = 0.0     # Current time

    # Draw resistances for susceptibles (Exp(1))
    resistances = BinaryMinHeap{Float64}(randexp(S₀))

    # Convert β and τ if necessary
    dβ = typeof(β) <: Distribution ? β : Dirac(β)
    dτ = typeof(τ) <: Distribution ? τ : Exponential(τ)

    # Draw transmission rate and recovery times for initial infecteds
    recovery_times = t .+ rand(dτ, I₀)
    transmission_rates = rand(dβ, I₀)
    recoveries = BinaryMinHeap{Tuple{Float64,Float64}}(collect(zip(recovery_times, transmission_rates)))
    slope = sum(transmission_rates) / N

    while I > 0

        # Work out whether infection or recovery occurs next
        if !isempty(resistances) && top(resistances) ≤ Λ + slope * (top(recoveries)[1] - t) # Infection event

            # Update number of infected and susceptibles
            I += 1; S -= 1

            # Update time
            t += (top(resistances) - Λ) / slope

            # Update cumulative infection pressure (and increment heap)
            Λ = pop!(resistances)

            # Draw recovery time and transmission rate for newly infected
            τ = rand(dτ)
            βᵢ = rand(dβ)

            # Add to recoveries heap
            push!(recoveries, (t + τ, β))

            # Update slope
            slope += (βᵢ / N)

        else    # Recovery event

            # Update number of infected and recovered
            I -= 1; R += 1

            # Update cumulative infection pressure
            Λ += slope * (top(recoveries)[1] - t)

            # Update time (and increment heap)
            t, βᵢ = pop!(recoveries)

            # Decrease slope
            slope -= (βᵢ / N)
        end
    end
    return (S=S, I=I, R=R)
end





function sellke(S₀::Int, 
                E₀::Int,
                I₀::Int, 
                β::Union{Float64, Distribution},    # Transmission rate
                τₑ::Union{Float64, Distribution},   # Incubation period
                τᵢ::Union{Float64, Distribution},   # Infectious period
                τₛ::Union{Float64, Distribution},    # Sampling time
                r::Float64                           # Probability of removal upon sampling
                )

    @assert I₀ ≥ 0 "Number of initial infected I₀ must be at least 0"
    @assert S₀ ≥ 0 "Number of susceptibles S₀ must be at least 0"

    # Initialize number of exposed, infected, susceptible and recovered
    E = E₀; I = I₀; S = S₀; R = 0; N = E₀ + I₀ + S₀

    Λ = 0.0     # Cumulative infection pressure
    t = 0.0     # Current time

    # Draw resistances for susceptibles (Exp(1))
    resistances = BinaryMinHeap{Float64}(randexp(S₀))

    # Convert parameters to distributions if necessary
    dβ = typeof(β) <: Distribution ? β : Dirac(β)
    dτₑ = typeof(τₑ) <: Distribution ? τₑ : Exponential(τₑ)
    dτᵢ = typeof(τᵢ) <: Distribution ? τᵢ : Exponential(τᵢ)
    dτₛ = typeof(τₛ) <: Distribution ? τₛ : Exponential(τₛ)

    @assert 0 < mean(dβ) < Inf "Transmission rate distribution must have a positive, finite mean"
    @assert 0 < mean(dτₑ) < Inf "Incubation time distribution must have a positive, finite mean"
    @assert 0 < mean(dτᵢ) < Inf "Infectious period distribution must have a positive, finite mean"
    @assert 0 < mean(dτₛ) < Inf "Sampling time distribution must have a positive, finite mean"

    # Draw transmission rates for all infected individuals
    transmission_rates = rand(dβ, E₀ + I₀)
    slope = sum(transmission_rates) / N

    # Draw activation, infectious and sampling periods for initial infecteds
    activation_periods = vcat(rand(dτₑ, E₀), fill(0.0, I₀))
    infectious_periods = rand(dτᵢ, E₀ + I₀)
    sampling_periods = rand(dτₛ, E₀ + I₀)

    # Convert periods to event times
    activation_times = t .+ activation_periods
    recovery_times = activation_times .+ infectious_periods
    sampling_times = activation_times .+ sampling_periods

    event_times = activation_times
    for i in 1:(E₀ + I₀)
        if recovery_times[i] < sampling_times[i]
            # Individual recovers before being sampled
            push!(event_times, recovery_times[i])
        elseif rand() < r
            # Individual is removed upon sampling
            push!(event_times, sampling_times[i])
        else
            # Individual is sampled but not removed
            push!(event_times, sampling_times[i])
            push!(event_times, recovery_times[i])
        end
    end

    # Create heaps for events
    activations = collect(zip(activation_times, 1:(E₀ + I₀), :activation))
    recoveries = collect(zip(recovery_times, transmission_rates, :recovery))
    samplings = collect(zip(sampling_times, transmission_rates, :sampling))
    events = vcat(activations, recoveries, samplings)
    event_heap = BinaryMinHeap{Tuple{Float64, Float64, Symbol}}(events)

    while (E > 0 || I > 0) && S > 0

        # Work out whether infection or temporal event occurs next
        if top(resistances) ≤ Λ + slope * (top(event_heap)[1] - t) # Infection event

            # Update number of infected and susceptibles
            I += 1; S -= 1

            # Update time
            t += (top(resistances) - Λ) / slope

            # Update cumulative infection pressure (and increment heap)
            Λ = pop!(resistances)

            # Draw times and transmission rate for newly infected
            βⱼ = rand(dβ)
            τₑⱼ = rand(dτₑ)
            τᵢⱼ = rand(dτᵢ)
            τₛⱼ = rand(dτₛ)

            # Add times events heap
            push!(event_heap, (t + τₑⱼ, βⱼ, :activation))   # Activation event
            push!(event_heap, (t + τₑⱼ + τᵢⱼ, βⱼ, :recovery)) # Recovery event
            push!(event_heap, (t + τₑⱼ + τₛⱼ, βⱼ, :sampling)) # Sampling event

            # Update slope
            slope += (βⱼ / N)

        elseif top(event_heap)[3] == :activation   # Activation event

            # Update number of exposed and infected
            E -= 1; I += 1

            # Update time


            # Update number of infected and recovered
            I -= 1; R += 1

            # Update cumulative infection pressure
            Λ += slope * (top(recoveries)[1] - t)

            # Update time (and increment heap)
            t, βᵢ = pop!(recoveries)

            # Decrease slope
            slope -= (βᵢ / N)
        end
    end
    return (S=S, I=I, R=R)
end