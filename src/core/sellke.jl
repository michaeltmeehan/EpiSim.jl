

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

    while I > 0 && S > 0

        # Work out whether infection or recovery occurs next
        if top(resistances) ≤ Λ + slope * (top(recoveries)[1] - t) # Infection event

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