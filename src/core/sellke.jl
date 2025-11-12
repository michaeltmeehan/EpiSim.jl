function sellke(I₀::Int, 
                S₀::Int, 
                β::Float64, 
                γ::Float64)

    @assert I₀ ≥ 0 "Number of initial infected I₀ must be at least 0"
    @assert S₀ ≥ 0 "Number of susceptibles S₀ must be at least 0"
    @assert β > 0.0 "Transmission rate β must be positive"
    @assert γ > 0.0 "Recovery rate γ must be positive"

    # Initialize number of infected, susceptible and recovered
    I = I₀; S = S₀; R = 0; N = I₀ + S₀

    # Rescale transmission rate
    β /= N

    Λ = 0.0     # Cumulative infection pressure
    t = 0.0     # Current time

    # Draw resistances for susceptibles (Exp(1))
    resistances = BinaryMinHeap{Float64}(randexp(S₀))

    # Draw recovery times for initial infecteds
    recoveries = BinaryMinHeap{Float64}(t .+ randexp(I₀) ./ γ)

    while I > 0 && S > 0

        # Work out whether infection or recovery occurs next
        if top(resistances) ≤ Λ + (β * I) * (top(recoveries)[1] - t) # Infection event

            # Update time
            t += (top(resistances) - Λ) / ((β * I))

            # Update number of infected and susceptibles
            I += 1; S -= 1

            # Update cumulative infection pressure (and increment heap)
            Λ = pop!(resistances)

            # Add to recoveries heap
            push!(recoveries, t + randexp() / γ)

        else    # Recovery event
            
            # Update cumulative infection pressure
            Λ += (β * I) * (top(recoveries)[1] - t)

            # Update number of infected and recovered
            I -= 1; R += 1

            # Update time (and increment heap)
            t = pop!(recoveries)
        end
    end
    return (S=S, I=I, R=R)
end



draw_recovery_time(γ::Float64) = randexp() / γ
draw_recovery_time(n::Int, γ::Float64) = randexp(n) / γ
draw_recovery_time(d::Distribution{Univariate,Continuous}) = rand(d)
draw_recovery_time(n::Int, d::Distribution{Univariate,Continuous}) = rand(d, n)

draw_transmission_rate(β::Float64) = β
draw_transmission_rate(n::Int, β::Float64) = fill(β, n)
draw_transmission_rate(d::Distribution{Univariate,Continuous}) = rand(d)
draw_transmission_rate(n::Int, d::Distribution{Univariate,Continuous}) = rand(d, n)

function sellke(I₀::Int, 
                S₀::Int, 
                transmission_rate::Union{Float64, Distribution{Univariate,Continuous}}, 
                infectious_period::Union{Float64, Distribution{Univariate,Continuous}})

    @assert I₀ ≥ 0 "Number of initial infected I₀ must be at least 0"
    @assert S₀ ≥ 0 "Number of susceptibles S₀ must be at least 0"

    # Initialize number of infected, susceptible and recovered
    I = I₀; S = S₀; R = 0; N = I₀ + S₀

    Λ = 0.0     # Cumulative infection pressure
    t = 0.0     # Current time

    # Draw resistances for susceptibles (Exp(1))
    resistances = BinaryMinHeap{Float64}(randexp(S₀))

    # Draw transmission rate and recovery times for initial infecteds
    recovery_times = t .+ draw_recovery_time(I₀, infectious_period)
    transmission_rates = draw_transmission_rate(I₀, transmission_rate)
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
            τ = draw_recovery_time(infectious_period)
            β = draw_transmission_rate(transmission_rate)

            # Add to recoveries heap
            push!(recoveries, (t + τ, β))

            # Update slope
            slope += (β / N)

        else    # Recovery event

            # Update number of infected and recovered
            I -= 1; R += 1

            # Update cumulative infection pressure
            Λ += slope * (top(recoveries)[1] - t)

            # Update time (and increment heap)
            t, β = pop!(recoveries)

            # Decrease slope
            slope -= (β / N)
        end
    end
    return (S=S, I=I, R=R)
end