

function sellke(m::Int, n::Int, β::Float64, d::Distribution{Univariate,Continuous})

    @assert m ≥ 1 "Number of initial infected m must be at least 1"
    @assert n ≥ 0 "Number of susceptibles n must be at least 0"

    # Initialize number of infected, susceptible and recovered
    I = m; S = n; R = 0

    # Draw resistances for susceptibles (Exp(1))
    resistances = BinaryMinHeap{Float64}()
    for _ in 1:n
        push!(resistances, randexp())
    end

    next_resistance = isempty(resistances) ? Inf : pop!(resistances)

    Λ = 0.0     # Cumulative infection pressure
    t = 0.0     # Current time

    # Draw recovery times for initial infected
    # recoveries = BinaryMinHeap{Tuple{Float64,Int}}()
    recoveries = BinaryMinHeap{Float64}()
    for _ in 1:m
        push!(recoveries, t + rand(d))
    end
    next_recovery = top(recoveries)

    while I > 0

        # Work out time to next recovery
        Δt = next_recovery - t

        # Check out if next infection occurs before next recovery
        if next_resistance ≤ Λ + β * I * Δt

            # Calculate time to infection
            δt = (next_resistance - Λ) / (β * I)

            # Update time
            t += δt

            # Update cumulative infection pressure
            Λ = next_resistance

            # Update number of infected and susceptibles
            I += 1; S -= 1

            # Get next resistance
            next_resistance = isempty(resistances) ? Inf : pop!(resistances)

            # Draw recovery time for newly infected
            push!(recoveries, t + rand(d))
            next_recovery = top(recoveries)

        else
            t = next_recovery
            Λ += β * I * Δt
            I -= 1
            R += 1
            pop!(recoveries)
            next_recovery = isempty(recoveries) ? Inf : top(recoveries)
        end

    end
    return (S=S, I=I, R=R)
end