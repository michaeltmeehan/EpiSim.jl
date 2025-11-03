

function sellke(m::Int, n::Int, β::Float64, γ::Float64)

    @assert m ≥ 1 "Number of initial infected m must be at least 1"
    @assert n ≥ 0 "Number of susceptibles n must be at least 0"

    # Initialize number of infected, susceptible and recovered
    I = m; S = n; R = 0; N = m + n

    # Draw resistances for susceptibles (Exp(1))
    resistances = BinaryMinHeap{Float64}()
    for _ in 1:n
        push!(resistances, randexp())
    end

    next_resistance = isempty(resistances) ? Inf : pop!(resistances)

    Λ = 0.0     # Cumulative infection pressure
    t = 0.0     # Current time

    # Draw recovery times for initial infected
    recoveries = BinaryMinHeap{Float64}()
    for _ in 1:m
        push!(recoveries, t + randexp() / γ)
    end
    next_recovery = top(recoveries)

    while I > 0

        # Work out time to next recovery
        Δt = next_recovery - t

        # Check out if next infection occurs before next recovery
        if next_resistance ≤ Λ + (β / N) * I * Δt   # Infection event

            # Calculate time to infection
            δt = (next_resistance - Λ) / ((β / N) * I)

            # Update time
            t += δt

            # Update cumulative infection pressure
            Λ = next_resistance

            # Update number of infected and susceptibles
            I += 1; S -= 1

            # Get next resistance
            next_resistance = isempty(resistances) ? Inf : pop!(resistances)

            # Draw recovery time for newly infected
            push!(recoveries, t + randexp() / γ)
            next_recovery = top(recoveries)

        else    # Recovery event
            # Update time
            t = next_recovery

            # Update cumulative infection pressure
            Λ += (β / N) * I * Δt

            # Update number of infected and recovered
            I -= 1; R += 1

            # Remove next recovery from heap
            pop!(recoveries)

            # Get time of next recovery
            next_recovery = isempty(recoveries) ? Inf : top(recoveries)
        end
    end
    return (S=S, I=I, R=R)
end



draw_recovery_time(γ::Float64) = randexp() / γ
draw_recovery_time(d::Distribution{Univariate,Continuous}) = rand(d)

draw_transmission_rate(β::Float64) = β
draw_transmission_rate(d::Distribution{Univariate,Continuous}) = rand(d)


function sellke(m::Int, n::Int, transmission_rate::Union{Float64, Distribution{Univariate,Continuous}}, infectious_period::Union{Float64, Distribution{Univariate,Continuous}})

    @assert m ≥ 1 "Number of initial infected m must be at least 1"
    @assert n ≥ 0 "Number of susceptibles n must be at least 0"

    # Initialize number of infected, susceptible and recovered
    I = m; S = n; R = 0; N = m + n

    # Draw resistances for susceptibles (Exp(1))
    resistances = BinaryMinHeap{Float64}()
    for _ in 1:n
        push!(resistances, randexp())
    end

    next_resistance = isempty(resistances) ? Inf : pop!(resistances)

    Λ = 0.0     # Cumulative infection pressure
    t = 0.0     # Current time

    # Draw transmission rate and recovery times for initial infected
    recoveries = BinaryMinHeap{Tuple{Float64,Float64}}()
    slope = 0.0
    for _ in 1:m
        τ = draw_recovery_time(infectious_period)
        β = draw_transmission_rate(transmission_rate)
        push!(recoveries, (t + τ, β))
        slope += (β / N)
    end
    next_recovery = top(recoveries)

    while I > 0

        # Work out time to next recovery
        Δt = next_recovery[1] - t

        # Check out if next infection occurs before next recovery
        if next_resistance ≤ Λ + slope * Δt   # Infection event

            # Calculate time to infection
            δt = (next_resistance - Λ) / slope

            # Update time
            t += δt

            # Update cumulative infection pressure
            Λ = next_resistance

            # Update number of infected and susceptibles
            I += 1; S -= 1

            # Get next resistance
            next_resistance = isempty(resistances) ? Inf : pop!(resistances)

            # Draw recovery time and transmission rate for newly infected
            τ = draw_recovery_time(infectious_period)
            β = draw_transmission_rate(transmission_rate)
            push!(recoveries, (t + τ, β))

            # Update slope
            slope += (β / N)

            # Get time of next recovery
            next_recovery = top(recoveries)

        else    # Recovery event
            # Update cumulative infection pressure
            Λ += slope * Δt

            t, β = next_recovery

            # Decrease slope
            slope -= (β / N)

            # Update number of infected and recovered
            I -= 1; R += 1

            # Remove next recovery from heap
            pop!(recoveries)

            # Get time of next recovery
            next_recovery = isempty(recoveries) ? Inf : top(recoveries)
        end
    end
    return (S=S, I=I, R=R)
end