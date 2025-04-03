struct SIRModel <: AbstractEpiModel
    transmission_rate::Float64
    recovery_rate::Float64
    sampling_rate::Float64
    N::Int
end


@inline function update_event_rates!(event_rates::Vector{Float64}, model::SIRModel, S::Int, I::Int)
    β = model.transmission_rate
    γ = model.recovery_rate
    ψ = model.sampling_rate
    N = model.N

    event_rates[1] = β * I * S / N  # Infection rate
    event_rates[2] = γ * I        # Recovery rate
    event_rates[3] = ψ * I        # Sampling rate
end


function simulate_chain(model::SIRModel; 
                        S_init::Int=9999, 
                        I_init::Int=1,
                        S_max::Int=100)

    # Initialize population parameters
    S = S_init
    I = I_init

    # Initialize the simulation parameters
    n_cumulative = I_init
    n_sampled = 0
    infected = Vector{Int}(undef, I_init)
    infected[1:I_init] .= 1:I_init

    chain = TransmissionChain(I_init)
    event_rates = Vector{Float64}(undef, 3)
    
    t = 0.0

    while !isempty(infected) && n_sampled < S_max

        update_event_rates!(event_rates, model, S, I)
        total_event_rate = sum(event_rates)

        rand_number = rand()
        t -= log(rand_number) / total_event_rate

        if rand_number ≤ event_rates[1] / total_event_rate
            # Infection event
            S -= 1
            I += 1
            n_cumulative += 1
            infection!(chain, sample(infected), t)
            push!(infected, n_cumulative)
        elseif rand_number ≤ (event_rates[1] + event_rates[2]) / total_event_rate
            # Recovery event
            I -= 1
            pop_random!(infected)
        else
            # Sampling event
            I -= 1
            n_sampled += 1
            sampling!(chain, pop_random!(infected), t)
        end
    end
    return chain
end