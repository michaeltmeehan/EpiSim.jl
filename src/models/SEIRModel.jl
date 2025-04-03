struct SEIRModel <: AbstractEpiModel
    transmission_rate::Float64
    activation_rate::Float64
    recovery_rate::Float64
    sampling_rate::Float64
    N::Int
end


function update_event_rates!(event_rates::Vector{Float64}, model::SEIRModel, S::Int, E::Int, I::Int)
    β = model.transmission_rate
    ν = model.activation_rate
    γ = model.recovery_rate
    ψ = model.sampling_rate
    N = model.N

    event_rates[1] = β * I * S / N  # Infection rate
    event_rates[2] = ν * E        # Activation rate
    event_rates[3] = γ * I        # Recovery rate
    event_rates[4] = ψ * I        # Sampling rate
end


function simulate_chain(model::SEIRModel;
                        S_init::Int=9999, 
                        E_init::Int=0,
                        I_init::Int=1,
                        S_max::Int=100)

    # Initialize population parameters
    S = S_init
    E = E_init
    I = I_init
     
    # Initialize the simulation parameters
    n_cumulative = I_init
    n_sampled = 0
    exposed = Vector{Int}(undef, E_init)
    infected = Vector{Int}(undef, I_init)
    infected[1:I_init] .= 1:I_init

    chain = TransmissionChain(I_init)
    event_rates = Vector{Float64}(undef, 4)
    
    t = 0.0

    while !isempty(infected) && n_sampled < S_max

        update_event_rates!(event_rates, model, S, E, I)
        total_event_rate = sum(event_rates)

        rand_number = rand()
        t -= log(rand_number) / total_event_rate

        if rand_number ≤ event_rates[1] / total_event_rate
            # Infection event
            S -= 1
            E += 1
            n_cumulative += 1
            infection!(chain, sample(infected), t)
            push!(exposed, n_cumulative)
        elseif rand_number ≤ (event_rates[1] + event_rates[2]) / total_event_rate
            # Activation event
            E -= 1
            I += 1
            push!(infected, pop_random!(exposed))
        elseif rand_number ≤ (event_rates[1] + event_rates[2] + event_rates[3]) / total_event_rate
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