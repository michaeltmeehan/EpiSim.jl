struct BirthDeathModel <: AbstractEpiModel
    birth_rate::Float64
    death_rate::Float64
    sampling_rate::Float64
end


function update_event_rates!(event_rates::Vector{Float64}, model::BirthDeathModel)
    λ = model.birth_rate
    μ = model.death_rate
    ψ = model.sampling_rate

    event_rates[1] = λ  # Birth rate
    event_rates[2] = μ  # Death rate
    event_rates[3] = ψ  # Sampling rate
end


function simulate_chain(model::BirthDeathModel;
                        N_max::Int=10_000, 
                        t_max::Float64=100.0, 
                        S_max::Int=100, 
                        I_init::Int=1)

    I = I_init
    n_cumulative = I_init
    n_sampled = 0

    infected = collect(1:I_init)

    chain = TransmissionChain(I_init)

    # Pre-calculate event rates
    event_rates = Vector{Float64}(undef, 3)
    update_event_rates!(event_rates, model)
    total_event_rate = sum(event_rates)
    
    t = 0.0

    while !isempty(infected) && n_cumulative < N_max && n_sampled < S_max
        rand_number = rand()
        t -= log(rand_number) / (total_event_rate * I)
        t > t_max && break

        if rand_number ≤ event_rates[1] / total_event_rate
            # Birth event
            I += 1
            n_cumulative += 1
            infector = sample(infected)
            push!(infected, n_cumulative)
            infection!(chain, infector, t)
        elseif rand_number ≤ (event_rates[1] + event_rates[2]) / total_event_rate
            # Death event
            I -= 1
            pop_random!(infected)
        else
            # Sampling event
            I -= 1
            n_sampled += 1
            sampled = pop_random!(infected)
            sampling!(chain, sampled, t)
        end
    end
    return chain
end