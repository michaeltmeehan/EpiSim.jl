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


function simulate_outbreak(rng::AbstractRNG,
                           model::BirthDeathModel;
                           N_max::Int=10_000, 
                           t_max::Float64=100.0, 
                           S_max::Int=100, 
                           I_init::Int=1)

    I = I_init
    n_cumulative = I_init
    n_sampled = 0

    currently_infected = collect(1:I_init)

    events = Vector{AbstractEpiEvent}()

    for seed_host in 1:I_init
        push!(events, Seed(seed_host, 0.0))
    end

    # Pre-calculate event rates
    event_rates = Vector{Float64}(undef, 3)
    update_event_rates!(event_rates, model)
    total_event_rate = sum(event_rates)
    
    t = 0.0

    while !isempty(currently_infected) && n_cumulative < N_max && n_sampled < S_max
        rand_number = rand(rng)
        t -= log(rand_number) / (total_event_rate * I)
        t > t_max && break

        if rand_number ≤ event_rates[1] / total_event_rate
            # Birth event
            I += 1
            n_cumulative += 1
            infectee = n_cumulative         # Label infected individuals sequentially
            infector = sample(rng, currently_infected)
            push!(currently_infected, infectee)
            transmission!(events, infector, infectee, t)
        elseif rand_number ≤ (event_rates[1] + event_rates[2]) / total_event_rate
            # Death event
            I -= 1
            recovered = pop_random!(rng, currently_infected)
            recovery!(events, recovered, t)
        else
            # Sampling event
            I -= 1
            n_sampled += 1
            sampled = pop_random!(rng, currently_infected)
            sampling!(events, sampled, t)
        end
    end
    return events
end


function simulate_outbreak(model::BirthDeathModel;
                           N_max::Int=10_000, 
                           t_max::Float64=100.0, 
                           S_max::Int=100, 
                           I_init::Int=1)
    return simulate_outbreak(Random.GLOBAL_RNG, model; N_max=N_max, t_max=t_max, S_max=S_max, I_init=I_init)
end