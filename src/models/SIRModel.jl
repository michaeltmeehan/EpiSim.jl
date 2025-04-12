struct SIRModel <: AbstractEpiModel
    transmission_rate::Float64
    recovery_rate::Float64
    sampling_rate::Float64
    N::Int
end


@inline function update_event_rates!(event_rates::Vector{Float64}, 
                                     model::SIRModel, 
                                     S::Int, 
                                     I::Int)
    β = model.transmission_rate
    γ = model.recovery_rate
    ψ = model.sampling_rate
    N = model.N

    event_rates[1] = β * I * S / N  # Infection rate
    event_rates[2] = γ * I        # Recovery rate
    event_rates[3] = ψ * I        # Sampling rate
end


function simulate_outbreak(model::SIRModel; 
                           S_init::Int=9999, 
                           I_init::Int=1,
                           S_max::Int=100)

    # Initialize population parameters
    S = S_init
    I = I_init

    # Initialize the simulation parameters
    n_cumulative = I_init
    n_sampled = 0
    currently_infected = collect(1:I_init)

    events = Vector{AbstractEpiEvent}()
    event_rates = Vector{Float64}(undef, 3)
    
    t = 0.0
    
    # Add initial infections as seeds
    for i in 1:I_init
        push!(events, Seed(i, 0.0))
    end

    while !isempty(currently_infected) && n_sampled < S_max

        update_event_rates!(event_rates, model, S, I)
        total_event_rate = sum(event_rates)

        rand_number = rand()
        t -= log(rand_number) / total_event_rate

        if rand_number ≤ event_rates[1] / total_event_rate
            # Infection event
            S -= 1
            I += 1
            n_cumulative += 1
            infectee = n_cumulative
            
            # Get a random infector from the infected pool
            infector = sample(currently_infected)
            transmission!(events, infector, infectee, t)
            push!(currently_infected, infectee)
        elseif rand_number ≤ (event_rates[1] + event_rates[2]) / total_event_rate
            # Recovery event
            I -= 1
            recovered = pop_random!(currently_infected)
            recovery!(events, recovered, t)
        else
            # Sampling event
            I -= 1
            sampled = pop_random!(currently_infected)
            sampling!(events, sampled, t)
            n_sampled += 1
        end
    end
    return events
end