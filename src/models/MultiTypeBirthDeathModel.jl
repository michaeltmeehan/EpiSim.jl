struct MultiTypeBirthDeathModel <: AbstractEpiModel
    birth_rate::Matrix{Float64} # Birth rate matrix (n_types x n_types): element (i,j) is the rate individuals in state i give birth to individuals in state j
    death_rate::Vector{Float64}
    sampling_rate::Vector{Float64}
    n_types::Int
    Λ::Vector{Float64} # Total birth rate for each type
end


function MultiTypeBirthDeathModel(birth_rate::Matrix{Float64}, 
                                  death_rate::Vector{Float64}, 
                                  sampling_rate::Vector{Float64})::MultiTypeBirthDeathModel
    typeof(death_rate)
    n_types = size(birth_rate, 1)
    Λ = sum(birth_rate, dims=2)[:]  # Total birth rate for each type
    return MultiTypeBirthDeathModel(birth_rate, death_rate, sampling_rate, n_types, Λ)
end


function update_event_rates!(event_rates::Vector{Float64}, model::MultiTypeBirthDeathModel, I::Vector{Int})
    Λ = model.Λ
    μ = model.death_rate
    ψ = model.sampling_rate

    event_rates[1] = Λ ⋅ I  # Birth rate
    event_rates[2] = μ ⋅ I  # Death rate
    event_rates[3] = ψ ⋅ I  # Sampling rate
end


function simulate_outbreak(model::MultiTypeBirthDeathModel;
                        I_init::Vector{Int}=[1, 0],
                        N_max::Int=10_000,
                        S_max::Int=100)

    I = I_init
    n_cumulative = sum(I_init)
    n_sampled = 0

    currently_infected = [collect(1:I_init[type]) for type in eachindex(I_init)]

    events = Vector{AbstractEpiEvent}()

    # Pre-calculate event rates
    event_rates = Vector{Float64}(undef, 3)
    
    t = 0.0
    
    # Add initial infections as seeds
    for type in eachindex(I_init)
        for i in 1:I_init[type]
            push!(events, Seed(i, 0.0))
        end
    end

    while !all(isempty.(currently_infected)) && n_cumulative < N_max && n_sampled < S_max
        
        update_event_rates!(event_rates, model, I)
        total_event_rate = sum(event_rates)
        
        rand_number = rand()
        t -= log(rand_number) / total_event_rate

        if rand_number ≤ event_rates[1] / total_event_rate
            # Birth event
            birth_weights = model.Λ .* I
            parent_type = wsample(1:model.n_types, birth_weights)
            child_type = wsample(1:model.n_types, model.birth_rate[parent_type, :])
            I[child_type] += 1
            n_cumulative += 1
            infectee = n_cumulative         # Label infected individuals sequentially
            infector = sample(currently_infected[parent_type])
            push!(currently_infected[child_type], infectee)
            transmission!(events, infector, infectee, t)
        elseif rand_number ≤ (event_rates[1] + event_rates[2]) / total_event_rate
            # Death event
            death_type = wsample(1:model.n_types, I)
            I[death_type] -= 1
            recovered = pop_random!(currently_infected[death_type])
            recovery!(events, recovered, t)
        else
            # Sampling event
            sampled_type = wsample(1:model.n_types, I)
            I[sampled_type] -= 1
            n_sampled += 1
            sampled = pop_random!(currently_infected[sampled_type])
            sampling!(events, sampled, t)
        end
    end
    return events
end