

function pop_random!(v::Vector)
    isempty(v) && throw(ArgumentError("Cannot pop from an empty vector."))
    i = rand(1:length(v))
    v[i], v[end] = v[end], v[i]   # swap with last element
    return pop!(v)
end


mutable struct Offspring
    child_ids::Vector{Int}
    infection_times::Vector{Float64}
    sampling_times::Vector{Float64}
end


Offspring() = Offspring(Vector{Int}(), Vector{Float64}(), Vector{Float64}())


mutable struct TransmissionChain
    infectors::Vector{Int}
    infection_times::Vector{Float64}
    sampling_times::Vector{Float64}
    sampling_ids::Vector{Int}
end


TransmissionChain(seed::Int) = TransmissionChain(
    fill(0, seed),
    fill(0.0, seed),
    fill(NaN, seed),
    Vector{Int}()
)


function Base.show(io::IO, chain::TransmissionChain)
    n_infected = length(chain.infectors)
    n_sampled = length(chain.sampling_ids)
    
    inf_times = chain.infection_times
    samp_times = filter(!isnan, chain.sampling_times)

    inf_range = isempty(inf_times) ? "N/A" : @sprintf("%.2f–%.2f", minimum(inf_times), maximum(inf_times))
    samp_range = isempty(samp_times) ? "N/A" : @sprintf("%.2f–%.2f", minimum(samp_times), maximum(samp_times))

    println(io, "TransmissionChain")
    println(io, "  Total infected       : $n_infected")
    println(io, "  Total sampled        : $n_sampled")
    println(io, "  Infection time range : $inf_range")
    println(io, "  Sampling time range  : $samp_range")
end


function infection!(chain::TransmissionChain, infector::Int, infection_time::Float64)
    push!(chain.infectors, infector)
    push!(chain.infection_times, infection_time)
    push!(chain.sampling_times, NaN)
end


function sampling!(chain::TransmissionChain, sampled::Int, sampling_time::Float64)
    push!(chain.sampling_ids, sampled)
    chain.sampling_times[sampled] = sampling_time
end


abstract type AbstractEpiModel end


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


function simulate_chain(model::MultiTypeBirthDeathModel;
                        I_init::Vector{Int}=[1, 0],
                        N_max::Int=10_000,
                        S_max::Int=100)

    I = I_init
    n_cumulative = sum(I_init)
    n_sampled = 0

    infected = [collect(1:I_init[type]) for type in eachindex(I_init)]

    chain = TransmissionChain(sum(I_init))

    # Pre-calculate event rates
    event_rates = Vector{Float64}(undef, 3)

    t = 0.0

    while !all(isempty.(infected)) && n_cumulative < N_max && n_sampled < S_max
        
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
            infector = sample(infected[parent_type])
            push!(infected[child_type], n_cumulative)
            infection!(chain, infector, t)
        elseif rand_number ≤ (event_rates[1] + event_rates[2]) / total_event_rate
            # Death event
            death_type = wsample(1:model.n_types, I)
            I[death_type] -= 1
            pop_random!(infected[death_type])
        else
            # Sampling event
            sampled_type = wsample(1:model.n_types, I)
            I[sampled_type] -= 1
            n_sampled += 1
            sampled = pop_random!(infected[sampled_type])
            sampling!(chain, sampled, t)
        end
    end
    return chain
end


function prune_tree(chain::TransmissionChain)
    return prune_tree(chain.sampling_ids, chain.infectors, chain.infection_times, chain.sampling_times)
end


function prune_tree(sampling_ids::Vector{Int}, 
                    infectors::Vector{Int}, 
                    infection_times::Vector{Float64}, 
                    sampling_times::Vector{Float64})
    sampled_ancestors = fill(false, length(infectors))
    transmission_tree = Dict{Int, Offspring}()
    for s in sampling_ids
        a = s
        while a > 0 && !sampled_ancestors[a]
            infector = infectors[a]
            offspring = get!(transmission_tree, infector, Offspring())
            push!(offspring.child_ids, a)
            push!(offspring.infection_times, infection_times[a])
            push!(offspring.sampling_times, sampling_times[a])
            sampled_ancestors[a] = true
            a = infector
        end
    end
    return transmission_tree
end