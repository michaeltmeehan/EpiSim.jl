mutable struct Offspring
    child_ids::Vector{Int}
    infection_times::Vector{Float64}
    sampling_times::Vector{Float64}
end


Offspring() = Offspring(Vector{Int}(), Vector{Float64}(), Vector{Float64}())


struct SampledTransmissionTree
    tree::Dict{Int, Offspring}
    root::Int
end


function extract_sampled_tree(chain::TransmissionChain)
    return extract_sampled_tree(chain.sampling_ids, chain.infectors, chain.infection_times, chain.sampling_times)
end


function extract_sampled_tree(sampling_ids::Vector{Int}, 
                              infectors::Vector{Int}, 
                              infection_times::Vector{Float64}, 
                              sampling_times::Vector{Float64})
    sampled_ancestors = fill(false, length(infectors))
    tree = Dict{Int, Offspring}()
    for s in sampling_ids
        a = s
        while a > 0 && !sampled_ancestors[a]
            infector = infectors[a]
            offspring = get!(tree, infector, Offspring())
            push!(offspring.child_ids, a)
            push!(offspring.infection_times, infection_times[a])
            push!(offspring.sampling_times, sampling_times[a])
            sampled_ancestors[a] = true
            a = infector
        end
    end
    return tree
end