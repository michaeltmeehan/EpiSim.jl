module Transmission

import Printf: @sprintf

export TransmissionChain, infection!, sampling!

mutable struct TransmissionChain
    infectors::Vector{Int}
    infection_times::Vector{Float64}
    sampling_times::Vector{Float64}
end


TransmissionChain(seed::Int) = TransmissionChain(
    fill(0, seed),
    fill(0.0, seed),
    fill(NaN, seed),
)


Base.getindex(chain::TransmissionChain, i::Int) = (chain.infectors[i], chain.infection_times[i])


function Base.show(io::IO, chain::TransmissionChain)
    n_infected = length(chain.infectors)
    
    inf_times = chain.infection_times
    samp_times = filter(!isnan, chain.sampling_times)
    n_sampled = length(samp_times)

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
    chain.sampling_times[sampled] = sampling_time
end


end