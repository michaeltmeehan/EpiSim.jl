using DataFrames
using EpiSim.Models


abstract type AbstractOutbreak end


struct Outbreak{T<:AbstractEpiModel} <: AbstractOutbreak
    model::T
    events::Vector{AbstractEpiEvent}
    linelist::DataFrame
    trajectory::Matrix{Float64}
end


# function enrich(out::Outbreak)
#     ismissing(out.linelist) || ismissing(out.trajectory) || return out
#     ll, traj = process_events(out.events)
#     return Outbreak(out.model, out.events; linelist=ll, trajectory=traj)
# end

