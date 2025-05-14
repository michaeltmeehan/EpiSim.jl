
# Default function for all models - can be overloaded if required
function get_prevalence_timeseries(df::DataFrame)
    if :I ∉ propertynames(df)
        throw(ArgumentError("Prevalence column `I` not found in DataFrame"))
    end
    (; t = df.t, prevalence = df.I)
end

function get_prevalence_timeseries(out::AbstractOutbreak)
    return get_prevalence_timeseries(out.state_log)
end


# An example overload for a specific model
# function get_prevalence_timeseries(out::Outbreak{<:SEIRModel})
#     df = out.state_log
#     (; t = df.t, prevalence = df.E + df.I)
# end


function get_prevalence(df::DataFrame, t::Real)
    idx = searchsortedlast(df.t, t)
    if isnothing(idx)
        throw(ArgumentError("No state recorded at or before t = $t"))
    end
    return get_prevalence_timeseries(df).prevalence[idx]
end


function get_prevalence(out::AbstractOutbreak, t::Real)
    return get_prevalence(out.state_log, t)
end


function get_prevalence(df::DataFrame, tvec::AbstractVector{<:Real})
    @assert issorted(tvec) "Query times must be sorted for get_prevalence(...)."

    ts_prev = get_prevalence_timeseries(df)
    ts = ts_prev.t
    prevalences = ts_prev.prevalence

    idxs = map(t -> searchsortedlast(ts, t), tvec)
    return [iszero(i) ? throw(ArgumentError("No state before t = $(tvec[j])")) : prevalences[i]
            for (j, i) in enumerate(idxs)]
end


function get_prevalence(out::AbstractOutbreak, t::AbstractVector{<:Real})
    return get_prevalence(out.state_log, t)
end


function get_prevalence(outbreaks::Vector{<:AbstractOutbreak}, t::Real)
    return [get_prevalence(out, t) for out in outbreaks]
end


function get_prevalence(outbreaks::Vector{<:AbstractOutbreak}, tvec::AbstractVector{<:Real})
    @assert issorted(tvec) "Query times must be sorted for get_prevalence(...)."

    n_outbreaks = length(outbreaks)
    n_t = length(tvec)
    prevalence_matrix = Matrix{Float64}(undef, n_t, n_outbreaks)

    for (j, out) in enumerate(outbreaks)
        prevalence_matrix[:, j] .= get_prevalence(out, tvec)
    end

    return prevalence_matrix  # rows = timepoints, cols = outbreak replicates
end


function get_prevalence(ens::AbstractEnsemble, tvec::AbstractVector{<:Real})
    return get_prevalence(ens.replicates, tvec)
end



function isextinct(df::DataFrame)::Bool
    if :I ∉ propertynames(df)
        throw(ArgumentError("Prevalence column `I` not found in DataFrame"))
    end
    return df.I[end] == 0.0
end


function isextinct(out::AbstractOutbreak)::Bool
    return isextinct(out.state_log)
end


function isextinct(outbreaks::Vector{<:AbstractOutbreak})::Vector{Bool}
    return [isextinct(out) for out in outbreaks]
end


function isextinct(ens::AbstractEnsemble)::Vector{Bool}
    return isextinct(ens.replicates)
end


function calc_extinction_probability(ens::AbstractEnsemble)::Float64
    n_extinct = count(isextinct(ens))
    return n_extinct / length(ens.replicates)
end


# TODO: Implement other processing functions
# state_log
# get_incidence_timeseries
# get_incidence
# get_cumulative_incidence
# get_removed
# get_peak_prevalence
# get_duration
# get_extinction_time
# get_trajectory
# get_R0_estimate

# event_log
# count_events
# get_sample_times
# get_transmission_times
# get_final_size

# Aggregate metrics
# was_detected
# reached_threshold(out, x)


