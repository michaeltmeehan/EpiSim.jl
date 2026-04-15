const _TRAJECTORY_COMPARTMENTS = (:S, :E, :I, :R)
const _HOST_SUMMARY_METRICS = (:transmissions_caused, :transmissions, :samples, :removals, :activations)


"""
    trajectory_series(traj::StateCountTrajectory, compartment::Symbol)

Return `(time, values)` for one compartment in a recovered trajectory.

`compartment` must be one of `:S`, `:E`, `:I`, or `:R`. Repeated event times are
preserved exactly and no interpolation or resampling is applied. The returned
vectors alias storage from `traj`; callers should treat them as read-only.
"""
function trajectory_series(traj::StateCountTrajectory, compartment::Symbol)
    values = _trajectory_values(traj, compartment)
    return (traj.time, values)
end


"""
    trajectory_series(traj::StateCountTrajectory)

Return a named tuple containing plotting-ready series for all compartments.
"""
function trajectory_series(traj::StateCountTrajectory)
    return (
        S=(traj.time, traj.S),
        E=(traj.time, traj.E),
        I=(traj.time, traj.I),
        R=(traj.time, traj.R),
    )
end


function _trajectory_values(traj::StateCountTrajectory, compartment::Symbol)
    compartment == :S && return traj.S
    compartment == :E && return traj.E
    compartment == :I && return traj.I
    compartment == :R && return traj.R
    throw(ArgumentError("compartment must be one of :S, :E, :I, :R"))
end


"""
    trajectory_final_sizes(trajs)

Return final total infected/removed outbreak size (`E + I + R`) for each
trajectory. No common time axis or aggregation is applied.
"""
function trajectory_final_sizes(trajs::AbstractVector{<:StateCountTrajectory})
    values = Vector{Int}(undef, length(trajs))
    @inbounds for i in eachindex(trajs)
        traj = trajs[i]
        values[i] = isempty(traj.time) ? 0 : traj.E[end] + traj.I[end] + traj.R[end]
    end
    return values
end


"""
    trajectory_final_times(trajs)

Return the final time from each trajectory.
"""
function trajectory_final_times(trajs::AbstractVector{<:StateCountTrajectory})
    values = Vector{Float64}(undef, length(trajs))
    @inbounds for i in eachindex(trajs)
        traj = trajs[i]
        values[i] = isempty(traj.time) ? 0.0 : traj.time[end]
    end
    return values
end


"""
    peak_infectious(trajs)

Return the maximum infectious count in each trajectory.
"""
function peak_infectious(trajs::AbstractVector{<:StateCountTrajectory})
    values = Vector{Int}(undef, length(trajs))
    @inbounds for i in eachindex(trajs)
        traj = trajs[i]
        values[i] = isempty(traj.I) ? 0 : maximum(traj.I)
    end
    return values
end


"""
    peak_infectious_time(trajs)

Return the first time of maximum infectious count in each trajectory.
"""
function peak_infectious_time(trajs::AbstractVector{<:StateCountTrajectory})
    values = Vector{Float64}(undef, length(trajs))
    @inbounds for i in eachindex(trajs)
        traj = trajs[i]
        if isempty(traj.I)
            values[i] = 0.0
        else
            peak_value = traj.I[1]
            peak_index = 1
            for j in 2:length(traj.I)
                if traj.I[j] > peak_value
                    peak_value = traj.I[j]
                    peak_index = j
                end
            end
            values[i] = traj.time[peak_index]
        end
    end
    return values
end


"""
    host_series(summary::HostEventSummary, metric::Symbol)

Return `(host_id, values)` for one host participation metric.

`metric` must be one of `:transmissions_caused`, `:transmissions`, `:samples`,
`:removals`, or `:activations`. Host ordering is the existing `summary.host_id`
ordering; no additional sorting is applied. The returned vectors alias storage
from `summary`; callers should treat them as read-only.
"""
function host_series(summary::HostEventSummary, metric::Symbol)
    values = _host_metric_values(summary, metric)
    return (summary.host_id, values)
end


function _host_metric_values(summary::HostEventSummary, metric::Symbol)
    (metric == :transmissions_caused || metric == :transmissions) && return summary.transmissions_caused
    metric == :samples && return summary.samples
    metric == :removals && return summary.removals
    metric == :activations && return summary.activations
    throw(ArgumentError("metric must be one of :transmissions_caused, :transmissions, :samples, :removals, :activations"))
end
