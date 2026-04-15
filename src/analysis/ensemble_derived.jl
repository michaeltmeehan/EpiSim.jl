function _retained_logs(summary::EnsembleSummary)
    summary.logs === nothing && throw(ArgumentError("EnsembleSummary does not retain raw event logs; rerun with retain_logs=true"))
    return summary.logs
end


"""
    ensemble_state_trajectories(summary; S0, E0, I0, R0=0, validate=true)

Recover one event-time [`StateCountTrajectory`](@ref) per retained event log in
an [`EnsembleSummary`](@ref).

Raw logs must have been retained by `run_ensemble(...; retain_logs=true)`.
Initial counts are explicit and are passed unchanged to
[`event_time_state_counts`](@ref) for every log. No common time grid,
interpolation, or aggregation is applied.
"""
function ensemble_state_trajectories(summary::EnsembleSummary; S0::Integer, E0::Integer, I0::Integer, R0::Integer=0, validate::Bool=true)
    logs = _retained_logs(summary)
    trajectories = Vector{StateCountTrajectory}(undef, length(logs))
    @inbounds for i in eachindex(logs)
        trajectories[i] = event_time_state_counts(logs[i]; S0=S0, E0=E0, I0=I0, R0=R0, validate=validate)
    end
    return trajectories
end


"""
    ensemble_host_event_summaries(summary)

Return one [`HostEventSummary`](@ref) per retained event log in an
[`EnsembleSummary`](@ref).

Raw logs must have been retained by `run_ensemble(...; retain_logs=true)`.
This is a thin ensemble-scale application of [`host_event_summary`](@ref); it
does not aggregate host summaries across replicates.
"""
function ensemble_host_event_summaries(summary::EnsembleSummary)
    logs = _retained_logs(summary)
    summaries = Vector{HostEventSummary}(undef, length(logs))
    @inbounds for i in eachindex(logs)
        summaries[i] = host_event_summary(logs[i])
    end
    return summaries
end
