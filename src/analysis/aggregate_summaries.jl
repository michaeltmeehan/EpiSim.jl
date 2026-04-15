"""
    ScalarSummary

Mean, minimum, and maximum for a finite scalar collection.
"""
struct ScalarSummary
    mean::Float64
    minimum::Float64
    maximum::Float64
end


function Base.show(io::IO, summary::ScalarSummary)
    print(io, "ScalarSummary(mean=", summary.mean, ", minimum=", summary.minimum, ", maximum=", summary.maximum, ")")
end


function _show_scalar_line(io::IO, label::AbstractString, summary::ScalarSummary)
    print(io, "  ", label, ": mean=", summary.mean, ", min=", summary.minimum, ", max=", summary.maximum)
end


"""
    EnsembleAggregateSummary

Interpreted aggregate summaries over per-replicate vectors in an
[`EnsembleSummary`](@ref).

Fields whose type is [`ScalarSummary`](@ref) summarize one value per replicate
using mean, minimum, and maximum. `total_samples` is different: it is the
ensemble-wide total number of fossilized plus serial sampling events across all
replicates.
"""
struct EnsembleAggregateSummary
    final_size::ScalarSummary
    final_time::ScalarSummary
    transmissions::ScalarSummary
    activations::ScalarSummary
    removals::ScalarSummary
    samples::ScalarSummary
    total_samples::Int
end


function Base.show(io::IO, summary::EnsembleAggregateSummary)
    print(io, "EnsembleAggregateSummary(")
    print(io, "\n")
    _show_scalar_line(io, "final outbreak size", summary.final_size)
    print(io, "\n")
    _show_scalar_line(io, "final event time", summary.final_time)
    print(io, "\n")
    _show_scalar_line(io, "transmission events per replicate", summary.transmissions)
    print(io, "\n")
    _show_scalar_line(io, "activation events per replicate", summary.activations)
    print(io, "\n")
    _show_scalar_line(io, "removal events per replicate", summary.removals)
    print(io, "\n")
    _show_scalar_line(io, "sampling events per replicate", summary.samples)
    print(io, "\n  total samples: ", summary.total_samples)
    print(io, "\n)")
end


"""
    TrajectoryAggregateSummary

Interpreted aggregate summaries over scalar values derived independently from
each [`StateCountTrajectory`](@ref). No common time axis is constructed.

`peak_infectious_time` is the first time at which the peak infectious count is
reached within each trajectory, then summarized across trajectories.
"""
struct TrajectoryAggregateSummary
    peak_infectious::ScalarSummary
    peak_infectious_time::ScalarSummary
    final_removed::ScalarSummary
    final_time::ScalarSummary
end


function Base.show(io::IO, summary::TrajectoryAggregateSummary)
    print(io, "TrajectoryAggregateSummary(")
    print(io, "\n")
    _show_scalar_line(io, "peak infectious count", summary.peak_infectious)
    print(io, "\n")
    _show_scalar_line(io, "first time of peak infectious count", summary.peak_infectious_time)
    print(io, "\n")
    _show_scalar_line(io, "final removed count", summary.final_removed)
    print(io, "\n")
    _show_scalar_line(io, "final trajectory time", summary.final_time)
    print(io, "\n)")
end


"""
    HostAggregateSummary

Interpreted aggregate summaries over per-replicate [`HostEventSummary`](@ref)
objects.

Per-host quantities are computed within each replicate first, then the resulting
replicate-level values are summarized across replicates. Host IDs are not
matched, merged, or pooled across the ensemble.
"""
struct HostAggregateSummary
    observed_hosts::ScalarSummary
    mean_transmissions_per_host::ScalarSummary
    max_transmissions_per_host::ScalarSummary
    mean_samples_per_host::ScalarSummary
    mean_removals_per_host::ScalarSummary
    mean_activations_per_host::ScalarSummary
end


function Base.show(io::IO, summary::HostAggregateSummary)
    print(io, "HostAggregateSummary(")
    print(io, "\n")
    _show_scalar_line(io, "observed hosts per replicate", summary.observed_hosts)
    print(io, "\n")
    _show_scalar_line(io, "mean transmissions caused per observed host", summary.mean_transmissions_per_host)
    print(io, "\n")
    _show_scalar_line(io, "maximum transmissions caused by one host", summary.max_transmissions_per_host)
    print(io, "\n")
    _show_scalar_line(io, "mean samples per observed host", summary.mean_samples_per_host)
    print(io, "\n")
    _show_scalar_line(io, "mean removals per observed host", summary.mean_removals_per_host)
    print(io, "\n")
    _show_scalar_line(io, "mean activations per observed host", summary.mean_activations_per_host)
    print(io, "\n)")
end


"""
    scalar_summary(values)

Return mean, minimum, and maximum for a non-empty scalar collection.
"""
function scalar_summary(values::AbstractVector{<:Real})
    isempty(values) && throw(ArgumentError("cannot summarize an empty collection"))
    total = 0.0
    min_value = Float64(values[firstindex(values)])
    max_value = min_value

    @inbounds for value in values
        x = Float64(value)
        total += x
        x < min_value && (min_value = x)
        x > max_value && (max_value = x)
    end

    return ScalarSummary(total / length(values), min_value, max_value)
end


"""
    ensemble_aggregate_summary(summary::EnsembleSummary)

Aggregate primitive per-replicate fields already stored in an ensemble summary.
Sample counts combine fossilized and serial sampling events per replicate.
`total_samples` is the ensemble-wide sum of fossilized plus serial sampling
events. Empty ensembles throw `ArgumentError`.
"""
function ensemble_aggregate_summary(summary::EnsembleSummary)
    length(summary) == 0 && throw(ArgumentError("cannot summarize an empty ensemble"))
    samples = summary.fossilized_samples .+ summary.serial_samples
    return EnsembleAggregateSummary(
        scalar_summary(summary.final_size),
        scalar_summary(summary.final_time),
        scalar_summary(summary.transmissions),
        scalar_summary(summary.activations),
        scalar_summary(summary.removals),
        scalar_summary(samples),
        sum(samples),
    )
end


"""
    trajectory_aggregate_summary(trajs)

Aggregate scalar trajectory summaries computed independently for each
trajectory. No interpolation, resampling, or common-axis alignment is performed.
Empty trajectory collections throw `ArgumentError`.
"""
function trajectory_aggregate_summary(trajs::AbstractVector{<:StateCountTrajectory})
    isempty(trajs) && throw(ArgumentError("cannot summarize an empty trajectory collection"))
    final_removed = Vector{Int}(undef, length(trajs))
    @inbounds for i in eachindex(trajs)
        traj = trajs[i]
        final_removed[i] = isempty(traj.R) ? 0 : traj.R[end]
    end

    return TrajectoryAggregateSummary(
        scalar_summary(peak_infectious(trajs)),
        scalar_summary(peak_infectious_time(trajs)),
        scalar_summary(final_removed),
        scalar_summary(trajectory_final_times(trajs)),
    )
end


"""
    host_aggregate_summary(summaries)

Aggregate replicate-level host participation summaries. Per-host quantities are
computed within each replicate, then summarized across replicates. Host IDs are
not matched across replicates. Replicates with zero observed hosts contribute
`0` to the per-host mean and maximum host-participation metrics. Empty
host-summary collections throw `ArgumentError`.
"""
function host_aggregate_summary(summaries::AbstractVector{<:HostEventSummary})
    isempty(summaries) && throw(ArgumentError("cannot summarize an empty host-summary collection"))

    observed = Vector{Int}(undef, length(summaries))
    mean_transmissions = Vector{Float64}(undef, length(summaries))
    max_transmissions = Vector{Int}(undef, length(summaries))
    mean_samples = Vector{Float64}(undef, length(summaries))
    mean_removals = Vector{Float64}(undef, length(summaries))
    mean_activations = Vector{Float64}(undef, length(summaries))

    @inbounds for i in eachindex(summaries)
        summary = summaries[i]
        n = length(summary)
        observed[i] = n
        if n == 0
            mean_transmissions[i] = 0.0
            max_transmissions[i] = 0
            mean_samples[i] = 0.0
            mean_removals[i] = 0.0
            mean_activations[i] = 0.0
        else
            mean_transmissions[i] = sum(summary.transmissions_caused) / n
            max_transmissions[i] = maximum(summary.transmissions_caused)
            mean_samples[i] = sum(summary.samples) / n
            mean_removals[i] = sum(summary.removals) / n
            mean_activations[i] = sum(summary.activations) / n
        end
    end

    return HostAggregateSummary(
        scalar_summary(observed),
        scalar_summary(mean_transmissions),
        scalar_summary(max_transmissions),
        scalar_summary(mean_samples),
        scalar_summary(mean_removals),
        scalar_summary(mean_activations),
    )
end
