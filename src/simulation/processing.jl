import ..Models: n_sampled, n_recovered, n_transmissions, n_activations, n_seeds, n_events

event_counts(sim::AbstractSimulation) = StatsBase.countmap(typeof.(sim.event_log))

@forward Simulation.event_log n_sampled, n_recovered, n_transmissions, n_activations, n_seeds, n_events


function get_state(sim::AbstractSimulation, t::Real)
    idx = searchsortedlast(sim.state_log.t, t)
    idx == 0 &&  throw(ArgumentError("No state recorded at or before t = $t"))
    return (; t=t, sim.state_log[idx, Not(:t)]...)
end


function get_state(sim::AbstractSimulation, tvec::AbstractVector{<:Real})
    states = [get_state(sim, t) for t in tvec]
    return Tables.columntable(states)
end


# struct StateTable
#     t::AbstractVector{<:Real}
#     vars::NamedTuple
# end



function get_state(ens::AbstractEnsemble, tvec::AbstractVector{<:Real})
    states = [get_state(sim, tvec) for sim in ens]
    states = Tables.columntable(states)[Not(:t)]
    return merge((t=tvec,), states)
end





function get_prevalence_quantiles(ens::AbstractEnsemble,
                                  tvec::AbstractVector{<:Real};
                                  probs::AbstractVector{<:Real} = [0.025, 0.25, 0.5, 0.75, 0.975])
    prevalence_matrix = get_prevalence(ens, tvec)  # size: (n_time, n_replicates)
    n_t = size(prevalence_matrix, 1)
    n_q = length(probs)

    quantile_matrix = Matrix{Float64}(undef, n_t, n_q)

    for i in 1:n_t
        quantile_matrix[i, :] .= quantile(view(prevalence_matrix, i, :), probs)
    end

    return quantile_matrix
end


function plot_prevalence(ens::AbstractEnsemble,
                         tvec::AbstractVector{<:Real};
                         probs::AbstractVector{<:Real} = [0.025, 0.25, 0.5, 0.75, 0.975],
                         labels::Tuple{String, String} = ("Time", "Prevalence"),
                         title::String = "Prevalence Ribbon Plot",
                         show_plot::Bool = true)
    Q = get_prevalence_quantiles(ens, tvec; probs)

    # Determine which quantiles to plot as ribbons
    # Here: outer (95% CI) and inner (IQR)
    median     = Q[:, findfirst(==(0.5), probs)]
    lower_outer = Q[:, findfirst(==(0.025), probs)]
    upper_outer = Q[:, findfirst(==(0.975), probs)]
    lower_inner = Q[:, findfirst(==(0.25), probs)]
    upper_inner = Q[:, findfirst(==(0.75), probs)]

    # Plot with two ribbons
    p = plot(tvec, median,
             ribbon = (median .- lower_outer, upper_outer .- median),
             label = "Median (95% CI)",
             color = :blue,
             alpha = 0.2,
             xlabel = labels[1],
             ylabel = labels[2],
             title = title)

    # Add IQR ribbon
    plot!(tvec, median,
          ribbon = (median .- lower_inner, upper_inner .- median),
          label = "Median (IQR)",
          color = :blue,
          alpha = 0.4)

    return show_plot ? display(p) : p
end



function isextinct(df::DataFrame)::Bool
    if :I ∉ propertynames(df)
        throw(ArgumentError("Prevalence column `I` not found in DataFrame"))
    end
    return df.I[end] == 0.0
end


function isextinct(sim::AbstractSimulation)::Bool
    return isextinct(sim.state_log)
end


function isextinct(sims::Vector{<:AbstractSimulation})::Vector{Bool}
    return [isextinct(sim) for sim in sims]
end


function isextinct(ens::AbstractEnsemble)::Vector{Bool}
    return isextinct(ens.simulations)
end


function get_extinction_probability(ens::AbstractEnsemble)::Float64
    n_extinct = count(isextinct(ens))
    return n_extinct / length(ens.simulations)
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


