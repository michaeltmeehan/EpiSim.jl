# Add plot function for outbreak type
using RecipesBase, RecipesPipeline


function plot_outbreak(outbreak::Outbreak)
    p = plot(outbreak.traj[1, :], permutedims(sum(outbreak.traj[2:end, :], dims=1)), label="Population size", xlabel="Time", ylabel="Population size", title="Outbreak Simulation")
    scatter!(p, outbreak.linelist.t_birth, outbreak.linelist.child_id, label="Births", marker=:diamond)
    scatter!(p, outbreak.linelist.t_death, outbreak.linelist.child_id, label="Deaths", marker=:xcross)
    display(p)
end

# plot.jl

# Define the recipe function for the Outbreak type
@recipe function f(outbreak::Outbreak; show_transmission=true)

    layout := @layout [prevalence transmission
                        gannt_chart]

    # Prevalence plot
    @series begin
        grid := false
        legend := nothing
        linetype := :steppost
        xlabel := "Time"
        ylabel := "Prevalence"
        subplot := 1
        outbreak.traj[1,:], permutedims(sum(outbreak.traj[2:end, :], dims=1))
    end

    t_final = ceil(maximum(vcat(outbreak.linelist.t_birth, outbreak.linelist.t_death, outbreak.linelist.t_sam)))
    id_final = maximum(outbreak.linelist.child_id)
    scale = minimum([1., 10 / id_final])

    # Incidence plot
    @series begin
        grid := false
        legend := nothing
        seriestype := histogram
        subplot := 2
        bins := 0:t_final
        ylabel := "Incidence"
        xlabel := "Time"
        outbreak.linelist.t_birth
    end

    # Gannt chart
    # Infection lifespan
    for row in eachrow(outbreak.linelist)
        @series begin
            grid := false
            legend := nothing
            lw := 8 * scale
            alpha := row.t_sam ≥ 0. ? 1. : 0.3
            # yticks := 1:id_final
            subplot := 3
            xlabel := "Time"
            id = row.child_id
            color := id
            t_start = row.t_birth
            t_end = minimum([row.t_death, t_final])
            [t_start, t_end], [id, id]
        end
    end

    # Infection time
    for row in eachrow(outbreak.linelist)
        @series begin
            grid := false
            id = row.child_id
            legend := nothing
            marker := :circle
            markersize := 6 * scale
            markerstrokewidth := 0
            color := id
            alpha := 0.6
            # yticks := 1:id_final
            subplot := 3
            color := id
            t_start = row.t_birth
            [t_start], [id]
        end
    end


    for row in eachrow(outbreak.linelist)
        @series begin
            grid := false
            id = row.child_id
            legend := nothing
            marker := :circle
            markersize := 3 * scale
            markerstrokewidth := 0
            alpha := 0.8
            # yticks := 1:id_final
            subplot := 3
            color := :white
            t_start = row.t_birth
            [t_start], [id]
        end
    end

    # Transmission layer
    if show_transmission
        for row in eachrow(outbreak.linelist)
            if row.parent_id != 0
                @series begin
                    grid := false
                    alpha := 0.3
                    legend := nothing
                    # yticks := 1:id_final
                    color := row.parent_id
                    ylabel := "ID"
                    xlabel := "Time"
                    subplot := 3
                    t_start = t_end = row.t_birth
                    y_start = row.parent_id
                    y_end = row.child_id
                    [t_start, t_end], [y_start, y_end-0.4]
                end
            end
        end
    end
end


# TODO: Incidence time series

# TODO: Transmission tree plot
# Notes: Use GraphPlot.jl and Graphs.jl

# TODO: Infection Gantt chart

