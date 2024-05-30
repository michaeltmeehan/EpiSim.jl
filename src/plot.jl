# Add plot function for outbreak type
using RecipesBase, RecipesPipeline

# Helper function to generate distinct colors
function color_from_id(id, scheme=:tab20)
    return id
    # palette = get(scheme)
    # return palette[id % length(palette) + 1]
end

@recipe function f(outbreak::Outbreak; show_transmission=true)
    layout := @layout [prevalence_incidence
                        gannt_chart]
    
    # Merged prevalence / incidence plot
    grid := false
    xlabel := "Time"
    subplot := 1
    color := RGB(0.2, 0.4, 0.6)

    # Prevalence plot
    @series begin
        linetype := :steppost
        label := "Prevalence"
        outbreak.traj[1,:], permutedims(sum(outbreak.traj[2:end, :], dims=1))
    end

    infections_times = outbreak.linelist.t_birth
    death_times = [t for t in outbreak.linelist.t_death if !isinf(t)]
    t_final = ceil(maximum(vcat(infections_times, death_times)))
    id_final = maximum(outbreak.linelist.child_id)
    scale = minimum([1., 10 / id_final])

    # Incidence plot
    @series begin
        alpha := 0.6
        seriestype := histogram
        bins := 0:t_final
        label := "Incidence"
        outbreak.linelist.t_birth
    end

    # Gannt chart
    grid := false
    legend := nothing
    subplot := 2
    xlabel := "Time"
    ylabel := "ID"

    for row in eachrow(outbreak.linelist)
        id = row.child_id
        parent = row.parent_id
        t_birth = row.t_birth
        sampled = row.t_sam ≥ 0.
        markerstrokewidth := 0.

        # Infection lifespan
        @series begin
            lw := 8 * scale
            alpha := sampled ? 1. : 0.3
            color := color_from_id(id)
            t_end = minimum([row.t_death, t_final])
            [t_birth, t_end], [id, id]
        end


        @series begin
            marker := :circle
            markersize := 6 * scale
            color := color_from_id(id)
            alpha := 0.6
            [t_birth], [id]
        end

        @series begin
            marker := :circle
            markersize := 3 * scale
            color := :white
            alpha := 0.6
            [t_birth], [id]
        end

        if show_transmission
            if parent != 0
                @series begin
                    lw := scale
                    alpha := 0.3
                    color := color_from_id(id)
                    [t_birth, t_birth], [parent, id-0.4]
                end

                @series begin
                    marker := :circle
                    markersize := 3 * scale
                    alpha := 0.6
                    color := color_from_id(id)
                    [t_birth], [parent]
                end
            end
        end
    end
end