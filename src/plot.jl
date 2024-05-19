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
@recipe function f(outbreak::Outbreak)
    xlabel := "Time"
    ylabel := "Population size"
    title := "Outbreak Simulation"
    outbreak.traj[1,:], permutedims(sum(outbreak.traj[2:end, :], dims=1))
end


# TODO: Incidence time series

# TODO: Transmission tree plot
# Notes: Use GraphPlot.jl and Graphs.jl

# TODO: Infection Gantt chart

