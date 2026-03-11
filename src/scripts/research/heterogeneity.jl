using EpiSim
using Random
using Distributions
using Tables

rng = MersenneTwister(1234)

############################################################
# Initial conditions
############################################################

S₀ = 9999
E₀ = 0
I₀ = 1

############################################################
# Fixed parameters
############################################################

R₀ = 2.0
removal_rate = 1.0
sampling_proportion = 0.1
r = 1.0

############################################################
# Grid parameters
############################################################

infectious_period_shape = 1:3
transmission_heterogeneity = 0.0:1.0:2.0

############################################################
# Simulation constraints
############################################################

target_valid = 10
min_samples = 100
stopping_criterion = StopWhenCumulativeSampled(min_samples)

############################################################
# Results container
############################################################

results = Dict{Tuple{Int,Float64}, Vector}()

############################################################
# Grid loop
############################################################

for k in infectious_period_shape
    for σβ in transmission_heterogeneity

        println("Running k=$k  σβ=$σβ")

        params = epidemic_parameter_mapping(
            R₀,
            removal_rate,
            sampling_proportion,
            k,
            σβ
        )

        exposed_distribution = DiscreteNonParametric([0.0], [1.0])

        model = SEIRModel(
            params.transmission_distribution,
            exposed_distribution,
            params.recovery_distribution,
            params.sampling_distribution,
            r
        )

        valid_sims = Vector{Any}()
        counter = 0

        while length(valid_sims) < target_valid && counter < 100 * target_valid
            counter += 1
            sim = simulate(
                SellkeEngine(),
                model,
                S₀,
                E₀,
                I₀;
                rng = rng,
                stopping_criterion = stopping_criterion
            )

            if nsampled(sim) >= min_samples
                push!(valid_sims, sim)
            end

        end

        results[(k, σβ)] = valid_sims

    end
end


############################################################
# Extract sampled trees
############################################################
trees = Dict{Tuple{Int,Float64}, Vector{Tree}}()
for (k, σβ) in keys(results)
    trees[(k, σβ)] = [extract_sampled_tree(sim) for sim in results[(k, σβ)]]
end



############################################################
# Calculate tree statistics
############################################################
tree_stats = Dict{Tuple{Int,Float64}, NamedTuple}()
for (k, σβ) in keys(trees)
    tstats = [tree_statistics(tree) for tree in trees[(k, σβ)]]
    tree_stats[(k, σβ)] = Tables.columntable(tstats)
end

using DataFrames

rows = DataFrame()

for ((k, σβ), stats) in tree_stats

    df = DataFrame(stats)
    df.k .= k
    df.σβ .= σβ

    append!(rows, df)

end

stats_long = stack(
    rows,
    Not([:k, :σβ]),
    variable_name = :statistic,
    value_name = :value
)

using StatsPlots

σβ_levels = sort(unique(stats_long.σβ))
stats = unique(stats_long.statistic)

nrows = length(σβ_levels)
ncols = length(stats)

plt = plot(layout = (nrows, ncols), size=(350*ncols, 250*nrows))

for (i, σβ) in enumerate(σβ_levels)
    for (j, stat) in enumerate(stats)

        df = stats_long[
            (stats_long.σβ .== σβ) .&
            (stats_long.statistic .== stat),
            :
        ]

        idx = (i - 1)*ncols + j

        violin!(
            plt[idx],
            df.k,
            df.value,
            xlabel = "k",
            ylabel = stat,
            legend = false
        )

        title!(plt[idx], "σβ = $σβ")
    end
end

plt

############################################################
# Fit single-type BD mles to each tree
############################################################
param = R0DeltaParameterization(BDFixedSpec(:p, 0.1))

mles = Dict{Tuple{Int, Float64},  Vector{NamedTuple{(:R₀, :δ, :p), Tuple{Float64, Float64, Float64}}}}()

for (k, σβ) in keys(trees)
    mles[(k, σβ)] = [fit_bd_pars(tree; param=param, r=r) for tree in trees[(k, σβ)]]
end



############################################################
# Plot mles
############################################################
using DataFrames

rows = DataFrame(
    k = Int[],
    σβ = Float64[],
    parameter = String[],
    estimate = Float64[]
)

for ((k, σβ), vec) in mles
    for est in vec
        push!(rows, (k, σβ, "R0", est.R₀))
        push!(rows, (k, σβ, "δ", est.δ))
    end
end

using StatsPlots

σβ_levels = sort(unique(rows.σβ))
params = ["R0", "δ"]

nrows = length(σβ_levels)
ncols = length(params)

plt = plot(layout = (nrows, ncols), size=(900, 300*nrows))

for (i, σβ) in enumerate(σβ_levels)
    for (j, param) in enumerate(params)

        df = rows[(rows.σβ .== σβ) .& (rows.parameter .== param), :]

        idx = (i - 1)*ncols + j

        violin!(
            plt[idx],
            df.k,
            df.estimate,
            xlabel = "k",
            ylabel = param,
            legend = false
        )

        if param == "R0"
            hline!(plt[idx], [R₀], linestyle=:dash)
        else
            hline!(plt[idx], [removal_rate], linestyle=:dash)
        end

        title!(plt[idx], "σβ = $σβ")
    end
end

plt