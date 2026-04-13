using EpiSim
using Random
using Distributions
using Tables
using DataFrames
using StatsPlots
using Statistics
using StatsBase
using CairoMakie
using PairPlots

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

removal_rate = 1.0
r = 1.0

############################################################
# Simulation settings
############################################################

k = 1
σβ = 0.0

min_samples = 100
target_valid = 1_000
stopping_criterion = StopWhenCumulativeSampled(min_samples)

############################################################
# Parameter ranges
############################################################

R0_range = (1.2, 3.5)
δ_range  = (0.2, 2.0)
p_range  = (0.05, 0.5)

############################################################
# Simulation storage
############################################################

sims = Vector{Any}()
pars = Vector{NamedTuple}()

############################################################
# Simulation loop
############################################################

counter = 0

while length(sims) < target_valid && counter < 100*target_valid

    counter += 1

    R₀ = rand(rng, Uniform(R0_range...))
    δ  = rand(rng, Uniform(δ_range...))
    p  = rand(rng, Uniform(p_range...))

    params = epidemic_parameter_mapping(
        R₀,
        δ,
        p,
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
        push!(sims, sim)
        push!(pars, (R0=R₀, δ=δ, p=p))
    end

end

println("Valid simulations = ", length(sims))

############################################################
# Extract trees
############################################################

trees = [extract_sampled_tree(sim) for sim in sims]

############################################################
# Compute tree statistics
############################################################

tstats = [tree_statistics(tree) for tree in trees]

tree_stats = Tables.columntable(tstats)

############################################################
# Fit birth–death MLEs
############################################################

mles = [
    fit_bd_pars(tree; param=R0DeltaParameterization(BDFixedSpec(:p, pars[i].p)), r=r)
    for (i,tree) in enumerate(trees)
]

############################################################
# Assemble dataset
############################################################

rows = DataFrame(tree_stats)

rows.R0_true = [p.R0 for p in pars]
rows.δ_true  = [p.δ for p in pars]
rows.p_true  = [p.p for p in pars]

rows.R0_hat = [m.R₀ for m in mles]
rows.δ_hat  = [m.δ for m in mles]
rows.p_hat  = [m.p for m in mles]

pairplot(rows)

############################################################
# Correlation analysis
############################################################

tree_stat_cols = setdiff(
    names(rows),
    [:R0_true, :δ_true, :p_true, :R0_hat, :δ_hat, :p_hat]
)

println("\nCorrelation with TRUE parameters")

for stat in tree_stat_cols
    println("\nStatistic: ", stat)
    println("R0_true  ", round(corspearman(rows.R0_true, rows[!,stat]), digits=3))
    println("δ_true   ", round(corspearman(rows.δ_true, rows[!,stat]), digits=3))
    println("p_true   ", round(corspearman(rows.p_true, rows[!,stat]), digits=3))
end

println("\nCorrelation with MLE parameters")

for stat in tree_stat_cols
    println("\nStatistic: ", stat)
    println("R0_hat  ", round(corspearman(rows.R0_hat, rows[!,stat]), digits=3))
    println("δ_hat   ", round(corspearman(rows.δ_hat, rows[!,stat]), digits=3))
    println("p_hat   ", round(corspearman(rows.p_hat, rows[!,stat]), digits=3))
end

params_true = [:R0_true, :δ_true, :p_true]

tree_stats = setdiff(
    names(rows),
    [:R0_true, :δ_true, :p_true, :R0_hat, :δ_hat, :p_hat]
)

plots = Plots.Plot[]

for stat in tree_stats

    p = plot(
        layout = (1, length(params_true)),
        size = (1200, 350),
        margin = 10Plots.mm
    )

    for (j, param) in enumerate(params_true)

        scatter!(
            p[j],
            rows[!,param],
            rows[!,stat],
            markersize = 3,
            alpha = 0.6,
            legend = false,
            xlabel = string(param),
            ylabel = string(stat)
        )

    end

    title!(p[2], string(stat))

    push!(plots, p)

end

plot(plots..., layout = (length(tree_stats), 1), size=(1200, 350*length(tree_stats)))