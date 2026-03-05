using EpiSim
using Random
using StatsPlots
using Plots
using PlotlyBase

plotly()

rng = MersenneTwister(1234)

S0 = 99999
E0 = 0
I0 = 1

λ = 2.
μ = 0.5
ψ = 0.5
ρ₀ = 0.
r = 1.

min_samples = 1_000

nrep = 10_000

model = BirthDeathModel(λ, 1. / μ, 1. / ψ, r)

stopping_criterion = StopWhenCumulativeSampled(min_samples)

ensemble = simulate_ensemble(GillespieEngine(), model, S0, E0, I0; rng=rng, stopping_criterion=stopping_criterion, nrep=nrep)

# Filter for valid simulations
valid_sims = filter(sim -> nsampled(sim) >= min_samples, ensemble)

trees = [extract_sampled_tree(sim) for sim in valid_sims]

mles = fit_bd_ensemble_mle(trees; fixed=(:ψ, ψ), r=r, ρ₀=ρ₀)

truth = Dict(:λ => λ, :μ => μ, :ψ => ψ, :R0 => λ/(μ+r*ψ), :δ => λ - (μ + r*ψ))

bd_pairs_plot(mles, r=r, truth=truth, bins=30)