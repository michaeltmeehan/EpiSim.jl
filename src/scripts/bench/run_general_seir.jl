using EpiSim
using Random
using Distributions
using Statistics
using StatsBase
using Plots

rng = MersenneTwister(1234)

############################################################
# Initial conditions
############################################################

S₀ = 9999
E₀ = 0
I₀ = 1

############################################################
# Model specification (direct distributions)
############################################################

βdist = LogNormal(log(2.0),0.5)

infectious_dist = Gamma(3,1/3)

sampling_dist = Gamma(3,10.0)

exposed_dist = DiscreteNonParametric([0.0],[1.0])

model = SEIRModel(
    βdist,
    exposed_dist,
    infectious_dist,
    sampling_dist,
    1.0
)

############################################################
# Simulation settings
############################################################

nsim = 500

############################################################
# Storage
############################################################

infectious_samples = Float64[]
# normalized_offspring_samples = Float64[]
early_offspring = Int[]

############################################################
# Simulate epidemics
############################################################

for i in 1:nsim

    sim = simulate(
        SellkeEngine(),
        model,
        S₀,E₀,I₀;
        rng=rng
    )

    ll = LineList(sim)

    append!(infectious_samples, infectious_periods(ll))

    # append!(
    #     normalized_offspring_samples,
    #     normalized_offspring(ll, sim, S₀, E₀, I₀)
    # )

    df = ll.hosts
    counts = secondary_cases(ll)

    mask = .!ismissing.(df.infection_time) .& (df.infection_time .< 1.0)
    early = counts[mask]

    append!(early_offspring, early)

end


histogram(
    infectious_samples,
    bins=50,
    normalize=true,
    label="Simulation",
)

plot!(
    x -> pdf(infectious_dist,x),
    0,
    maximum(infectious_samples),
    lw=3,
    label="Theoretical"
)



theoretical_secondary = Int[]

for i in 1:100000

    β = rand(βdist)
    T = rand(infectious_dist)

    push!(theoretical_secondary, rand(Poisson(β*T)))

end


p1 = fit(Histogram, early_offspring, 0:20)
p2 = fit(Histogram, theoretical_secondary, 0:20)

bar(
    p1.edges[1][1:end-1],
    p1.weights/sum(p1.weights),
    label="Simulation",
)

bar!(
    p2.edges[1][1:end-1],
    p2.weights/sum(p2.weights),
    alpha=0.5,
    label="Theory"
)