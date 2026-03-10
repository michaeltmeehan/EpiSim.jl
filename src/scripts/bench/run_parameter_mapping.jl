using EpiSim
using Random
using Distributions
using Statistics
using DataFrames


rng = MersenneTwister(1234)

############################################################
# PARAMETERS TO VALIDATE
############################################################

R0 = 2.0
δ  = 1.0
p  = 0.1
k  = 2
σβ = 0.3

############################################################
# SIMULATION SETTINGS
############################################################

nsim = 500
min_samples = 100

S₀ = 9999
E₀ = 0
I₀ = 1

# stopping_criterion = StopWhenCumulativeSampled(min_samples)
stopping_criterion = NoStopping()

############################################################
# DERIVE MODEL PARAMETERS
############################################################

params = epidemic_parameter_mapping(R0, δ, p, k, σβ)

exposed_distribution = DiscreteNonParametric([0.0],[1.0])

model = SEIRModel(
    params.transmission_distribution,
    exposed_distribution,
    params.recovery_distribution,
    params.sampling_distribution,
    1.0
)

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
    x -> pdf(params.recovery_distribution,x),
    0,
    maximum(infectious_samples),
    lw=3,
    label="Theoretical"
)



theoretical_secondary = Int[]

for i in 1:100000

    β = rand(params.transmission_distribution)
    T = rand(params.recovery_distribution)

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