using EpiSim
using Random
using Distributions

rng = MersenneTwister(1234)

############################################################
# Initial conditions
############################################################

S₀ = 99999
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

infectious_period_shape = 1:5
transmission_heterogeneity = 0.0:1.0:5.0

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

        while length(valid_sims) < target_valid && counter < 1000
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