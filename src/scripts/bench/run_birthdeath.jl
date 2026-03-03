using EpiSim
using Random

rng = Random.MersenneTwister(1234)

S0 = 99999
E0 = 0
I0 = 1
t_grid = collect(0.:0.1:1.)

model = BirthDeathModel(2., 2., 2., 1.)

summary_fn = log -> reconstruct_on_grid(                                                                                                                 
           log, model, S0, E0, I0, t_grid                                                                                                       
           )

stopping_criterion = EpiSim.StopWhenTimeReached(1.)

ensemble = simulate_ensemble(GillespieEngine(), model, S0, E0, I0; nrep=100_000, rng=MersenneTwister(1234), stopping_criterion=stopping_criterion, summary_fn=summary_fn)

infected_distribution_matrix(ensemble, 10)

pₙ(collect(0:5), tᵢ, collect(0.:0.1:0.5), λ, μ, ψ, r)