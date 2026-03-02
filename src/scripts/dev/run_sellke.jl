N = 100
I0 = 1
E0 = 0
S0 = N - I0 - E0

model = SEIRModel(1.5, 2.0, 3.0, 1., 1.)

rng = MersenneTwister(1234)

log = simulate(SellkeEngine(), model, S0, E0, I0; rng)

t_grid = collect(0.:0.1:20.)

states = reconstruct_on_grid(log, model, S0, E0, I0, t_grid)

summary_fn = log -> reconstruct_on_grid(
    log, model, S0, E0, I0, t_grid
    )


ensemble = simulate_ensemble(SellkeEngine(), model, S0, E0, I0; nrep=100, rng=rng, summary_fn = summary_fn)


distr = infected_distribution_matrix(ensemble, N)