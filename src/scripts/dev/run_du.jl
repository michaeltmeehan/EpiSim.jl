using EpiSim

N  = 5
β  = 2.0
α  = 1.0
I0 = 1
S0 = N - I0

gen = SIRGenerator(N, α, β)

dim = (N+1)^2
p0  = zeros(dim)
idx = I0 + S0*(N+1) + 1
p0[idx] = 1.0

p = EpiSim.DU.uniformize(gen, p0, 1.)


m = EpiSim.DU.infected_marginal(p, gen)
sum(m)

model = SIRModel(β, 1. / α, 1e9, 0.)

summary_fn = log -> reconstruct_on_grid(
    log, model, S0, E0, I0, [1.]
    )

E0 = 0

sellke_ensemble = simulate_ensemble(SellkeEngine(), model, S0, E0, I0; nrep=100_000, rng=MersenneTwister(1234), summary_fn = summary_fn)

sellke_distr = infected_distribution_matrix(sellke_ensemble, N)

gillespie_ensemble = simulate_ensemble(GillespieEngine(), model, S0, E0, I0; nrep=100_000, rng=MersenneTwister(1234), summary_fn = summary_fn)

gillespie_distr = infected_distribution_matrix(gillespie_ensemble, N)



#########################################
# Try and work out why new uniformization is wrong
#########################################
v = rand(dim)
v ./= sum(v)

out_legacy = similar(v)
Q_mul!(out_legacy, v, N, α, β, Q_bands)


out_new = similar(v)
EpiSim.DU.Q_mul!(out_new, v, gen)


p_legacy = uniformize(p0, 1., γ, Q_mul!, N, α, β, Q_bands)

p_new = EpiSim.DU.uniformize(gen, p0, 1.0)


m_new = EpiSim.DU.infected_marginal(p_new, gen)
m_legacy = aggregate_p(p_legacy, N; dims=2) |> vec

maximum(abs.(m_new .- m_legacy))