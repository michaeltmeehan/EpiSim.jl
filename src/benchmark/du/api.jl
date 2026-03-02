# ------------------------------------------------------------
# Public DU API for SIR
# ------------------------------------------------------------

function du_infected_distribution(model,
                                   S0::Int,
                                   I0::Int,
                                   t_grid::Vector{Float64})

    N = model.N
    α = model.recovery_rate
    β = model.transmission_rate

    gen = SIRGenerator(N, α, β)

    # Initial distribution
    dim = (N + 1)^2
    p0 = zeros(Float64, dim)

    # Index mapping: I + S*(N+1) + 1
    idx = I0 + S0 * (N + 1) + 1
    p0[idx] = 1.0

    n_t = length(t_grid)
    result = zeros(Float64, n_t, N+1)

    for (i, t) in enumerate(t_grid)
        p = uniformize(gen, p0, t)
        result[i, :] .= infected_marginal(p, gen)
    end

    return result
end