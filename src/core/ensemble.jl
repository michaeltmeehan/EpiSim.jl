# =========================================
# Ensemble simulation framework
# =========================================

function simulate_ensemble(engine,
                           model,
                           S0::Int,
                           E0::Int,
                           I0::Int;
                           nrep::Int,
                           rng::AbstractRNG,
                           summary_fn = identity)

    # Run first simulation to determine return type
    log1 = simulate(engine, model, S0, E0, I0; rng=rng)
    result1 = summary_fn(log1)

    T = typeof(result1)
    results = Vector{T}(undef, nrep)
    results[1] = result1

    # Remaining simulations
    for r in 2:nrep
        log = simulate(engine, model, S0, E0, I0; rng=rng)
        results[r] = summary_fn(log)
    end

    return results
end