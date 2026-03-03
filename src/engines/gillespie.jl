# =========================================
# Gillespie engine (Markovian count-based)
# =========================================

using Distributions
using Random

function simulate(::GillespieEngine,
                  model::AbstractEpidemicModel,
                  S0::Int,
                  E0::Int,
                  I0::Int;
                  rng::AbstractRNG = Random.default_rng(),
                  stopping_criterion = default_stopping(model))

    # Enforce exponential assumptions
    if !(model.dτi isa Exponential)
        error("GillespieEngine requires exponential infectious period")
    end

    if has_latent_stage(model) &&
       !(model.dτe isa Exponential)
        error("GillespieEngine requires exponential incubation period")
    end

    N = S0 + E0 + I0

    S = S0
    E = E0
    I = I0

    t = 0.0
    log = EventLog()

    # Log seeding
    for id in 1:(E0 + I0)
        push_event!(log, t, Seeding, id, 0)
    end

    β̄ = mean(model.dβ)

    state = SimulationState(t, E0 + I0, 0)

    # For ID assignment
    next_id = N + 1

    while E + I > 0 && !should_stop(stopping_criterion, state)

        λ_inf = β̄ * S * I / N
        λ_act = has_latent_stage(model) ? (E / mean(model.dτe)) : 0.0
        λ_rem = I / mean(model.dτi)

        λ_total = λ_inf + λ_act + λ_rem

        if λ_total == 0.0
            break
        end

        # Time step
        Δt = randexp(rng) / λ_total
        t += Δt

        u = rand(rng) * λ_total

        if u < λ_inf
            # Infection
            if S > 0
                S -= 1
                if has_latent_stage(model)
                    E += 1
                else
                    I += 1
                end

                push_event!(log, t, Transmission, next_id, 0)
                next_id += 1
                update_state!(state, t, Transmission)
            end

        elseif u < λ_inf + λ_act
            # Activation
            if E > 0
                E -= 1
                I += 1
                push_event!(log, t, Activation, 0, 0)
                update_state!(state, t, Activation)
            end

        else
            # Removal
            if I > 0
                I -= 1
                push_event!(log, t, Removal, 0, 0)
                update_state!(state, t, Removal)
            end
        end
    end

    return log
end