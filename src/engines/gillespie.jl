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

    model_has_latent = has_latent_stage(model)
    if model_has_latent && !(model.dτe isa Exponential)
        error("GillespieEngine requires exponential incubation period")
    end

    if !(model.dτs isa Exponential)
        error("GillespieEngine requires exponential sampling interval")
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

    exposed_ids = collect(1:E0)
    active_ids = collect((E0 + 1):(E0 + I0))

    β_bar = mean(model.dβ)
    σ_bar = model_has_latent ? 1 / mean(model.dτe) : 0.0
    α_bar = 1 / mean(model.dτi)
    ψ_bar = 1 / mean(model.dτs)

    state = SimulationState(t, E0 + I0, 0)

    # For ID assignment
    next_id = E0 + I0 + 1

    while E + I > 0 && !should_stop(stopping_criterion, state)

        if model isa BirthDeathModel
            λ_inf = β_bar * I
        else
            λ_inf = β_bar * S * I / N
        end
        λ_act = σ_bar * E
        λ_rem = α_bar * I
        λ_sam = ψ_bar * I

        λ_total = λ_inf + λ_act + λ_rem + λ_sam

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
                infector = rand(rng, active_ids)
                push_event!(log, t, Transmission, next_id, infector)
                S -= 1
                if model_has_latent
                    E += 1
                    push!(exposed_ids, next_id)
                else
                    I += 1
                    push!(active_ids, next_id)
                end
                update_state!(state, t, Transmission)
                next_id += 1
            end

        elseif u < λ_inf + λ_act
            # Activation
            if E > 0
                E -= 1
                I += 1
                activated_id = popr!(rng, exposed_ids)
                push!(active_ids, activated_id)
                push_event!(log, t, Activation, activated_id, 0)
                update_state!(state, t, Activation)
            end

        elseif u < λ_inf + λ_act + λ_rem
            # Removal
            if I > 0
                I -= 1
                removed_id = popr!(rng, active_ids)
                push_event!(log, t, Removal, removed_id, 0)
                update_state!(state, t, Removal)
            end
        else
            # Sampling
            if I > 0
                if rand(rng) < model.r
                    I -= 1
                    sampled_id = popr!(rng, active_ids)
                    ev = SerialSampling
                else
                    sampled_id = rand(rng, active_ids)
                    ev = FossilisedSampling
                end
                push_event!(log, t, ev, sampled_id, 0)
                update_state!(state, t, ev)
            end
        end
    end

    return log
end