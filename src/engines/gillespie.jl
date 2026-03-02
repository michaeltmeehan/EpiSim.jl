# =========================================
# Gillespie engine (Markovian, count-based)
# =========================================

using Random

function simulate(::GillespieEngine,
                  model::AbstractEpidemicModel,
                  S0::Int,
                  E0::Int,
                  I0::Int)

    # Enforce Markovian assumptions
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
    R = 0

    t = 0.0
    next_host_id = N + 1

    log = EventLog()

    # Seed initial infected
    for id in 1:(E0 + I0)
        push_event!(log, t, Seeding, id, 0)
    end

    while E + I > 0

        # Infection rate
        β̄ = mean(model.dβ)  # Approximate average β
        λ_inf = β̄ * S * I / N

        # Activation rate (if latent stage)
        λ_act = has_latent_stage(model) ? (E / mean(model.dτe)) : 0.0

        # Removal rate
        λ_rem = I / mean(model.dτi)

        λ_total = λ_inf + λ_act + λ_rem

        if λ_total == 0.0
            break
        end

        # Time step
        Δt = rand(Exponential(λ_total))
        t += Δt

        u = rand() * λ_total

        if u < λ_inf
            # Infection event

            if S > 0
                S -= 1
                if has_latent_stage(model)
                    E += 1
                else
                    I += 1
                end

                parent = rand(1:max(I,1))  # placeholder infector
                push_event!(log, t, Transmission, next_host_id, parent)
                next_host_id += 1
            end

        elseif u < λ_inf + λ_act
            # Activation event
            if E > 0
                E -= 1
                I += 1
                push_event!(log, t, Activation, 0, 0)
            end

        else
            # Removal event
            if I > 0
                I -= 1
                R += 1
                push_event!(log, t, Removal, 0, 0)
            end
        end
    end

    return log
end