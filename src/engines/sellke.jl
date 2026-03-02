# =========================================
# Minimal Sellke engine (clean heap design)
# =========================================

using DataStructures: BinaryMinHeap
using Distributions

function simulate(::SellkeEngine,
                  model::AbstractEpidemicModel,
                  S0::Int,
                  E0::Int,
                  I0::Int;
                  rng::AbstractRNG = Random.default_rng())

    N = S0 + E0 + I0

    # ---------------------------------
    # Heaps
    # ---------------------------------

    resistances = BinaryMinHeap{Tuple{Float64,Int}}()
    temporal    = BinaryMinHeap{Tuple{Float64,EventKind,Int}}()

    # ---------------------------------
    # Infection pressure bookkeeping
    # ---------------------------------

    β = zeros(Float64, N)
    active_ids = Int[]
    total_beta = 0.0

    # ---------------------------------
    # Time & cumulative pressure
    # ---------------------------------

    t = 0.0
    Λ = 0.0

    log = EventLog()

    # ---------------------------------
    # Initialise susceptible resistances
    # ---------------------------------

    for id in (E0 + I0 + 1):N
        push!(resistances, (randexp(rng), id))
    end

    # ---------------------------------
    # Seed exposed
    # ---------------------------------

    for id in 1:E0
        push_event!(log, t, Seeding, id, 0)

        if has_latent_stage(model)
            τe = draw_tau_e(rng, model)
            push!(temporal, (t + τe, Activation, id))
        else
            activate_host!(rng, temporal, model, id, t, β, active_ids, total_beta)
        end
    end

    # ---------------------------------
    # Seed infected
    # ---------------------------------

    for id in (E0 + 1):(E0 + I0)
        push_event!(log, t, Seeding, id, 0)
        total_beta = activate_host!(rng, temporal, model, id, t, β, active_ids, total_beta)
    end

    # ---------------------------------
    # Main loop
    # ---------------------------------

    while !isempty(temporal)

        # Next temporal event
        t_next, ev, id = top(temporal)

        dΛ = total_beta / N

        # ---------------------------------
        # Infection event?
        # ---------------------------------

        if !isempty(resistances) && dΛ > 0.0

            r_threshold, sid = top(resistances)
            Λ_next = Λ + dΛ * (t_next - t)

            if r_threshold ≤ Λ_next
                # Infection occurs first
                pop!(resistances)

                t += (r_threshold - Λ) / dΛ
                Λ = r_threshold

                parent = sample_infector(rng, active_ids, β, total_beta)

                push_event!(log, t, Transmission, sid, parent)

                if has_latent_stage(model)
                    τe = draw_tau_e(rng, model)
                    push!(temporal, (t + τe, Activation, sid))
                else
                    activate_host!(rng, temporal, model, sid, t,
                                   β, active_ids, total_beta)
                end

                continue
            end
        end

        # ---------------------------------
        # Temporal event happens
        # ---------------------------------

        pop!(temporal)

        Λ += dΛ * (t_next - t)
        t = t_next

        push_event!(log, t, ev, id, 0)

        if ev == Activation

            total_beta = activate_host!(rng, temporal, model, id, t,
                                        β, active_ids, total_beta)

        elseif ev == Removal || ev == SerialSampling

            total_beta -= β[id]
            β[id] = 0.0
            remove_active!(active_ids, id)

        elseif ev == FossilisedSampling
            # Nothing structural changes
        end
    end

    return log
end


# ====================================================
# Host activation (adds infectious pressure + schedules events)
# ====================================================

function activate_host!(rng::AbstractRNG,
                        temporal,
                        model,
                        id,
                        t,
                        β,
                        active_ids,
                        total_beta)

    β[id] = draw_beta(rng, model)
    total_beta += β[id]
    push!(active_ids, id)

    schedule_infectious_events!(rng, temporal, model, id, t)

    return total_beta
end


# ====================================================
# Schedule all infectious-phase events at activation
# ====================================================

function schedule_infectious_events!(rng::AbstractRNG,
                                     temporal,
                                     model,
                                     id,
                                     t_activation)

    τi = draw_tau_i(rng, model)

    # Sampling process
    τs = draw_tau_s(rng, model)

    while τs < τi
        if rand(rng) < model.r
            push!(temporal, (t_activation + τs, SerialSampling, id))
            return
        else
            push!(temporal, (t_activation + τs, FossilisedSampling, id))
        end
        τs += draw_tau_s(rng, model)
    end

    # Final removal
    push!(temporal, (t_activation + τi, Removal, id))
end


# ====================================================
# Sample infector proportional to β
# ====================================================

function sample_infector(rng::AbstractRNG, active_ids, β, total_beta)

    u = rand(rng) * total_beta
    s = 0.0

    for id in active_ids
        s += β[id]
        if s ≥ u
            return id
        end
    end

    return active_ids[end]
end


# ====================================================
# Remove from active list
# ====================================================

function remove_active!(active_ids, id)
    idx = findfirst(==(id), active_ids)
    idx !== nothing && deleteat!(active_ids, idx)
end