# =========================================
# EventLog → state trajectory reconstruction
# =========================================

export StateTrajectory,
       reconstruct_trajectory,
       reconstruct_on_grid

struct StateTrajectory
    t::Vector{Float64}
    S::Vector{Int}
    E::Vector{Int}
    I::Vector{Int}
    R::Vector{Int}
end


function reconstruct_trajectory(log::EventLog,
                                model::AbstractEpidemicModel,
                                S0::Int,
                                E0::Int,
                                I0::Int)

    n = length(log.t)

    t_vec = Vector{Float64}(undef, n)
    S_vec = Vector{Int}(undef, n)
    E_vec = Vector{Int}(undef, n)
    I_vec = Vector{Int}(undef, n)
    R_vec = Vector{Int}(undef, n)

    S = S0
    E = E0
    I = I0
    R = 0

    for i in 1:n

        t_vec[i] = log.t[i]
        ev = log.kind[i]

        if ev == Transmission
            S -= 1
            # If no latent stage this will be corrected externally
            if has_latent_stage(model)
                E += 1
            else
                I += 1
            end

        elseif ev == Activation
            E -= 1
            I += 1

        elseif ev == Removal || ev == SerialSampling
            I -= 1
            R += 1

        end

        S_vec[i] = S
        E_vec[i] = E
        I_vec[i] = I
        R_vec[i] = R
    end

    return StateTrajectory(t_vec, S_vec, E_vec, I_vec, R_vec)
end


function reconstruct_on_grid(log::EventLog,
                             model::AbstractEpidemicModel,
                             S0::Int,
                             E0::Int,
                             I0::Int,
                             t_grid::Vector{Float64})

    n_grid = length(t_grid)

    S_vec = Vector{Int}(undef, n_grid)
    E_vec = Vector{Int}(undef, n_grid)
    I_vec = Vector{Int}(undef, n_grid)
    R_vec = Vector{Int}(undef, n_grid)

    S = S0
    E = E0
    I = I0
    R = 0

    event_index = 1
    n_events = length(log.t)

    for j in 1:n_grid
        t_g = t_grid[j]

        # Process all events up to t_g
        while event_index <= n_events &&
              log.t[event_index] ≤ t_g

            ev = log.kind[event_index]

            if ev == Transmission
                S -= 1
                if has_latent_stage(model)
                    E += 1
                else
                    I += 1
                end

            elseif ev == Activation
                E -= 1
                I += 1

            elseif ev == Removal || ev == SerialSampling
                I -= 1
                R += 1
            end

            event_index += 1
        end

        S_vec[j] = S
        E_vec[j] = E
        I_vec[j] = I
        R_vec[j] = R
    end

    return StateTrajectory(t_grid, S_vec, E_vec, I_vec, R_vec)
end