# =========================================
# EventLog → state trajectory reconstruction
# =========================================

export StateTrajectory,
       reconstruct_trajectory

struct StateTrajectory
    t::Vector{Float64}
    S::Vector{Int}
    E::Vector{Int}
    I::Vector{Int}
    R::Vector{Int}
end


function reconstruct_trajectory(log::EventLog,
                                S0::Int,
                                E0::Int,
                                I0::Int)

    n = length(log.time)

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

        t_vec[i] = log.time[i]
        ev = log.kind[i]

        if ev == Transmission
            S -= 1
            # If no latent stage this will be corrected externally
            E += 1

        elseif ev == Activation
            E -= 1
            I += 1

        elseif ev == Removal
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