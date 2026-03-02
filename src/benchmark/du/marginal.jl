# ------------------------------------------------------------
# Marginal extraction
# ------------------------------------------------------------

function infected_marginal(p::Vector{Float64},
                           gen::SIRGenerator)::Vector{Float64}

    N = gen.N
    P = reshape(p, N+1, N+1)

    # Sum over S dimension (rows)
    return vec(sum(P, dims=2))
end


# -----------------------------------------
# Infected marginal
# -----------------------------------------

function infected_marginal(p::Vector{Float64},
                           gen::SEIRGenerator)

    N = gen.N
    m = zeros(Float64, N+1)

    @inbounds for (k, (_,_,I)) in enumerate(gen.index_state)
        m[I+1] += p[k]
    end

    return m
end