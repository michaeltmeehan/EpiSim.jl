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