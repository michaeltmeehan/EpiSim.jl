# ------------------------------------------------------------
# Uniformization solver (model-agnostic)
# ------------------------------------------------------------

function uniformize(gen::AbstractGenerator,
                    p0::Vector{Float64},
                    t::Float64;
                    ε::Float64 = 1e-12,
                    n_max::Int = 10000)::Vector{Float64}

    γ = uniformization_rate(gen)

    p = zeros(Float64, length(p0))
    q = copy(p0)

    temp = similar(p0)

    w = 1.0
    n = 0

    expfactor = exp(-γ * t)
    poisson_mass = 0.0

    while poisson_mass < 1.0 - ε && n < n_max

        term_weight = w * expfactor
        poisson_mass += term_weight

        @inbounds p .+= term_weight .* q

        Q_mul!(temp, q, gen)
        @inbounds q .+= (1.0 / γ) .* temp

        n += 1
        w *= γ * t / n
    end

    return p
end