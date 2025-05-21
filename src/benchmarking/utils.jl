import StatsBase: countmap

function get_distribution(a::Matrix{Int}, max::Int)
    n = size(a, 1)
    m = size(a, 2)
    dist = zeros(n, max+1)
    for i in 1:n
        freqs = countmap(a[i, :])
        for (k, v) in freqs
            dist[i, k + 1] = v / m
        end
    end
    return dist
end


function uniformize(p0::Vector{Float64}, t::Float64, γ::Float64,
                    Q_mul!::Function, Q_args...; ε=1e-10, n_max=1000)

    p = zeros(length(p0))
    q = copy(p0)
    w = 1.0
    n = 0
    norm_deficit = 1.0

    while norm_deficit > ε && n < n_max
        p .+= exp(-γ * t) * w * q

        # Apply Q·q → temp, then P·q = q + (1/γ) * Q·q
        temp = similar(q)
        Q_mul!(temp, q, Q_args...)
        q .+= (1 / γ) .* temp

        # Update Poisson weight
        n += 1
        w *= γ * t / n

        norm_deficit = 1 - sum(p)  # L1 mass deficit
    end

    return p
end


function uniformize_leaky(p0::Vector{Float64}, t::Float64, γ::Float64,
                          Q_mul!::Function, Q_args...; ε=1e-12, n_max=1000)

    p = zeros(length(p0))
    q = copy(p0)
    n = 0
    w = 1.0  # Poisson weight = (γt)^n / n!
    expfactor = exp(-γ * t)

    while w * expfactor > ε && n < n_max
        p .+= w * expfactor * q

        temp = similar(q)
        Q_mul!(temp, q, Q_args...)
        q .+= (1 / γ) .* temp

        n += 1
        w *= γ * t / n
    end

    return p
end


function aggregate_p(p::Vector{Float64}, N::Int; dims::Int=2)
    # Reshape p into a matrix
    P = reshape(p, N+1, N+1)

    # Sum over the rows to get the aggregate distribution
    aggregate_dist = sum(P, dims=dims)

    return aggregate_dist
end