@inline function _bd_coefficients(w, λ, μ, ψ, r)
    a = muladd(r * ψ, w, μ)
    b = muladd((1 - r) * ψ, w, -(λ + μ + ψ))
    disc = muladd(-4a, λ, b*b)
    disc = ifelse(disc < zero(disc), zero(disc), disc)
    Δ = sqrt(disc)
    return a, b, Δ
end


@inline function γ(w::T, tᵢ::T, tⱼ::T, λ::T, μ::T, ψ::T, r::T) where {T<:AbstractFloat}

    τ = tⱼ - tᵢ
    τ ≤ zero(T) && return zero(T)

    a, b, Δ = _bd_coefficients(w, λ, μ, ψ, r)

    x = Δ * τ

    # Small Δ τ branch (series expansion)
    if abs(x) ≤ sqrt(eps(T))

        # Use first-order expansion of exp(-x)
        # 1 - exp(-x) ≈ x - x^2/2
        num = 2λ * (x - x*x/2)

        den = (Δ - b) + (Δ + b) * (1 - x + x*x/2)

        # stabilise denominator
        scale = max(abs(den), one(T))
        den = abs(den) < eps(T)*scale ?
              copysign(eps(T)*scale, den) :
              den

        return num / den
    end

    # General branch
    em = exp(-x)
    one_minus_em = -expm1(-x)   # 1 - exp(-x), stable

    num = 2λ * one_minus_em
    den = (Δ - b) + (Δ + b) * em

    scale = max(abs(den), one(T))
    den = abs(den) < eps(T)*scale ?
          copysign(eps(T)*scale, den) :
          den

    return num / den
end


@inline function α(w::T, tᵢ::T, tⱼ::T, λ::T, μ::T, ψ::T, r::T) where {T<:AbstractFloat}
    a = muladd(r * ψ, w, μ)
    return a / λ * γ(w, tᵢ, tⱼ, λ, μ, ψ, r)
end


@inline function β(w::T, tᵢ::T, tⱼ::T, λ::T, μ::T, ψ::T, r::T) where {T<:AbstractFloat}
    a, b, _ = _bd_coefficients(w, λ, μ, ψ, r)
    γij = γ(w, tᵢ, tⱼ, λ, μ, ψ, r)
    return one(T) + b / λ * γij + a / λ * γij * γij
end


@inline function pₙ(n::S,tᵢ::T, tⱼ::T, λ::T, μ::T, ψ::T, r::T) where {S<:Integer, T<:AbstractFloat}
    n == 0 && return α(one(T), tᵢ, tⱼ, λ, μ, ψ, r)
    # return (one(T) - α(one(T), tᵢ, tⱼ, λ, μ, ψ, r)) * (one(T) - γ(one(T), tᵢ, tⱼ, λ, μ, ψ, r)) * γ(one(T), tᵢ, tⱼ, λ, μ, ψ, r)^(n-1)
    return β(one(T), tᵢ, tⱼ, λ, μ, ψ, r) * γ(one(T), tᵢ, tⱼ, λ, μ, ψ, r)^(n-1)
end


@inline function pₙ(n::Vector{S}, tᵢ::T, tⱼ::T, λ::T, μ::T, ψ::T, r::T) where {S<:Integer, T<:AbstractFloat}
    return [pₙ(nᵢ, tᵢ, tⱼ, λ, μ, ψ, r) for nᵢ in n]
end


@inline function pₙ(n::Vector{S}, tᵢ::T, tⱼ::Vector{T}, λ::T, μ::T, ψ::T, r::T) where {S<:Integer, T<:AbstractFloat}
    return transpose(reduce(hcat,[pₙ(n, tᵢ, tⱼᵢ, λ, μ, ψ, r) for tⱼᵢ in tⱼ]))
end