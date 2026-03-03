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




# @inline function γ(w::Number, tᵢ::Number, tⱼ::Number, λ::Float64, μ::Float64, ψ::Float64, r::Float64)
#     aw = muladd(r * ψ, w, μ)
#     bw = muladd((one(w) - r) * ψ, w, -(λ + μ + ψ))

#     disc = bw*bw - 4 * aw * λ
#     disc = ifelse(disc < zero(disc), zero(disc), disc)
#     Δ    = sqrt(disc)

#     τ = tⱼ - tᵢ
#     x = Δ * τ

#     if abs(x) ≤ oftype(x, 1e-8)
#         denom = oftype(τ, 2) - (Δ + bw) * τ
#         e = eps(denom) * one(denom)
#         denom = ifelse(abs(denom) < e, copysign(e, denom), denom)
#         return (oftype(τ, 2) * λ * τ) / denom
#     else
#         num = oftype(τ, 2) * λ * (-expm1(-x))         # 1 - exp(-x)
#         den = (Δ - bw) + (Δ + bw) * exp(-x)
#         ed  = eps(den) * one(den)
#         den = ifelse(abs(den) < ed, copysign(ed, den), den)
#         return num / den
#     end
# end




# @inline function γ₀(tᵢ::Number, tⱼ::Number, λ::Float64, μ::Float64, ψ::Float64)
#     bw = -(λ + μ + ψ)

#     disc = bw*bw - 4 * μ * λ
#     disc = ifelse(disc < zero(disc), zero(disc), disc)
#     Δ    = sqrt(disc)

#     τ = tⱼ - tᵢ
#     x = Δ * τ

#     if abs(x) ≤ oftype(x, 1e-8)
#         denom = oftype(τ, 2) - (Δ + bw) * τ
#         e = eps(denom) * one(denom)
#         denom = ifelse(abs(denom) < e, copysign(e, denom), denom)
#         return (oftype(τ, 2) * λ * τ) / denom
#     else
#         num = oftype(τ, 2) * λ * (-expm1(-x))         # 1 - exp(-x)
#         den = (Δ - bw) + (Δ + bw) * exp(-x)
#         ed  = eps(den) * one(den)
#         den = ifelse(abs(den) < ed, copysign(ed, den), den)
#         return num / den
#     end
# end


# @inline function γ₁(tᵢ::Number, tⱼ::Number, λ::Float64, μ::Float64, ψ::Float64, r::Float64)
#     a = μ + r * ψ
#     b = (1. - r) * ψ - (λ + μ + ψ)

#     disc = b*b - 4 * a * λ
#     disc = ifelse(disc < zero(disc), zero(disc), disc)
#     Δ    = sqrt(disc)

#     τ = tⱼ - tᵢ
#     x = Δ * τ

#     if abs(x) ≤ oftype(x, 1e-8)
#         denom = oftype(τ, 2) - (Δ + b) * τ
#         e = eps(denom) * one(denom)
#         denom = ifelse(abs(denom) < e, copysign(e, denom), denom)
#         return (oftype(τ, 2) * λ * τ) / denom
#     else
#         num = oftype(τ, 2) * λ * (-expm1(-x))         # 1 - exp(-x)
#         den = (Δ - b) + (Δ + b) * exp(-x)
#         ed  = eps(den) * one(den)
#         den = ifelse(abs(den) < ed, copysign(ed, den), den)
#         return num / den
#     end
# end


# function p₀(tᵢ::Float64, tⱼ::Float64, λ::Float64, μ::Float64, ψ::Float64)
#     r = λ - μ - ψ
#     c₁ = sqrt(r^2 + 4 * λ * ψ)
#     c₂ = -r / c₁
#     τ = tⱼ - tᵢ
#     return ((λ + μ + ψ) + c₁ * (exp(-c₁ * τ) * (1. - c₂) - (1. + c₂)) / (exp(-c₁ * τ) * (1 - c₂) + (1 + c₂))) / (2 * λ)
# end


# function β(tᵢ::Number, tⱼ::Number, λ::Float64, μ::Float64, ψ::Float64, r::Float64)
#     a = μ + r * ψ
#     b = (1. - r) * ψ - (λ + μ + ψ)
#     γij = γ₁(tᵢ, tⱼ, λ, μ, ψ, r)
#     return (λ + b * γij + a * γij * γij) / λ
# end


# function p₀(tᵢ::Float64, tⱼ::Float64, λ::Float64, μ::Float64, ψ::Float64)
#     return 1. - ψ / λ * γ₀(tᵢ, tⱼ, λ, μ, ψ) / (1. - γ₀(tᵢ, tⱼ, λ, μ, ψ))
# end


# TODO: Add root conditioning
# @inline function likelihood(tree::Vector{Node}, λ::Float64, μ::Float64, ψ::Float64, r::Float64)
#     t_tip = tree[1].time
#     logΦ(t) = log(β(t, t_tip, λ, μ, ψ, r))
#     E(t) = p₀(t, t_tip, λ, μ, ψ)
#     logλ = log(λ)
#     logψ = log(ψ)
#     log1mr = log(1 - r)
#     logL = 0.0
#     for node in tree
#         if node isa Binary
#             logL += logλ + logΦ(node.time)
#         elseif node isa SampledLeaf
#             logL += logψ + log((1 - r) * E(node.time) + r) - logΦ(node.time)
#         elseif node isa SampledUnary
#             logL += logψ + log1mr
#         end
#     end
#     return logL
# end


# @inline function likelihood(tree::Tree, λ::Float64, μ::Float64, ψ::Float64, r::Float64)
#     t_tip = tree.time[1]
#     logΦ(t) = log(β(t, t_tip, λ, μ, ψ, r))
#     E(t) = p₀(t, t_tip, λ, μ, ψ)
#     logλ = log(λ)
#     logψ = log(ψ)
#     log1mr = log(1 - r)
#     logL = 0.0
#     for idx in eachindex(tree.time)
#         if tree.kind[idx] == K_Binary
#             logL += logλ + logΦ(tree.time[idx])
#         elseif tree.kind[idx] == K_SampledLeaf
#             logL += logψ + log((1 - r) * E(tree.time[idx]) + r) - logΦ(tree.time[idx])
#         elseif tree.kind[idx] == K_SampledUnary
#             logL += logψ + log1mr
#         end
#     end
#     return logL
# end