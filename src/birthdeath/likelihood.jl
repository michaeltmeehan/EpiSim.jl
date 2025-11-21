@inline function γ(w::Number, tᵢ::Number, tⱼ::Number, λ::Float64, μ::Float64, ψ::Float64, r::Float64)
    aw = muladd(r * ψ, w, μ)
    bw = muladd((one(w) - r) * ψ, w, -(λ + μ + ψ))

    disc = bw*bw - 4 * aw * λ
    disc = ifelse(disc < zero(disc), zero(disc), disc)
    Δ    = sqrt(disc)

    τ = tⱼ - tᵢ
    x = Δ * τ

    if abs(x) ≤ oftype(x, 1e-8)
        denom = oftype(τ, 2) - (Δ + bw) * τ
        e = eps(denom) * one(denom)
        denom = ifelse(abs(denom) < e, copysign(e, denom), denom)
        return (oftype(τ, 2) * λ * τ) / denom
    else
        num = oftype(τ, 2) * λ * (-expm1(-x))         # 1 - exp(-x)
        den = (Δ - bw) + (Δ + bw) * exp(-x)
        ed  = eps(den) * one(den)
        den = ifelse(abs(den) < ed, copysign(ed, den), den)
        return num / den
    end
end


@inline function γ₀(tᵢ::Number, tⱼ::Number, λ::Float64, μ::Float64, ψ::Float64)
    bw = -(λ + μ + ψ)

    disc = bw*bw - 4 * μ * λ
    disc = ifelse(disc < zero(disc), zero(disc), disc)
    Δ    = sqrt(disc)

    τ = tⱼ - tᵢ
    x = Δ * τ

    if abs(x) ≤ oftype(x, 1e-8)
        denom = oftype(τ, 2) - (Δ + bw) * τ
        e = eps(denom) * one(denom)
        denom = ifelse(abs(denom) < e, copysign(e, denom), denom)
        return (oftype(τ, 2) * λ * τ) / denom
    else
        num = oftype(τ, 2) * λ * (-expm1(-x))         # 1 - exp(-x)
        den = (Δ - bw) + (Δ + bw) * exp(-x)
        ed  = eps(den) * one(den)
        den = ifelse(abs(den) < ed, copysign(ed, den), den)
        return num / den
    end
end


@inline function γ₁(tᵢ::Number, tⱼ::Number, λ::Float64, μ::Float64, ψ::Float64, r::Float64)
    a = μ + r * ψ
    b = (1. - r) * ψ - (λ + μ + ψ)

    disc = b*b - 4 * a * λ
    disc = ifelse(disc < zero(disc), zero(disc), disc)
    Δ    = sqrt(disc)

    τ = tⱼ - tᵢ
    x = Δ * τ

    if abs(x) ≤ oftype(x, 1e-8)
        denom = oftype(τ, 2) - (Δ + b) * τ
        e = eps(denom) * one(denom)
        denom = ifelse(abs(denom) < e, copysign(e, denom), denom)
        return (oftype(τ, 2) * λ * τ) / denom
    else
        num = oftype(τ, 2) * λ * (-expm1(-x))         # 1 - exp(-x)
        den = (Δ - b) + (Δ + b) * exp(-x)
        ed  = eps(den) * one(den)
        den = ifelse(abs(den) < ed, copysign(ed, den), den)
        return num / den
    end
end


# function p₀(tᵢ::Float64, tⱼ::Float64, λ::Float64, μ::Float64, ψ::Float64)
#     r = λ - μ - ψ
#     c₁ = sqrt(r^2 + 4 * λ * ψ)
#     c₂ = -r / c₁
#     τ = tⱼ - tᵢ
#     return ((λ + μ + ψ) + c₁ * (exp(-c₁ * τ) * (1. - c₂) - (1. + c₂)) / (exp(-c₁ * τ) * (1 - c₂) + (1 + c₂))) / (2 * λ)
# end


function β(tᵢ::Number, tⱼ::Number, λ::Float64, μ::Float64, ψ::Float64, r::Float64)
    a = μ + r * ψ
    b = (1. - r) * ψ - (λ + μ + ψ)
    γij = γ₁(tᵢ, tⱼ, λ, μ, ψ, r)
    return (λ + b * γij + a * γij * γij) / λ
end


function p₀(tᵢ::Float64, tⱼ::Float64, λ::Float64, μ::Float64, ψ::Float64)
    return 1. - ψ / λ * γ₀(tᵢ, tⱼ, λ, μ, ψ) / (1. - γ₀(tᵢ, tⱼ, λ, μ, ψ))
end


# TODO: Add root conditioning
# TODO: Fix tree structure (transmission events do not necessarily correspond to binary nodes in the sampled tree)
@inline function likelihood(tree::Vector{Node}, λ::Float64, μ::Float64, ψ::Float64, r::Float64)
    t_tip = tree[1].time
    logΦ(t) = log(β(t, t_tip, λ, μ, ψ, r))
    E(t) = p₀(t, t_tip, λ, μ, ψ)
    logλ = log(λ)
    logψ = log(ψ)
    log1mr = log(1 - r)
    logL = 0.0
    for node in tree
        if node isa Binary
            logL += logλ + logΦ(node.time)
        elseif node isa SampledLeaf
            logL += logψ + log((1 - r) * E(node.time) + r) - logΦ(node.time)
        elseif node isa SampledUnary
            logL += logψ + log1mr
        end
    end
    return logL
end