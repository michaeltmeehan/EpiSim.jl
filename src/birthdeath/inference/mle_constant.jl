using Optim
using ForwardDiff
using LinearAlgebra
using ADTypes: AutoForwardDiff

# ============================================================
# 1. Fixed-Parameter Specification (Identifiability Handling)
# ============================================================

struct BDFixedSpec{T}
    fixed_symbol::Symbol   # :λ, :μ, or :ψ
    fixed_value::T
end

function expand_rates(
    θ::AbstractVector{T},
    spec::BDFixedSpec
) where {T<:Real}

    if spec.fixed_symbol === :λ
        λ = spec.fixed_value
        μ = exp(θ[1])
        ψ = exp(θ[2])

    elseif spec.fixed_symbol === :μ
        λ = exp(θ[1])
        μ = spec.fixed_value
        ψ = exp(θ[2])

    elseif spec.fixed_symbol === :ψ
        λ = exp(θ[1])
        μ = exp(θ[2])
        ψ = spec.fixed_value

    else
        error("Invalid fixed parameter. Must be :λ, :μ, or :ψ")
    end

    return λ, μ, ψ
end


# ============================================
# 2. Negative Log-Likelihood (AD Compatible)
# ============================================

function bd_negloglikelihood_constant(
    θ::AbstractVector{T},
    tree::Tree,
    spec::BDFixedSpec;
    r,
    ρ₀ = zero(eltype(θ))
) where {T}

    λ, μ, ψ = expand_rates(θ, spec)

    ll = bd_loglikelihood_constant(tree, λ, μ, ψ, r; ρ₀=ρ₀)

    if !isfinite(ll)
        return 1e12   # large penalty
    end

    return -ll
end


# ============================================
# 3. Objective Struct (Avoid Closure Issues)
# ============================================

struct BDConstantObjective{T}
    tree::Tree
    spec::BDFixedSpec{T}
    r::T
    ρ₀::T
end

function (obj::BDConstantObjective)(θ::AbstractVector{T}) where {T<:Real}
    return bd_negloglikelihood_constant(
        θ,
        obj.tree,
        obj.spec;
        r = obj.r,
        ρ₀ = obj.ρ₀
    )
end


# ============================================
# 4. Main MLE Routine
# ============================================

function fit_bd_constant(
    tree::Tree;
    fixed::Tuple{Symbol,Float64},
    r::Float64,
    ρ₀::Float64 = 0.0,
    θ_init = log.([0.5, 0.5])
)

    spec = BDFixedSpec(fixed[1], fixed[2])
    obj  = BDConstantObjective(tree, spec, r, ρ₀)

    # Optimisation
    lower = log.([1e-6, 1e-6])  # avoid zero rates
    upper = log.([1e3, 1e3])    # reasonable upper

    result = optimize(
        obj,
        lower,
        upper,
        θ_init,
        Fminbox(LBFGS());
        autodiff = AutoForwardDiff()
    )

    θ̂ = Optim.minimizer(result)

    # Gradient and Hessian at optimum
    grad = ForwardDiff.gradient(obj, θ̂)
    H    = ForwardDiff.hessian(obj, θ̂)

    # Invert observed information
    vcov = inv(H)

    se_log = sqrt.(diag(vcov))

    # Recover rates
    λ̂, μ̂, ψ̂ = expand_rates(θ̂, spec)

    # Delta-method SEs on rate scale
    if spec.fixed_symbol === :λ
        se_rate = (
            λ = 0.0,
            μ = μ̂ * se_log[1],
            ψ = ψ̂ * se_log[2]
        )

    elseif spec.fixed_symbol === :μ
        se_rate = (
            λ = λ̂ * se_log[1],
            μ = 0.0,
            ψ = ψ̂ * se_log[2]
        )

    else  # fixed ψ
        se_rate = (
            λ = λ̂ * se_log[1],
            μ = μ̂ * se_log[2],
            ψ = 0.0
        )
    end

    return (
        result   = result,
        θ̂        = θ̂,
        gradient = grad,
        hessian  = H,
        vcov     = vcov,
        se_log   = se_log,
        rates    = (λ = λ̂, μ = μ̂, ψ = ψ̂),
        se_rate  = se_rate
    )
end


function mle_bd_constant(
    tree::Tree;
    fixed::Tuple{Symbol,Float64},
    r::Float64,
    ρ₀::Float64 = 0.0,
    θ_init = log.([1.0, 0.5])
)

    spec = BDFixedSpec(fixed[1], fixed[2])
    obj  = BDConstantObjective(tree, spec, r, ρ₀)

    lower = log.([1e-6, 1e-6])  # avoid zero rates
    upper = log.([1e3, 1e3])    # reasonable upper

    result = optimize(
        obj,
        lower,
        upper,
        θ_init,
        Fminbox(LBFGS());
        autodiff = AutoForwardDiff()
    )

    θ̂ = Optim.minimizer(result)

    return expand_rates(θ̂, spec)
end