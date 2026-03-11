using Optim
using ForwardDiff
using LinearAlgebra
using ADTypes: AutoForwardDiff

function tuple_to_vector(nt::NamedTuple)
    collect(values(nt))
end

# ============================================================
# 1. Fixed-Parameter Specification (Identifiability Handling)
# ============================================================

struct BDFixedSpec{T}
    fixed_symbol::Symbol   # :λ, :μ, or :ψ
    fixed_value::T
end

abstract type AbstractBDParameterization end

struct RateParameterization{T} <: AbstractBDParameterization
    spec::BDFixedSpec{T}
end

struct R0DeltaParameterization{T} <: AbstractBDParameterization
    spec::BDFixedSpec{T}
end

backtransform(::RateParameterization, λ, μ, ψ) =
    (λ = λ, μ = μ, ψ = ψ)

function backtransform(::R0DeltaParameterization, λ, μ, ψ)

    δ = μ + ψ
    p = ψ / δ
    R₀ = λ / δ

    return (R₀ = R₀, δ = δ, p = p)
end

function parameter_transform(param, θ)

    λ, μ, ψ = expand_rates(param, θ)

    tuple_to_vector(backtransform(param, λ, μ, ψ))
end

function expand_rates(
    param::RateParameterization,
    θ::AbstractVector{T}
    ) where {T<:Real}

    spec = param.spec

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

# TODO: Expand to include r and ρ₀
function expand_rates(
    param::R0DeltaParameterization, 
    θ::AbstractVector{T}
    ) where {T<:Real}

    spec = param.spec

    if spec.fixed_symbol === :R0
        R₀ = spec.fixed_value
        δ  = exp(θ[1])
        p  = 1/(1+exp(-θ[2]))

    elseif spec.fixed_symbol === :δ
        R₀ = exp(θ[1])
        δ  = spec.fixed_value
        p  = 1/(1+exp(-θ[2]))

    elseif spec.fixed_symbol === :p
        R₀ = exp(θ[1])
        δ  = exp(θ[2])
        p  = spec.fixed_value

    else
        error("Invalid fixed parameter")
    end

    λ = R₀ * δ
    ψ = p * δ
    μ = (1 - p) * δ

    return λ, μ, ψ
end


# ============================================
# 2. Negative Log-Likelihood (AD Compatible)
# ============================================

function bd_negloglikelihood(
    θ::AbstractVector{T},
    tree::Tree,
    param::AbstractBDParameterization;
    r,
    ρ₀ = zero(eltype(θ))
) where {T}

    λ, μ, ψ = expand_rates(param, θ)

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
    param::AbstractBDParameterization
    r::T
    ρ₀::T
end

function (obj::BDConstantObjective)(θ::AbstractVector{T}) where {T<:Real}
    return bd_negloglikelihood(
        θ,
        obj.tree,
        obj.param;
        r = obj.r,
        ρ₀ = obj.ρ₀
    )
end


# ============================================
# 4. Main MLE Routine
# ============================================

function fit_bd_full(
    tree::Tree;
    param::AbstractBDParameterization,
    r::Float64,
    ρ₀::Float64 = 0.0,
    θ_init = zeros(2)
)

    obj  = BDConstantObjective(tree, param, r, ρ₀)

    lower = log.([1e-6, 1e-6])
    upper = log.([1e3, 1e3])

    result = optimize(
        obj,
        lower,
        upper,
        θ_init,
        Fminbox(LBFGS());
        autodiff = AutoForwardDiff()
    )

    θ̂ = Optim.minimizer(result)

    # derivatives in θ-space
    grad = ForwardDiff.gradient(obj, θ̂)
    H    = ForwardDiff.hessian(obj, θ̂)

    vcov_θ = inv(H)

    se_θ = sqrt.(diag(vcov_θ))

    # recover rates
    λ̂, μ̂, ψ̂ = expand_rates(param, θ̂)

    rates = (λ = λ̂, μ = μ̂, ψ = ψ̂)

    # transform parameters
    params = backtransform(param, λ̂, μ̂, ψ̂)

    # delta method for transformed parameters
    J = ForwardDiff.jacobian(θ -> parameter_transform(param, θ), θ̂)

    vcov_param = J * vcov_θ * J'
    se_param = sqrt.(diag(vcov_param))

    return (
        result = result,
        θ̂ = θ̂,
        gradient = grad,
        hessian = H,
        vcov_θ = vcov_θ,
        se_θ = se_θ,
        rates = rates,
        parameters = params,
        vcov_parameters = vcov_param,
        se_parameters = se_param
    )
end


function fit_bd_pars(
    tree::Tree;
    param::AbstractBDParameterization,
    r::Float64,
    ρ₀::Float64 = 0.0,
    θ_init = zeros(2)
)

    fit = fit_bd_full(
        tree;
        param = param,
        r = r,
        ρ₀ = ρ₀,
        θ_init = θ_init
    )

    return fit.parameters
end