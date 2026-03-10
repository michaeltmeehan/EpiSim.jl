using DifferentialEquations


@inline safe_log(x) =
    x > zero(x) ? log(x) : -1e12

@inline clamp01(x) =
    clamp(x, eps(x), one(x) - eps(x))

@inline function expm1_safe(x)

    abs(x) < sqrt(eps(x)) ?
        x + x*x/2 :
        exp(x) - 1

end


@inline function bd_coefficients(λ, μ, ψ)

    δ = λ + μ + ψ

    disc = muladd(δ, δ, -4λ*μ)
    disc = ifelse(disc < zero(λ), zero(λ), disc)

    Δ = sqrt(disc)

    E_plus  = (δ + Δ) / (2λ)
    E_minus = (δ - Δ) / (2λ)

    return Δ, E_plus, E_minus
end


@inline function E_constant(t, λ, μ, ψ; ρ₀ = zero(λ))

    t ≤ zero(λ) && return one(λ) - ρ₀

    Δ, E_plus, E_minus = bd_coefficients(λ, μ, ψ)

    E0 = one(λ) - ρ₀

    numC = E0 - E_minus
    denC = E0 - E_plus

    if abs(denC) < eps(λ)
        return E_minus
    end

    C = numC / denC

    x = Δ * t
    em = exp(-x)

    num = E_minus - C * E_plus * em
    den = one(λ) - C * em

    if abs(den) < eps(λ)
        return E_minus
    end

    E = num / den

    return clamp01(E)
end

@inline function E_constant_coeff(t, Δ, E_plus, E_minus; ρ₀)

    t ≤ zero(t) && return one(t) - ρ₀

    E0 = one(t) - ρ₀

    numC = E0 - E_minus
    denC = E0 - E_plus

    if abs(denC) < eps(t)
        return E_minus
    end

    C = numC / denC

    x  = Δ * t
    em = exp(-x)

    num = E_minus - C * E_plus * em
    den = one(t) - C * em

    if abs(den) < eps(t)
        return E_minus
    end

    E = num / den

    return clamp01(E)
end

# This is the log of Φ(t) (eq 10. MacPherson et al. 2022)
@inline function g_constant(t, λ, μ, ψ; ρ₀ = zero(λ))

    t ≤ zero(λ) && return zero(λ)

    Δ, E_plus, E_minus = bd_coefficients(λ, μ, ψ)

    E0 = one(λ) - ρ₀

    numC = E0 - E_minus
    denC = E0 - E_plus

    if abs(denC) < eps(λ)
        return -Δ * t
    end

    C = numC / denC

    x = Δ * t
    em = exp(-x)

    num = one(λ) - C * em
    den = one(λ) - C

    if abs(num) < eps(λ) || abs(den) < eps(λ)
        return -Δ * t
    end

    log_ratio = log(num) - log(den)

    return -x - 2 * log_ratio
end

@inline function g_constant_coeff(t, Δ, E_plus, E_minus; ρ₀)

    t ≤ zero(t) && return zero(t)

    E0 = one(t) - ρ₀

    numC = E0 - E_minus
    denC = E0 - E_plus

    if abs(denC) < eps(t)
        return -Δ * t
    end

    C = numC / denC

    x  = Δ * t
    em = exp(-x)

    num = one(t) - C * em
    den = one(t) - C

    if abs(num) < eps(t) || abs(den) < eps(t)
        return -Δ * t
    end

    log_ratio = log(num) - log(den)

    return -x - 2 * log_ratio
end


"""
Solve the coupled backward-time ODE system for:

E(t)  = extinction probability
g(t)  = log Φ(t)

Time convention:
t = 0   present
t > 0   into the past
"""
function solve_Eg(
    λ::Function,
    μ::Function,
    ψ::Function;
    ρ₀::Real = 0.0,
    tspan::Tuple{<:Real,<:Real} = (0.0, 5.0),
    abstol::Real = 1e-12,
    reltol::Real = 1e-12
)

    # ODE system
    function Eg_ode!(du, u, p, t)
        E = u[1]
        g = u[2]

        λt = λ(t)
        μt = μ(t)
        ψt = ψ(t)

        # dE/dt
        du[1] =
            λt * E^2 -
            (λt + μt + ψt) * E +
            μt

        # dg/dt
        du[2] =
            2 * λt * E -
            (λt + μt + ψt)
    end

    # Initial conditions at present
    E0 = 1 - ρ₀
    g0 = 0.0

    u0 = [E0, g0]

    prob = ODEProblem(Eg_ode!, u0, tspan, nothing)

    sol = solve(
        prob,
        Vern9(),
        abstol = abstol,
        reltol = reltol
    )

    return sol
end


@inline backward_time(t, Tfinal) = Tfinal - t

@inline logaddexp(a, b) =
    max(a,b) + log1p(exp(-abs(a-b)))


function bd_loglikelihood_constant(tree::Tree,
                                   λ, μ, ψ, r;
                                   ρ₀ = zero(λ))
    try
        Tfinal = maximum(tree.time)

        log_λ   = safe_log(λ)
        log_ψ   = safe_log(ψ)
        log_r   = safe_log(r)
        log_1mr = log1p(-r)

        Δ, E_plus, E_minus = bd_coefficients(λ, μ, ψ)

        # survival conditioning
        E_T = clamp01(E_constant_coeff(Tfinal, Δ, E_plus, E_minus; ρ₀=ρ₀))
        E_T ≥ 1 - sqrt(eps(E_T)) && return -1e12

        ll = log1p(-E_T)

        for node in tree

            τ = Tfinal - node.time

            gτ = g_constant_coeff(τ, Δ, E_plus, E_minus; ρ₀=ρ₀)
            Eτ = clamp01(E_constant_coeff(τ, Δ, E_plus, E_minus; ρ₀=ρ₀))

            if node.kind === Binary

                ll += log_λ + gτ

            elseif node.kind === SampledLeaf

                # stable log(r + (1-r)Eτ)
                log_term = logaddexp(log_r,
                                     log_1mr + log(Eτ))

                ll += log_ψ + log_term - gτ

            elseif node.kind === SampledUnary

                ll += log_ψ + log_1mr
            end
        end

        return ll
    catch
        return -1e12
    end
end


function bd_loglikelihood_ode(
    tree::Tree,
    λ::Function,
    μ::Function,
    ψ::Function,
    r::T;
    ρ₀::T = zero(T),
    abstol::Real = 1e-12,
    reltol::Real = 1e-12
) where {T<:AbstractFloat}

    Tfinal = maximum(tree.time)

    # Solve ODE once from 0 → Tfinal
    sol = solve_Eg(
        λ, μ, ψ;
        ρ₀ = ρ₀,
        tspan = (zero(T), Tfinal),
        abstol = abstol,
        reltol = reltol
    )

    log_r   = log(r)
    log_1mr = log1p(-r)

    # Survival conditioning
    E_T = sol(Tfinal)[1]
    ll  = log1p(-E_T)

    for node in tree

        τ = Tfinal - node.time

        u = sol(τ)
        Eτ = u[1]
        gτ = u[2]

        if node.kind === Binary

            ll += log(λ(τ)) + gτ

        elseif node.kind === SampledLeaf

            # stable log(r + (1-r)Eτ)
            log_term = logaddexp(log_r,
                                 log_1mr + log(Eτ))

            ll += log(ψ(τ)) + log_term - gτ

        elseif node.kind === SampledUnary

            ll += log(ψ(τ)) + log_1mr
        end
    end

    return ll
end