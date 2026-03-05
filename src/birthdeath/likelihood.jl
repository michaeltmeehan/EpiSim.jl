using DifferentialEquations

@inline function E_constant(t, λ, μ, ψ; ρ₀ = zero(λ))

    t ≤ zero(λ) && return one(λ) - ρ₀

    δ = λ + μ + ψ
    disc = muladd(δ, δ, -4λ*μ)
    disc = ifelse(disc < zero(λ), zero(λ), disc)
    Δ = sqrt(disc)

    E_plus  = (δ + Δ) / (2λ)
    E_minus = (δ - Δ) / (2λ)

    E0 = one(λ) - ρ₀

    # C = (E0 - E_minus)/(E0 - E_plus)
    numC = E0 - E_minus
    denC = E0 - E_plus

    scale = max(abs(denC), one(λ))
    denC = abs(denC) < eps(λ)*scale ?
           copysign(eps(λ)*scale, denC) :
           denC

    C = numC / denC

    x = Δ * t

    # small Δ t branch
    if abs(x) ≤ sqrt(eps(λ))
        # exp(-x) ≈ 1 - x + x^2/2
        em = one(λ) - x + x*x/2
    else
        em = exp(-x)
    end

    num = E_minus - C * E_plus * em
    den = one(λ) - C * em

    scale = max(abs(den), one(λ))
    den = abs(den) < eps(λ)*scale ?
          copysign(eps(λ)*scale, den) :
          den

    return num / den
end

# @inline function E(t, T, λ, μ, ψ, r)
#     γ₀ = γ(zero(λ), t, T, λ, μ, ψ, r)
#     return one(λ) - ψ / λ * γ₀ / (one(λ) - γ₀)
# end

# This is the log of Φ(t) (eq 10. MacPherson et al. 2022)
@inline function g_constant(t, λ, μ, ψ; ρ₀ = zero(λ))

    t ≤ zero(λ) && return zero(λ)

    δ = λ + μ + ψ

    disc = muladd(δ, δ, -4λ*μ)
    disc = ifelse(disc < zero(λ), zero(λ), disc)
    Δ = sqrt(disc)

    # Riccati roots
    E_plus  = (δ + Δ) / (2λ)
    E_minus = (δ - Δ) / (2λ)

    E0 = one(λ) - ρ₀

    # C coefficient
    numC = E0 - E_minus
    denC = E0 - E_plus

    scale = max(abs(denC), one(λ))
    denC = abs(denC) < eps(λ)*scale ?
           copysign(eps(λ)*scale, denC) :
           denC

    C = numC / denC

    x = Δ * t

    em = abs(x) ≤ sqrt(eps(λ)) ?
         one(λ) - x + x*x/2 :
         exp(-x)

    num = one(λ) - C * em
    den = one(λ) - C

    scale = max(abs(den), one(λ))
    den = abs(den) < eps(λ)*scale ?
          copysign(eps(λ)*scale, den) :
          den

    return -x - 2 * log(num / den)
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

    Tfinal = maximum(tree.time)

    log_λ   = log(λ)
    log_ψ   = log(ψ)
    log_r   = log(r)
    log_1mr = log1p(-r)

    # survival conditioning
    E_T = E_constant(Tfinal, λ, μ, ψ; ρ₀=ρ₀)
    E_T = clamp(E_T, zero(E_T), one(E_T))
    ll = log1p(-E_T)

    for node in tree

        τ = Tfinal - node.time

        gτ = g_constant(τ, λ, μ, ψ; ρ₀=ρ₀)
        Eτ = E_constant(τ, λ, μ, ψ; ρ₀=ρ₀)
        Eτ = clamp(Eτ, zero(Eτ), one(Eτ))

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