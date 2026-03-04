using DifferentialEquations

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
    ρ0::Real = 0.0,
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
    E0 = 1 - ρ0
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


λ(t) = 2.0
μ(t) = 0.5
ψ(t) = 0.5

sol = solve_Eg(λ, μ, ψ; ρ0=0.0, tspan=(0.0, 3.0))

E1 = sol(1.0)[1]
g1 = sol(1.0)[2]
Φ1 = exp(g1)

println("E(1) = ", E1)
println("Φ(1) = ", Φ1)