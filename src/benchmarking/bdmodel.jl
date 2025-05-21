import StatsBase: countmap

# Full probability distribution
function ξ(t::Float64; λ::Float64=birth_rate, μ::Float64=removal_rate)
    r = λ - μ
    r == 0. && return λ * t / (1. + λ * t)
    return μ * (exp(r * t) - 1.) / (λ * exp(r * t) - μ)
end


function η(t::Float64; λ::Float64=birth_rate, μ::Float64=removal_rate)
    return birth_rate / removal_rate * ξ(t, λ=λ, μ=μ)
end


# P(X(t) = n | X(0) = 1)
function p(n::Int, t::Float64; λ::Float64=birth_rate, μ::Float64=removal_rate)
    n == 0 && return ξ(t, λ=λ, μ=μ)
    return (1. - ξ(t, λ=λ, μ=μ)) * (1. - η(t, λ=λ, μ=μ)) * η(t, λ=λ, μ=μ)^(n-1)
end


function extinction_probability(t::Float64, n_0::Int; λ::Float64=birth_rate, μ::Float64=removal_rate)
    return ξ(t, λ=λ, μ=μ)^n_0
end

birth_rate = 0.5
death_rate = 0.1
sampling_rate = 0.1
initial_infected = 1

removal_rate = death_rate + sampling_rate

model = BDModel(birth_rate=birth_rate, 
                death_rate=death_rate, 
                sampling_rate=sampling_rate, 
                I=initial_infected, 
                agentic=false)

rng = Random.MersenneTwister(1234)
n_sim = 100_000
ens = simulate(rng, model, n_sim, stop_condition = (state) -> state.t > 5.)

t_grid = collect(0.:1.:5.)
empirical_prevalence = get_prevalence(ens, t_grid)
n_max = maximum(empirical_prevalence)
empirical_pop_distribution = get_distribution(empirical_prevalence, n_max)
theoretical_pop_distribution = reduce(vcat, [p.(0:n_max, t, λ=birth_rate, μ=removal_rate)' for t in t_grid])