rng = Random.MersenneTwister(1234)

birth_rate = 0.5
death_rate = 0.1
sampling_rate = 0.1
initial_infected = 1

removal_rate = death_rate + sampling_rate

model = BirthDeathModel(birth_rate=birth_rate, 
                        death_rate=death_rate, 
                        sampling_rate=sampling_rate, 
                        I=initial_infected, 
                        agentic=false)

ens = simulate(rng, model, 10_000, stop_condition = (state) -> state.t > 5.)

theoretical_extinction_probability = calc_extinction_probability(model)
empirical_extinction_probability = get_extinction_probability(ens)

println("Theoretical extinction probability: ", theoretical_extinction_probability)
println("Empirical extinction probability  : ", empirical_extinction_probability)

tvec = collect(0.:0.1:5.)
empirical_prevalence = get_prevalence(ens, tvec)

empirical_mean_prevalence = mean(empirical_prevalence, dims=2)
theoretical_mean_prevalence = exp.((birth_rate - (removal_rate)) * tvec) * initial_infected

plot(tvec, empirical_mean_prevalence, label="Empirical Mean Prevalence", color=:blue)
plot!(tvec, theoretical_mean_prevalence, label="Theoretical Mean Prevalence", color=:red)

empirical_std_prevalence = [std(row) for row in eachrow(empirical_prevalence)]
theoretical_std_prevalence = initial_infected * (birth_rate + removal_rate) / (birth_rate - removal_rate) * exp.((birth_rate - removal_rate) * tvec) .* (exp.((birth_rate - removal_rate) * tvec) .- 1.)

plot(tvec, empirical_std_prevalence, label="Empirical Std Prevalence", color=:blue)
plot!(tvec, theoretical_std_prevalence, label="Theoretical Std Prevalence", color=:red)

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


function simulate_birth_death_process(rng::AbstractRNG,
                                      λ::Float64,
                                      μ::Float64,
                                      n_0::Int;
                                      t_max::Float64=5.0)

    n = n_0
    t = 0.0
    tvec = [t]
    nvec = [n]
    while t < t_max && n > 0
        # Calculate the time to the next event
        dt = -log(rand(rng)) / (λ * n + μ * n)
        t += dt

        # Determine the type of event (birth or death)
        if rand(rng) < λ * n / (λ * n + μ * n)
            n += 1  # Birth
        else
            n -= 1  # Death
        end

        push!(tvec, t)
        push!(nvec, n)
    end
    return tvec, nvec
end


n_sim = 100_000
t_grid = collect(0.:0.1:5.)
n = fill(0, length(t_grid), n_sim)

for j in 1:n_sim
    tvec, nvec = simulate_birth_death_process(rng, birth_rate, removal_rate, initial_infected)
    for i in eachindex(t_grid)
        idx = searchsortedlast(tvec, t_grid[i])
        n[i, j] = nvec[idx]
    end
end

n_max = maximum(n)
n_dist = zeros(length(t_grid), n_max + 1)
for i in eachindex(t_grid)
    freqs = countmap(n[i, :])
    for (k, v) in freqs
        n_dist[i, k + 1] = v / n_sim
    end
end


function get_distribution(a::Matrix{T}, max::Int) where T <: Real
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


initial_infected = 1
model = BirthDeathModel(birth_rate=birth_rate, 
                        death_rate=death_rate, 
                        sampling_rate=sampling_rate, 
                        I=initial_infected, 
                        agentic=true) # Adjust to false to check the aggregate model
ens = simulate(rng, model, 10_000, stop_condition = (state) -> state.t > 5.)
empirical_prevalence = get_prevalence(ens, t_grid)
n_max = maximum(empirical_prevalence)
pop_distribution = get_distribution(empirical_prevalence, n_max)
[p.(0:5, t_grid[1], λ=birth_rate, μ=removal_rate)';
p.(0:5, t_grid[2], λ=birth_rate, μ=removal_rate)';
p.(0:5, t_grid[3], λ=birth_rate, μ=removal_rate)';
p.(0:5, t_grid[4], λ=birth_rate, μ=removal_rate)';
p.(0:5, t_grid[5], λ=birth_rate, μ=removal_rate)';
p.(0:5, t_grid[end], λ=birth_rate, μ=removal_rate)']


# More than 1 initial infected
n_sim = 100_000
initial_infected = 3
model = BirthDeathModel(birth_rate=birth_rate, 
                        death_rate=death_rate, 
                        sampling_rate=sampling_rate, 
                        I=initial_infected, 
                        agentic=false)
ens = simulate(rng, model, n_sim, stop_condition = (state) -> state.t > 5.)
empirical_prevalence = get_prevalence(ens, t_grid)
extinction_prob = zeros(length(t_grid))
for i in eachindex(t_grid)
    count = 0
    for j in 1:n_sim
        count += empirical_prevalence[i, j] == 0 ? 1 : 0
    end
    extinction_prob[i] = count / n_sim 
end

[extinction_prob ξ.(t_grid).^initial_infected]
