using Distributions
using Roots

############################################################
# Mean infectious duration for Erlang competition
# E[min(T_r, T_s)] where T_r ~ Erlang(k, λr), T_s ~ Erlang(k, λs)
#
# Derived from ∫₀^∞ S_r(t)S_s(t)dt, expanding the incomplete
# gamma survival functions and integrating term-by-term:
#
#   E[T] = Σᵢ₌₀^{k-1} Σⱼ₌₀^{k-1} C(i+j, i) / (λr + λs)^{i+j+1}
#
# Note: the individual rates λr, λs cancel out completely.
############################################################

function mean_infectious_duration(λ, k)

    total = 0.0

    for i in 0:k-1
        for j in 0:k-1
            total += binomial(i+j, i) / λ^(i+j+1)
        end
    end

    return total
end

############################################################
# Sampling probability
# P(T_s < T_r) where both are Erlang(k)
# Uses the negative binomial CDF representation
############################################################

function sampling_probability(λr, λs, k)

    q = λs / (λr + λs)
    total = 0.0

    for i in 0:k-1
        total += binomial(k+i-1, i) * q^k * (1-q)^i
    end

    return total
end

############################################################
# Solve for λr, λs given:
#   δ  - removal rate (mean infectious period = 1/δ)
#   p  - sampling proportion P(T_s < T_r)
#   k  - shared Erlang shape parameter
#
# Strategy:
#   1. Solve for q = λs/(λr+λs) from sampling_probability = p
#      (this depends only on the ratio of rates, not their scale)
#   2. Solve for total rate λ = λr+λs using mean duration = 1/δ
#      (the mean depends only on λ, not the individual rates)
#   3. Split: λs = q*λ, λr = (1-q)*λ
############################################################

function solve_removal_rates(δ, p, k)

    # Step 1: solve for q = λs/(λr+λs)
    function f(q)
        total = 0.0
        for i in 0:k-1
            total += binomial(k+i-1, i) * q^k * (1-q)^i
        end
        return total - p
    end

    q = find_zero(f, (1e-10, 1-1e-10))

    # Step 2: solve for total rate λ = λr+λs directly
    h(λ) = mean_infectious_duration(λ, k) - 1/δ
    λ_total = find_zero(h, (1e-6, 1e6))

    # Step 3: split by q
    λs = q * λ_total
    λr = (1-q) * λ_total

    return λr, λs
end

############################################################
# Transmission rate distribution
# If R0 = E[β]/δ and β ~ LogNormal, parameterise accordingly
# σβ is the log-scale standard deviation of β
############################################################

function transmission_lognormal(R0, δ, σβ)
    meanβ = R0 * δ
    μβ = log(meanβ) - σβ^2 / 2
    return LogNormal(μβ, σβ)
end

############################################################
# Main parameter mapping
# Returns Erlang rate parameters and distributions for a
# two-pathway removal model (recovery + sampling)
############################################################

function epidemic_parameter_mapping(R0, δ, p, k, σβ)

    λr, λs = solve_removal_rates(δ, p, k)

    recovery_distribution     = Gamma(k, 1/λr)   # Erlang(k, λr)
    sampling_distribution     = Gamma(k, 1/λs)   # Erlang(k, λs)
    transmission_distribution = transmission_lognormal(R0, δ, σβ)

    return (
        λr = λr,
        λs = λs,
        recovery_distribution     = recovery_distribution,
        sampling_distribution     = sampling_distribution,
        transmission_distribution = transmission_distribution
    )
end