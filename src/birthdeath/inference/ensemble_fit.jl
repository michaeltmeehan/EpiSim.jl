function fit_bd_ensemble_mle(
    trees::Vector{Tree};
    fixed::Tuple{Symbol,Float64},
    r::Float64,
    ρ₀::Float64 = 0.0
)

    n = length(trees)

    λ̂ = Vector{Float64}(undef, n)
    μ̂ = Vector{Float64}(undef, n)
    ψ̂ = Vector{Float64}(undef, n)

    for i in eachindex(trees)

        λ, μ, ψ = mle_bd_constant(
            trees[i];
            fixed = fixed,
            r = r,
            ρ₀ = ρ₀
        )

        λ̂[i] = λ
        μ̂[i] = μ
        ψ̂[i] = ψ
    end

    return (λ = λ̂, μ = μ̂, ψ = ψ̂)
end