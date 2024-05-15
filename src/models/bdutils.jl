function initialize_linelist(N::Int64; type::Union{Nothing, Vector{Int64}}=nothing)::DataFrame
    type = isnothing(type) ? fill(1,N) : type
    return DataFrame(
        event = fill(1, N),
        child_id = collect(1:N),
        child_type = type,
        parent_id = fill(0, N),
        parent_type = fill(0, N),
        t_birth = fill(0., N),
        t_death = fill(Inf, N),
        t_sam = fill(-1., N))
end


function initialize_active(N::Vector{Int64})::Tuple{Vector{Vector{Int64}}, Vector{Int64}}
    N_tot = 0
    active = [Vector{Int64}() for _ in 1:length(N)]
    type = Vector{Int64}()
    for t in eachindex(N)
        for _ in 1:N[t]
            N_tot += 1
            push!(active[t], N_tot)
            push!(type, t)
        end
    end
    return active, type
end


function check_parameters(params::CRBDParameters)
    @unpack N₀, λ, μ, ψ, ρ₀, r, t_max = params
    # Ensure input parameters are valid
    @assert N₀ ≥ 0 "Initial population size (N₀) must be non-negative"
    @assert λ ≥ 0. "Birth rate (λ) must be non-negative"
    @assert μ ≥ 0. "Death rate (μ) must be non-negative"
    @assert ψ ≥ 0. "Extinct sampling rate (ψ) must be non-negative"
    @assert 0. ≤ ρ₀ ≤ 1. "Extant sampling rate (ρ₀) must be between 0 and 1"
    @assert 0. ≤ r ≤ 1. "Removal probability (r) must be between 0 and 1"
    @assert t_max > 0 "Maximum simulation time (t_max) must be positive"
end


function check_parameters(params::MTBDParameters)
    @unpack n_types, N₀, λ, μ, γ, ψ, ρ₀, r, t_max = params
    @assert size(λ)[1] == size(λ)[2] == length(μ) == size(γ)[1] == size(γ)[2] == length(ψ) == length(ρ₀) == length(r) "Parameters must have consistent dimensions (=n_types)"
    @assert all(N₀ .≥ 0) "Initial population size (N₀) must be non-negative"
    @assert all(λ .≥ 0.) "Birth rate (λ) must be non-negative"
    @assert all(μ .≥ 0.) "Death rate (μ) must be non-negative"
    @assert all(ψ .≥ 0.) "Extinct sampling rate (ψ) must be non-negative"
    @assert all(0. .≤ ρ₀ .≤ 1.) "Extant sampling rate (ρ₀) must be between 0 and 1"
    @assert all(0. .≤ r .≤ 1.) "Removal probability (r) must be between 0 and 1"
    @assert t_max > 0 "Maximum simulation time (t_max) must be positive"
end