@inline function popr!(rng::AbstractRNG, v::Vector)
    isempty(v) && throw(ArgumentError("Cannot pop from an empty vector."))
    i = rand(rng, 1:length(v))
    v[i], v[end] = v[end], v[i]   # swap with last element
    return pop!(v)
end


"""
    wsample(rng, x, w) -> el
    wsample(x, w)      -> el               # uses default RNG
    wsampleindex(rng, w) -> i::Int
    wsampleindex(w)      -> i::Int

Weighted sampling without replacement-of-mass: picks one element of `x` with
probabilities proportional to weights `w`. `x` and `w` must have same length.

- `w` may be any `AbstractVector{<:Real}`; values must be finite and ≥ 0.
- Weights are *not* normalized; scaling `w` by a positive constant is OK.
- Throws `ArgumentError` if `w` is empty, sums to 0, or contains invalid values.
"""
@inline function wsample(rng::AbstractRNG,
                        x::AbstractVector{T},
                        w::AbstractVector{<:Real})::T where {T}
    i = wsampleindex(rng, w)
    @inbounds return x[i]
end

@inline wsample(x::AbstractVector{T}, w::AbstractVector{<:Real}) where {T} =
    wsample(Random.default_rng(), x, w)

@inline function wsampleindex(rng::AbstractRNG,
                             w::AbstractVector{<:Real})::Int
    n = length(w)
    n == 0 && throw(ArgumentError("weights must be non-empty"))
    # Validate + compute total (single pass)
    tot = zero(Float64)
    @inbounds @simd for i in 1:n
        wi = float(w[i])
        (!isfinite(wi) || wi < 0.0) &&
            throw(ArgumentError("weights must be finite and ≥ 0 (bad at i=$i: $wi)"))
        tot += wi
    end
    tot <= 0.0 && throw(ArgumentError("sum(weights) must be > 0"))

    # Draw and scan (no allocations)
    r = rand(rng) * tot
    acc = 0.0
    @inbounds @simd for i in 1:n
        acc += float(w[i])
        if r <= acc
            return i
        end
    end
    # Fallback for possible FP rounding (return last index)
    return n
end

@inline wsampleindex(w::AbstractVector{<:Real})::Int =
    wsampleindex(Random.default_rng(), w)


@inline function wsampleindex(rng::AbstractRNG,
                              W::AbstractArray{<:Real})
    isempty(W) && throw(ArgumentError("weights must be non-empty"))

    # 1) Validate + sum in one pass (linear indexing for SIMD)
    tot = 0.0
    @inbounds @simd for I in eachindex(W)
        wi = float(W[I])
        (!isfinite(wi) || wi < 0.0) &&
            throw(ArgumentError("weights must be finite and ≥ 0 (bad at linidx=$I: $wi)"))
        tot += wi
    end
    tot <= 0.0 && throw(ArgumentError("sum(weights) must be > 0"))

    # 2) Draw and single scan (no allocations)
    r = rand(rng) * tot
    acc = 0.0
    @inbounds @simd for I in eachindex(W)
        acc += float(W[I])
        if r <= acc
            # Map linear index -> CartesianIndex in O(1)
            return CartesianIndices(W)[I]
        end
    end
    # FP fallback
    return CartesianIndices(W)[lastindex(W)]
end

@inline wsampleindex(W::AbstractArray{<:Real}) =
    wsampleindex(Random.default_rng(), W)

# 2-D convenience returning (i, j)
@inline function wsampleindex2(rng::AbstractRNG,
                               W::AbstractMatrix{<:Real})::Tuple{Int,Int}
    CI = wsampleindex(rng, W)
    return (CI.I[1], CI.I[2])
end

@inline wsampleindex2(W::AbstractMatrix{<:Real}) =
    wsampleindex2(Random.default_rng(), W)


# Draw (i,j) from weights proportional to w[i,j] * I[j], without forming W .* I'
@inline function wsampleindex_cols(rng::AbstractRNG,
                                   W::AbstractMatrix{<:Real},
                                   I::AbstractVector{<:Real})::Tuple{Int,Int}
    nrow, ncol = size(W)
    length(I) == ncol || throw(ArgumentError("length(I) must equal size(W,2)"))
    nrow == 0 && throw(ArgumentError("weights must be non-empty"))
    ncol == 0 && throw(ArgumentError("weights must be non-empty"))

    # -------- pass 1: validate + total column-weight sum --------
    tot = 0.0
    @inbounds for j in 1:ncol
        Ij = float(I[j])
        (!isfinite(Ij) || Ij < 0.0) &&
            throw(ArgumentError("I[$j] must be finite and ≥ 0 (got $Ij)"))
        colsum = 0.0
        @simd for i in 1:nrow
            wij = float(W[i,j])
            (!isfinite(wij) || wij < 0.0) &&
                throw(ArgumentError("W[$i,$j] must be finite and ≥ 0 (got $wij)"))
            colsum += wij
        end
        tot += Ij * colsum
    end
    tot <= 0.0 && throw(ArgumentError("sum(weights) must be > 0"))

    # -------- pass 2: draw column j --------
    r = rand(rng) * tot
    acc = 0.0
    jsel = ncol
    @inbounds for j in 1:ncol
        Ij = float(I[j])
        colsum = 0.0
        @simd for i in 1:nrow
            colsum += float(W[i,j])
        end
        acc += Ij * colsum
        if r <= acc
            jsel = j
            break
        end
    end

    # -------- conditional draw of row i within column jsel --------
    # Note: I[jsel] cancels here, so use the raw column W[:, jsel].
    r2 = rand(rng)
    acc2 = 0.0
    coltot = 0.0
    @inbounds @simd for i in 1:nrow
        coltot += float(W[i, jsel])
    end
    coltot <= 0.0 && throw(ArgumentError("Selected column has zero weight unexpectedly"))
    r2 *= coltot
    @inbounds @simd for i in 1:nrow
        acc2 += float(W[i, jsel])
        if r2 <= acc2
            return (i, jsel)
        end
    end
    return (nrow, jsel)  # FP fallback
end

@inline wsampleindex_cols(W::AbstractMatrix{<:Real}, I::AbstractVector{<:Real}) =
    wsampleindex_cols(Random.default_rng(), W, I)