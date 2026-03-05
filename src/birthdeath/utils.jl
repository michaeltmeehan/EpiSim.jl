using StatsPlots
using Statistics
using StatsBase

# --------------------------------------------------
# Derived epidemiological parameters
# --------------------------------------------------

compute_R0(λ, μ, ψ, r) = λ ./ (μ .+ ψ .* r)
compute_delta(λ, μ, ψ, r) = λ .- (μ .+ ψ .* r)

# --------------------------------------------------
# Summary diagnostics
# --------------------------------------------------

bias(est, truth) = mean(est) - truth
rmse(est, truth) = sqrt(mean((est .- truth).^2))

# --------------------------------------------------
# Correlation-based font scaling
# --------------------------------------------------

function corr_fontsize(r)
    # scale between 8 and 22 roughly
    return 8 + 14 * abs(r)
end


function bd_pairs_plot(mles; r=1.0, truth=nothing, bins=30, figsize=(900,900))
    λ  = mles.λ
    μ  = mles.μ
    ψ  = mles.ψ
    R0 = compute_R0(λ, μ, ψ, r)
    δ  = compute_delta(λ, μ, ψ, r)

    data = [
        (:λ,  λ),
        (:μ,  μ),
        (:ψ,  ψ),
        (:R0, R0),
        (:δ,  δ)
    ]

    k = length(data)
    C = [cor(data[i][2], data[j][2]) for i in 1:k, j in 1:k]

    p = plot(layout=(k, k), size=figsize)

    for i in 1:k
        for j in 1:k
            idx = (i - 1) * k + j
            name_row, vec_row = data[i]
            name_col, vec_col = data[j]

            xlim = extrema(vec_col)
            ylim = extrema(vec_row)

            if i == j
                lo, hi = extrema(vec_col)
                if hi ≈ lo
                    lo -= 0.5; hi += 0.5
                end
                edges = range(lo, hi; length=bins + 1)   # keep as Range
edge_width = step(edges)                  # works fine on Range

histogram!(
    p[idx],
    vec_col;
    bins      = collect(edges),           # convert to Vector here
    normalize = :pdf,
    alpha     = 0.4,
    legend    = false,
    xlim      = (lo, hi),
    title     = string(name_col)
)

if truth !== nothing && haskey(truth, name_col)
    tv = truth[name_col]

    vline!(p[idx], [tv]; linewidth=2, color=:red, legend=false)

    b_val = bias(vec_col, tv)
    r_val = rmse(vec_col, tv)

    counts = fit(StatsBase.Histogram, vec_col, collect(edges))   # Vector here too
    peak_y = maximum(counts.weights) / (length(vec_col) * edge_width) * 0.80

    annotate!(
        p[idx],
        mean(vec_col), peak_y,
        text("bias=$(round(b_val,digits=3)) RMSE=$(round(r_val,digits=3))", 7, :left)
    )
end

            elseif i < j
                scatter!(
                    p[idx],         # ← index directly
                    vec_col, vec_row;
                    alpha      = 0.3,
                    markersize = 2,
                    legend     = false,
                    xlim       = xlim,
                    ylim       = ylim
                )
            else
                c_val = C[i, j]
                label = isnan(c_val) ? "NaN" : string(round(c_val, digits=2))
                fsize = isnan(c_val) ? 10 : corr_fontsize(c_val)

                plot!(
                    p[idx];         # ← index directly
                    framestyle = :none,
                    xticks     = false,
                    yticks     = false,
                    legend     = false
                )
                annotate!(
                    p[idx],         # ← index directly
                    0.5, 0.5,
                    text(label, fsize, :center)
                )
            end
        end
    end

    return p
end