function infected_distribution_matrix(
    trajs::Vector{StateTrajectory},
    N::Int
)

    n_logs = length(trajs)
    n_grid = length(trajs[1].t)

    # Matrix: rows = time index, cols = I = 0:N
    M = zeros(Float64, n_grid, N+1)

    for tr in trajs
        @assert length(tr.t) == n_grid "Grid mismatch"

        for j in 1:n_grid
            i = tr.I[j]
            M[j, i+1] += 1.0
        end
    end

    # Convert counts → probabilities
    M ./= n_logs

    return M
end