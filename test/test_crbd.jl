# test/test_crbd.jl

@testset "CRBD Model Tests" begin
    params = CRBDParameters(N₀=10, λ=0.5, μ=0.2, ψ=0.1, ρ₀=0.1, r=0.1, t_max=100.0)
    outbreak = simulate_outbreak(params, N_max=1000, S_max=100)
    
    test_linelist_not_empty(outbreak)
    test_linelist_columns(outbreak)
    test_linelist_ids_non_negative(outbreak)
    test_traj_non_negative(outbreak)
    test_parms_match(outbreak, params)
    test_unique_child_ids(outbreak)
    test_valid_times(outbreak)
    test_valid_sampling_times(outbreak)
    test_non_negative_parameters(params)
end
