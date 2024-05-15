# test/test_mtbd.jl

@testset "MTBD Model Tests" begin
    params = MTBDParameters(
        n_types=2,
        N₀=[1, 0],
        λ=[0.5 0.6; 1. 2.],
        μ=[0.2, 0.3],
        γ=[0. 0.1; 0.2 0.],
        ψ=[0.1, 0.1],
        ρ₀=[0.1, 0.2],
        r=[0.1, 0.1],
        t_max=100.0
    )
    outbreak = simulate(params, N_max=1000, S_max=100)
    
    test_linelist_not_empty(outbreak)
    test_linelist_columns(outbreak)
    test_linelist_ids_non_negative(outbreak)
    test_traj_non_negative(outbreak)
    test_parms_match(outbreak, params)
    test_unique_child_ids(outbreak)
    test_valid_times(outbreak)
    test_valid_sampling_times(outbreak)
end
