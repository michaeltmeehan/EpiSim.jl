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
    outbreak = simulate_outbreak(params, N_max=1000, S_max=100)
    
    test_linelist_not_empty(outbreak)
    test_linelist_columns(outbreak)
    test_linelist_ids_non_negative(outbreak)
    test_traj_non_negative(outbreak)
    test_parms_match(outbreak, params)
    test_unique_child_ids(outbreak)
    test_valid_times(outbreak)
    test_valid_sampling_times(outbreak)

   # Additional Tests

    # Test population sizes
    @test all(outbreak.traj[2:end, :] .≥ 0)  # Population sizes should be non-negative at all time points

    # Test event consistency
    @test all(outbreak.linelist.child_type .≤ params.n_types)  # Ensure valid types
    @test all(outbreak.linelist.parent_type .≤ params.n_types)  # Ensure valid types
    
    # Test boundary conditions
    params_zero_pop = MTBDParameters(
        n_types=2,
        N₀=[0, 0],
        λ=[0.5 0.6; 1. 2.],
        μ=[0.2, 0.3],
        γ=[0. 0.1; 0.2 0.],
        ψ=[0.1, 0.1],
        ρ₀=[0.1, 0.2],
        r=[0.1, 0.1],
        t_max=100.0
    )
    outbreak_zero_pop = simulate_outbreak(params_zero_pop, N_max=1000, S_max=100)
    @test isempty(outbreak_zero_pop.linelist)  # No events should occur with zero initial population

    params_max_pop = MTBDParameters(
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
    outbreak_max_pop = simulate_outbreak(params_max_pop, N_max=1, S_max=100)
    @test all(outbreak_max_pop.traj[2:end, :] .<= 1)  # Population size should not exceed N_max=1 at any time point

    params_max_sample = MTBDParameters(
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
    outbreak_max_sample = simulate_outbreak(params_max_sample, N_max=1000, S_max=1)
    @test count(x -> x >= 0, outbreak_max_sample.linelist.t_sam) == 1  # Sampled individuals should not exceed S_max=1

    # Test time progression
    @test all(diff(outbreak.traj[1, :]) .>= 0)  # Ensure time does not decrease
end
