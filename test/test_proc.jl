@testset "n_sampled Tests" begin
    # Test with a mix of sampled and non-sampled individuals
    linelist = DataFrame(
        t_sam = [0.0, -1.0, 1.0, 2.0, -1.0]
    )
    @test n_sampled(linelist) == 2

    # Test with an empty DataFrame
    linelist_empty = DataFrame(
        t_sam = Float64[]
    )
    @test n_sampled(linelist_empty) == 0

    # Test with no sampled individuals
    linelist_no_samples = DataFrame(
        t_sam = [-1.0, -1.0, -1.0]
    )
    @test n_sampled(linelist_no_samples) == 0

    # Test with all individuals sampled
    linelist_all_samples = DataFrame(
        t_sam = [1.0, 2.0, 3.0]
    )
    @test n_sampled(linelist_all_samples) == 3

    linelist_for_outbreak = DataFrame(
        t_sam = [0.0, 1.0, 2.0, -1.0]
    )
    params = CRBDParameters(λ=1., μ=0.5, ψ=0.1, ρ₀=0.5, r=1., t_max = 20.)
    outbreak = Outbreak(params, [1. 2.; 0. 1.], linelist_for_outbreak)

    @test n_sampled(outbreak) == 2
end



### Unit Tests for `n_deceased` Function

@testset "n_deceased Tests" begin
    # Test with a mix of deceased and non-deceased individuals
    linelist = DataFrame(
        t_death = [0.0, -1.0, 1.0, 2.0, Inf]
    )
    @test n_deceased(linelist) == 2

    # Test with an empty DataFrame
    linelist_empty = DataFrame(
        t_death = Float64[]
    )
    @test n_deceased(linelist_empty) == 0

    # Test with no deceased individuals
    linelist_no_deceased = DataFrame(
        t_death = [-1.0, -1.0, -1.0, Inf]
    )
    @test n_deceased(linelist_no_deceased) == 0

    # Test with all individuals deceased
    linelist_all_deceased = DataFrame(
        t_death = [1.0, 2.0, 3.0]
    )
    @test n_deceased(linelist_all_deceased) == 3

    linelist_for_outbreak = DataFrame(
        t_death = [0.0, 1.0, 2.0, -1.0],
        t_sam = [-1., 2., -1., 3.]
    )
    params = CRBDParameters(λ=1., μ=0.5, ψ=0.1, ρ₀=0.5, r=1., t_max = 20.)
    outbreak = Outbreak(params, [1. 2.; 0. 1.], linelist_for_outbreak)

    @test n_deceased(outbreak) == 2
end



### Unit Tests for `offspring_dist` Function
@testset "offspring_dist Tests" begin
    # Test with a mix of deceased and non-deceased individuals, and different parent counts
    linelist = DataFrame(
        parent_id = [0, 1, 1, 2, 3],
        child_id = [1, 2, 3, 4, 5],
        t_death = [Inf, 0.0, 1.0, 2.0, Inf]
    )

    @test offspring_dist(linelist, height=1.0, deceased_only=true) == [1]
    @test offspring_dist(linelist, height=2.0, deceased_only=true) == [1, 1]
    @test offspring_dist(linelist, height=Inf, deceased_only=true) == [1, 1, 0]
    @test offspring_dist(linelist, height=Inf, deceased_only=false) == [2, 1, 1, 0, 0]

    # Test with an empty DataFrame
    linelist_empty = DataFrame(
        parent_id = Int[],
        child_id = Int[],
        t_death = Float64[]
    )
    @test offspring_dist(linelist_empty) == Int[]

    # Test with no offspring
    linelist_no_offspring = DataFrame(
        parent_id = [0, 0, 0, 0, 0],
        child_id = [1, 2, 3, 4, 5],
        t_death = [Inf, Inf, Inf, Inf, Inf]
    )
    @test offspring_dist(linelist_no_offspring) == []

    linelist_for_outbreak = DataFrame(
        parent_id = [0, 1, 1, 2, 3],
        child_id = [1, 2, 3, 4, 5],
        t_death = [Inf, 0.0, 1.0, 2.0, Inf],
        t_sam = [-1., 2., -1., 3., 4.]
    )
    params = CRBDParameters(λ=1., μ=0.5, ψ=0.1, ρ₀=0.5, r=1., t_max = 20.)
    outbreak = Outbreak(params, [1. 2.; 0. 1.], linelist_for_outbreak)

    @test offspring_dist(outbreak.linelist, height=1.0, deceased_only=true) == [1]
end


@testset "type_dist Tests" begin
    # Define a sample MTBDParameters
    mtbd_params = MTBDParameters(
        n_types=3,
        N₀=[5, 0, 0],
        λ=[0.5 0.6 0.0; 0.0 0.5 0.0; 0.0 0.0 0.5],
        μ=[0.2, 0.3, 0.4],
        γ=[0.0 0.1 0.0; 0.0 0.0 0.1; 0.1 0.0 0.0],
        ψ=[0.1, 0.1, 0.1],
        ρ₀=[0.1, 0.1, 0.1],
        r=[0.1, 0.1, 0.1],
        t_max=100.0
    )

    # Test with a mix of types
    linelist = DataFrame(
        child_type = [1, 2, 1, 2, 1],
        t_sam = [0.0, -1.0, 1.0, 2.0, -1.0]
    )
    traj = [0.0 1.0 2.0; 10.0 20.0 30.0]
    outbreak = Outbreak(mtbd_params, traj, linelist)
    @test type_dist(outbreak) == [3, 2, 0]

    # Test with a single type (but multi-type parameters)
    linelist_single_type = DataFrame(
        child_type = [1, 1, 1, 1, 1],
        t_sam = [0.0, -1.0, 1.0, 2.0, -1.0]
    )
    outbreak_single_type = Outbreak(mtbd_params, traj, linelist_single_type)
    @test type_dist(outbreak_single_type) == [5, 0, 0]

    # Test with multiple types but only one present
    linelist_one_present = DataFrame(
        child_type = [1, 1, 1, 1, 1],
        t_sam = [0.0, -1.0, 1.0, 2.0, -1.0]
    )
    outbreak_one_present = Outbreak(mtbd_params, traj, linelist_one_present)
    @test type_dist(outbreak_one_present) == [5, 0, 0]

    # Test with no types (empty linelist)
    linelist_empty = DataFrame(
        child_type = Int[],
        t_sam = Float64[]
    )
    outbreak_empty = Outbreak(mtbd_params, traj, linelist_empty)
    @test type_dist(outbreak_empty) == [0, 0, 0]
end


@testset "summarize Tests" begin
    # Define a sample linelist and current time
    linelist = DataFrame(
        child_id = [1, 2, 3, 4, 5],
        parent_id = [0, 1, 1, 2, 3],
        child_type = [1, 2, 1, 2, 1],
        t_death = [Inf, 0.0, 1.0, 2.0, Inf],
        t_sam = [0.0, -1.0, 1.0, 2.0, -1.0],
        t_birth = [0.0, 0.0, 0.0, 0.0, 0.0]
    )
    t = 3.0

    # Expected outputs
    expected_type = [1, 2, 1, 2, 1]
    expected_deceased = [false, true, true, true, false]
    expected_sampled = [false, false, true, true, false]
    expected_lifespan = [3.0, 0.0, 1.0, 2.0, 3.0]
    expected_offspring = [2, 1, 1, 0, 0]

    type, deceased, sampled, lifespan, offspring = summarize(linelist, t)

    @test type == expected_type
    @test deceased == expected_deceased
    @test sampled == expected_sampled
    @test lifespan == expected_lifespan
    @test offspring == expected_offspring

end
