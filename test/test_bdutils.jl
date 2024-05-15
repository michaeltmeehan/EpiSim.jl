# Test cases for CRBDParameters
@testset "CRBDParameters Tests" begin
    # Valid parameters
    valid_crbd_params = CRBDParameters(N₀=10, λ=0.5, μ=0.2, ψ=0.1, ρ₀=0.1, r=0.1, t_max=100.0)
    @test check_parameters(valid_crbd_params) == nothing

    # Invalid parameters
    invalid_crbd_params1 = CRBDParameters(N₀=-1, λ=0.5, μ=0.2, ψ=0.1, ρ₀=0.1, r=0.1, t_max=100.0)
    @test_throws AssertionError check_parameters(invalid_crbd_params1)

    invalid_crbd_params2 = CRBDParameters(N₀=10, λ=-0.5, μ=0.2, ψ=0.1, ρ₀=0.1, r=0.1, t_max=100.0)
    @test_throws AssertionError check_parameters(invalid_crbd_params2)

    invalid_crbd_params3 = CRBDParameters(N₀=10, λ=0.5, μ=-0.2, ψ=0.1, ρ₀=0.1, r=0.1, t_max=100.0)
    @test_throws AssertionError check_parameters(invalid_crbd_params3)

    invalid_crbd_params4 = CRBDParameters(N₀=10, λ=0.5, μ=0.2, ψ=-0.1, ρ₀=0.1, r=0.1, t_max=100.0)
    @test_throws AssertionError check_parameters(invalid_crbd_params4)

    invalid_crbd_params5 = CRBDParameters(N₀=10, λ=0.5, μ=0.2, ψ=0.1, ρ₀=1.5, r=0.1, t_max=100.0)
    @test_throws AssertionError check_parameters(invalid_crbd_params5)

    invalid_crbd_params6 = CRBDParameters(N₀=10, λ=0.5, μ=0.2, ψ=0.1, ρ₀=0.1, r=1.5, t_max=100.0)
    @test_throws AssertionError check_parameters(invalid_crbd_params6)

    invalid_crbd_params7 = CRBDParameters(N₀=10, λ=0.5, μ=0.2, ψ=0.1, ρ₀=0.1, r=0.1, t_max=-100.0)
    @test_throws AssertionError check_parameters(invalid_crbd_params7)
end

# Test cases for MTBDParameters
@testset "MTBDParameters Tests" begin
    # Valid parameters
    valid_mtbd_params = MTBDParameters(n_types=2, N₀=[10, 5], λ=[0.5 0.3; 0.2 0.4], μ=[0.2, 0.1], γ=[0.1 0.2; 0.3 0.4], ψ=[0.1, 0.05], ρ₀=[0.1, 0.2], r=[0.1, 0.1], t_max=100.0)
    @test check_parameters(valid_mtbd_params) == nothing

    # Invalid parameters
    invalid_mtbd_params1 = MTBDParameters(n_types=2, N₀=[-10, 5], λ=[0.5 0.3; 0.2 0.4], μ=[0.2, 0.1], γ=[0.1 0.2; 0.3 0.4], ψ=[0.1, 0.05], ρ₀=[0.1, 0.2], r=[0.1, 0.1], t_max=100.0)
    @test_throws AssertionError check_parameters(invalid_mtbd_params1)

    invalid_mtbd_params2 = MTBDParameters(n_types=2, N₀=[10, 5], λ=[-0.5 0.3; 0.2 0.4], μ=[0.2, 0.1], γ=[0.1 0.2; 0.3 0.4], ψ=[0.1, 0.05], ρ₀=[0.1, 0.2], r=[0.1, 0.1], t_max=100.0)
    @test_throws AssertionError check_parameters(invalid_mtbd_params2)

    invalid_mtbd_params3 = MTBDParameters(n_types=2, N₀=[10, 5], λ=[0.5 0.3; 0.2 0.4], μ=[-0.2, 0.1], γ=[0.1 0.2; 0.3 0.4], ψ=[0.1, 0.05], ρ₀=[0.1, 0.2], r=[0.1, 0.1], t_max=100.0)
    @test_throws AssertionError check_parameters(invalid_mtbd_params3)

    invalid_mtbd_params4 = MTBDParameters(n_types=2, N₀=[10, 5], λ=[0.5 0.3; 0.2 0.4], μ=[0.2, 0.1], γ=[0.1 0.2; 0.3 0.4], ψ=[-0.1, 0.05], ρ₀=[0.1, 0.2], r=[0.1, 0.1], t_max=100.0)
    @test_throws AssertionError check_parameters(invalid_mtbd_params4)

    invalid_mtbd_params5 = MTBDParameters(n_types=2, N₀=[10, 5], λ=[0.5 0.3; 0.2 0.4], μ=[0.2, 0.1], γ=[0.1 0.2; 0.3 0.4], ψ=[0.1, 0.05], ρ₀=[1.5, 0.2], r=[0.1, 0.1], t_max=100.0)
    @test_throws AssertionError check_parameters(invalid_mtbd_params5)

    invalid_mtbd_params6 = MTBDParameters(n_types=2, N₀=[10, 5], λ=[0.5 0.3; 0.2 0.4], μ=[0.2, 0.1], γ=[0.1 0.2; 0.3 0.4], ψ=[0.1, 0.05], ρ₀=[0.1, 0.2], r=[1.5, 0.1], t_max=100.0)
    @test_throws AssertionError check_parameters(invalid_mtbd_params6)

    invalid_mtbd_params7 = MTBDParameters(n_types=2, N₀=[10, 5], λ=[0.5 0.3; 0.2 0.4], μ=[0.2, 0.1], γ=[0.1 0.2; 0.3 0.4], ψ=[0.1, 0.05], ρ₀=[0.1, 0.2], r=[0.1, 0.1], t_max=-100.0)
    @test_throws AssertionError check_parameters(invalid_mtbd_params7)
end


@testset "initialize_linelist Tests" begin
    # Test with default type (nothing)
    @testset "Default Type" begin
        N = 10
        linelist = initialize_linelist(N)

        @test size(linelist, 1) == N
        @test all(linelist.event .== 1)
        @test all(linelist.child_id .== collect(1:N))
        @test all(linelist.child_type .== 1)
        @test all(linelist.parent_id .== 0)
        @test all(linelist.parent_type .== 0)
        @test all(linelist.t_birth .== 0.0)
        @test all(linelist.t_death .== Inf)
        @test all(linelist.t_sam .== -1.0)
    end

    # Test with specified type
    @testset "Specified Type" begin
        N = 10
        type = repeat([1, 2], N÷2)
        linelist = initialize_linelist(N, type=type)

        @test size(linelist, 1) == N
        @test all(linelist.event .== 1)
        @test all(linelist.child_id .== collect(1:N))
        @test all(linelist.child_type .== type)
        @test all(linelist.parent_id .== 0)
        @test all(linelist.parent_type .== 0)
        @test all(linelist.t_birth .== 0.0)
        @test all(linelist.t_death .== Inf)
        @test all(linelist.t_sam .== -1.0)
    end

    # Test with non-default size
    @testset "Non-Default Size" begin
        N = 15
        linelist = initialize_linelist(N)

        @test size(linelist, 1) == N
        @test all(linelist.event .== 1)
        @test all(linelist.child_id .== collect(1:N))
        @test all(linelist.child_type .== 1)
        @test all(linelist.parent_id .== 0)
        @test all(linelist.parent_type .== 0)
        @test all(linelist.t_birth .== 0.0)
        @test all(linelist.t_death .== Inf)
        @test all(linelist.t_sam .== -1.0)
    end
end


@testset "initialize_active Tests" begin
    # Test with single type
    @testset "Single Type" begin
        N = [10]
        active, type = initialize_active(N)

        @test length(active) == length(N)
        @test length(type) == sum(N)
        @test all(type .== 1)
        for i in 1:N[1]
            @test active[1][i] == i
        end
    end

    # Test with multiple types
    @testset "Multiple Types" begin
        N = [5, 3]
        active, type = initialize_active(N)

        @test length(active) == length(N)
        @test length(type) == sum(N)
        @test all(type[1:5] .== 1)
        @test all(type[6:8] .== 2)
        for i in 1:N[1]
            @test active[1][i] == i
        end
        for i in 1:N[2]
            @test active[2][i] == N[1] + i
        end
    end

    # Test with non-default sizes
    @testset "Non-Default Sizes" begin
        N = [7, 2, 4]
        active, type = initialize_active(N)

        @test length(active) == length(N)
        @test length(type) == sum(N)
        @test all(type[1:7] .== 1)
        @test all(type[8:9] .== 2)
        @test all(type[10:13] .== 3)
        for i in 1:N[1]
            @test active[1][i] == i
        end
        for i in 1:N[2]
            @test active[2][i] == N[1] + i
        end
        for i in 1:N[3]
            @test active[3][i] == N[1] + N[2] + i
        end
    end
end
