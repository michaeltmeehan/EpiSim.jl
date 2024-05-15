@testset "Outbreak Tests" begin
    # Example parameters
    params = CRBDParameters(N₀=10, λ=0.5, μ=0.2, ψ=0.1, ρ₀=0.1, r=0.1, t_max=100.0)
    traj = [0.0 10.0; 1.0 9.0; 2.0 8.0]
    linelist = DataFrame(
        event = [1, 1, 1],
        child_id = [1, 2, 3],
        child_type = [1, 1, 1],
        parent_id = [0, 1, 1],
        parent_type = [1, 1, 1],
        t_birth = [0.0, 0.5, 1.0],
        t_death = [1.0, 1.5, 2.0],
        t_sam = [-1.0, -1.0, -1.0]
    )
    outbreak = Outbreak(params, traj, linelist)
    
    # Test creation of Outbreak instance
    @test typeof(outbreak) == Outbreak
    @test typeof(outbreak.parms) <: EpiParameters
    @test typeof(outbreak.traj) == Matrix{Float64}
    @test typeof(outbreak.linelist) == DataFrame

    # Test n_sampled function
    linelist.t_sam[3] = 1.5
    outbreak = Outbreak(params, traj, linelist)
    @test n_sampled(outbreak) == 1  # One individual sampled

    # Test the show method
    io = IOBuffer()
    show(io, outbreak)
    output = String(take!(io))
    @test occursin("Outbreak Summary", output)
    @test occursin("Initial pop.", output)
    @test occursin("Final pop. (cum)", output)
    @test occursin("Sampled", output)
    @test occursin("Time span", output)
end
