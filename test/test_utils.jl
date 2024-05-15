# test/test_utils.jl
"""
Test that the linelist is not empty.
"""
function test_linelist_not_empty(outbreak::Outbreak)
    @test !isempty(outbreak.linelist)
end

"""
Test that the linelist contains all required columns.
"""
function test_linelist_columns(outbreak::Outbreak)
    required_columns = [:event, :child_id, :child_type, :parent_id, :parent_type, :t_birth, :t_death, :t_sam]
    for col in required_columns
        @test hasproperty(outbreak.linelist, col)
    end
end

"""
Test that all IDs in the linelist are non-negative.
"""
function test_linelist_ids_non_negative(outbreak::Outbreak)
    @test all(outbreak.linelist.child_id .>= 0)
    @test all(outbreak.linelist.parent_id .>= 0)
end

"""
Test that the traj array contains only non-negative numbers.
"""
function test_traj_non_negative(outbreak::Outbreak)
    @test all(outbreak.traj .>= 0)
end

"""
Test that the parms field matches the input parameters.
"""
function test_parms_match(outbreak::Outbreak, input_params::EpiParameters)
    @test outbreak.parms == input_params
end

"""
Test that each child_id in the linelist is unique.
"""
function test_unique_child_ids(outbreak::Outbreak)
    @test length(outbreak.linelist.child_id) == length(unique(outbreak.linelist.child_id))
end

"""
Test that t_death is greater than or equal to t_birth.
"""
function test_valid_times(outbreak::Outbreak)
    @test all(outbreak.linelist.t_death .>= outbreak.linelist.t_birth)
end

"""
Test that t_sam is either negative or greater than or equal to t_birth.
"""
function test_valid_sampling_times(outbreak::Outbreak)
    @test all((outbreak.linelist.t_sam .< 0) .| (outbreak.linelist.t_sam .>= outbreak.linelist.t_birth))
end

"""
Test that all parameters in parms are non-negative.
"""
function test_non_negative_parameters(params::EpiParameters)
    for field in fieldnames(typeof(params))
        @test getfield(params, field) >= 0
    end
end
