using Patches
using Test
using Random
using Logging

### helper_functions.jl
@testset "Testing non_zero_pos" begin
    @test Patches.non_zero_pos([1, 0, 3]) == [(1,), (3,)]
    @test Patches.non_zero_pos([0, 0, 0]) == []
    @test Patches.non_zero_pos([1, 2, 3]) == [(1,), (2,), (3,)]
    @test Patches.non_zero_pos([0 1; 0 3]) == [(1, 2), (2, 2)]
    
    # Edge cases
    @test Patches.non_zero_pos([]) == []
end

@testset "Testing normalize" begin
    @test Patches.normalize([1, 2, 3]) == [1/6, 1/3, 1/2]
    @test Patches.normalize([1, 2, 3], h=1) == [1/6, 1/3, 1/2]
    @test Patches.normalize([1, 2, 3], h=2) == [1/14, 2/14, 3/14]
    @test Patches.normalize([0, 6, 0], h=2) == [0, 6/36, 0]
    
    # Edge cases
    @test_throws DomainError Patches.normalize([])
    @test_throws DomainError Patches.normalize([0])
    @test_throws DomainError Patches.normalize([1, 0, -1])
end

@testset "Testing sample_reaction_indices" begin
    rng = MersenneTwister(1234)
    D = [0 1 0; 0 0 1; 1 0 0]
    @test Patches.sample_reaction_indices(rng, D, 2) == [(1, 2), (2, 3)]
    @test Patches.sample_reaction_indices(rng, D, 3) == [(2, 3), (1, 2), (3, 1)]
    
    # Edge cases
    @test_throws DomainError Patches.sample_reaction_indices(rng, D, 4)
end


### generative_functions.jl
@testset "Testing create_metabolism" begin
    n_complex = 3
    n_simple = 6
    n_monomer = 3
    gapsize = 10
    monomer_content, limited_pathways = create_metabolism(n_complex, n_simple, n_monomer, gapsize)
    @test limited_pathways isa Union{Nothing, Array{Int64, 2}}
    @test length(monomer_content) == n_complex + n_simple + n_monomer
    @test all(monomer_content .>= 1)
    @test all(monomer_content .<= 10 + gapsize)
    
    # Edge cases
    @test_throws DomainError create_metabolism(0, 0, 0, gapsize)

    pool = create_species_pool(100, n_complex, n_simple, n_monomer, monomer_content)
    @test length(pool.consumption_rates) == 100
    @test length(pool.energy) == 100
    @test length(pool.production_matrices) == 100
    @test size(pool.production_matrices[1]) == (n_complex + n_simple + n_monomer, n_complex + n_simple + n_monomer)
    @test all(pool.m .> 0)
    @test all(pool.n_reactions .> 0)
    @test all(pool.n_splits .> 0)

    sample = sample_pool(pool, 10, 5)
    @test sample.n_species == 10
    @test sample.n_invaders == 5
    @test length(sample.production_matrices) == 10 + 5
    @test length(sample.consumption_rates) == 10 + 5
    @test length(sample.energy) == 10 + 5
end

