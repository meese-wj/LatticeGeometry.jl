using Test
using LatticeGeometry
import StaticArrays as SA

@testset "CellIndexing.jl" begin
    
    @testset "AbstractCellIndexer" begin
        struct BoringCellIndexer <: AbstractCellIndexer end

        @test_throws LatticeGeometry.AbstractionError dimension(BoringCellIndexer())
        @test_throws LatticeGeometry.AbstractionError size(BoringCellIndexer())
    end

    @testset "DefaultCellIndexer" begin
        dims = SA.@SVector [4, 8, 16]

        indexer = CellIndexer(dims)

        @test dimension(indexer) == length(dims)
        @test length(indexer) == prod(dims)
        @test size(indexer) == Tuple(dims)
        @test products(indexer) == [ 1, dims[1], dims[1] * dims[2] ]
    end

    @testset "CellIndices" begin
        inds = CellIndices( (1,2,3) )

        @test dimension(inds) == 3
        @test basis(inds) === ReducedBasis()
        @test coordinates(inds) == [1, 2, 3]
    end

    @testset "Indexing Cells" begin
        dims = SA.@SVector [9, 16, 25, 36]
        Ncells = prod(dims)

        indexer = CellIndexer(dims)
        @test length(indexer) == Ncells
        
        bad_site1 = CellIndices( (3, 4, 5) )
        bad_site2 = CellIndices( (3, 4, 5, 6, 7) )
        @test_throws AssertionError index(indexer, bad_site1)
        @test_throws AssertionError index(indexer, bad_site2)
        
        ste = (3, 4, 5, 6)
        cell = 18606 # hard-coded, check it and see🎶
        good_site = CellIndices( ste )
        @test index(indexer, good_site) == ste[1] + (ste[2] - 1) * dims[1] + (ste[3] - 1) * dims[1] * dims[2] + (ste[4] - 1) * dims[1] * dims[2] * dims[3]
        @show @btime CellIndices($indexer, $cell)
    end
end