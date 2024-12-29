using Test
using LatticeGeometry

@testset "LatticeIndexing.jl" begin
    
    @testset "Single Atom Unitcell" begin
        #=
            Atom_1: 1 2 3 4 5
        =#
        N = 5
        NA = 1

        indexer = LatticeIndexer(; Natoms = NA, Ncells = N)
        
        @test atoms_per_cell(indexer) == NA
        @test num_cells(indexer) == N
        @test num_atoms(indexer) == N * NA

        for cidx ∈ 1:N
            @test index(indexer, LatticeIndices(; atom = 1, cell = cidx)) == cidx
            @test index(indexer, LatticeIndices(; atom = 1, cell = cidx), CellOrdered()) == cidx
            @test index(indexer, LatticeIndices(; atom = 1, cell = cidx), AtomOrdered()) == cidx
        end
    end

    @testset "Two Atoms per Unitcell" begin
        #=
            Atom_1: 1 2 3 4 5
            Atom_2: 1 2 3 4 5
        =#
        NA = 2
        N  = 5

        indexer = LatticeIndexer(; Natoms = NA, Ncells = N)

        @test atoms_per_cell(indexer) == NA
        @test num_cells(indexer) == N
        @test num_atoms(indexer) == N * NA
    end
    
end