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

        for cidx ∈ 1:N
            @test atom(LatticeIndices(indexer, cidx)) == 1
            @test atom(LatticeIndices(indexer, cidx, AtomOrdered())) == 1
            @test cell(LatticeIndices(indexer, cidx)) == cidx
            @test cell(LatticeIndices(indexer, cidx, AtomOrdered())) == cidx
        end
    end

    @testset "Two Atoms per Unitcell" begin
        #=
            CellOrdered()
                Atom_1: 1 2 3 4 5
                Atom_2: 6 7 8 9 10

            AtomOrdered()
                Cell_1: 1 2
                Cell_2: 3 4
                Cell_3: 5 6
                Cell_4: 7 8
                Cell_5: 9 10
        =#
        NA = 2
        N  = 5

        indexer = LatticeIndexer(; Natoms = NA, Ncells = N)

        @test atoms_per_cell(indexer) == NA
        @test num_cells(indexer) == N
        @test num_atoms(indexer) == N * NA

        @test DefaultIndexing() === CellOrdered() # do this one
        @test index( indexer, LatticeIndices(; atom = 1, cell = 1) ) == 1
        @test index( indexer, LatticeIndices(; atom = 1, cell = 2) ) == 2
        @test index( indexer, LatticeIndices(; atom = 2, cell = 3) ) == 8
        @test index( indexer, LatticeIndices(; atom = 2, cell = 5) ) == 10
        @test LatticeIndices( indexer, 1 ) |> atom == 1
        @test LatticeIndices( indexer, 1 ) |> cell == 1
        @test LatticeIndices( indexer, 7 ) |> atom == 2
        @test LatticeIndices( indexer, 7 ) |> cell == 2
        @test LatticeIndices( indexer, 10 ) |> atom == 2
        @test LatticeIndices( indexer, 10 ) |> cell == 5
        
        @test DefaultIndexing() !== AtomOrdered()
        @test index( indexer, LatticeIndices(; atom = 1, cell = 1), AtomOrdered() ) == 1
        @test index( indexer, LatticeIndices(; atom = 1, cell = 2), AtomOrdered() ) == 3
        @test index( indexer, LatticeIndices(; atom = 2, cell = 3), AtomOrdered() ) == 6
        @test index( indexer, LatticeIndices(; atom = 2, cell = 5), AtomOrdered() ) == 10
        @test LatticeIndices( indexer, 1, AtomOrdered() ) |> atom == 1
        @test LatticeIndices( indexer, 1, AtomOrdered() ) |> cell == 1
        @test LatticeIndices( indexer, 7, AtomOrdered() ) |> atom == 1
        @test LatticeIndices( indexer, 7, AtomOrdered() ) |> cell == 4
        @test LatticeIndices( indexer, 10, AtomOrdered() ) |> atom == 2
        @test LatticeIndices( indexer, 10, AtomOrdered() ) |> cell == 5
        
    end
    
end