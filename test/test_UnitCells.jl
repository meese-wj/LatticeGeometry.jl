using LatticeGeometry
using Test
using StaticArrays

@testset "UnitCells.jl" begin
    
    @testset "UnitCell" begin
        
        @testset "Constructor Errors" begin
            @test_throws "There must be a name and position" UnitCell( SVector("Fe"), SVector( SVector(0.0, 0.0), SVector(0.0, 1.0) ) )
            
            @test_throws "Each atomic position" UnitCell( SVector("Fe"), (@SVector [SVector(0.0, 2.0)] ) )
            @test_throws "Each atomic position" UnitCell( (@SVector [SVector(0.0, 2.0)] ) )
        end

        @testset "Getter Functionality" begin
            labels = SVector( "Fe", "As" )
            positions = SVector( SVector(0.0, 0.0), SVector(0.5, 0.5) )
            cell = UnitCell(labels, positions)
            
            @test dimension(cell) == length(positions[begin])
            @test (names(cell) .== labels) |> all
        end
    end
end