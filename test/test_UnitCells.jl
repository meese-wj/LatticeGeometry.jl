using LatticeGeometry
using Test
using StaticArrays

@testset "UnitCells.jl" begin

    @testset "AbstractUnitCell" begin
        @testset "Enforced Interface Definitions" begin
            struct _EmptyCell <: AbstractUnitCell end

            @test_throws LatticeGeometry.AbstractionError names(_EmptyCell())
            @test_throws LatticeGeometry.AbstractionError dimension(_EmptyCell())
            @test_throws LatticeGeometry.AbstractionError atomic_positions(_EmptyCell())
        end

        @testset "Method Abstractions" begin
            # Working Example
            struct _BoringCell <: AbstractUnitCell end
            LatticeGeometry.names(::_BoringCell) = @SVector ["I'm boring"]
            LatticeGeometry.dimension(::_BoringCell) = zero(Int)
            LatticeGeometry.atomic_positions(::_BoringCell) = @SVector [zero(Int)]
            
            # Broken Example
            struct _BrokenCell <: AbstractUnitCell end
            LatticeGeometry.names(::_BrokenCell) = @SVector ["I'm broken"]
            LatticeGeometry.dimension(::_BrokenCell) = zero(Int)
            LatticeGeometry.atomic_positions(::_BrokenCell) = zero(Int) # This breaks mutual iterabilility with names
 
            @testset "Check EID Implementations" begin
                @test length(names(_BoringCell())) == 1
                @test dimension(_BoringCell()) == 0
                @test (atomic_positions(_BoringCell()) .== [0]) |> all
                
                @test length(names(_BrokenCell())) == 1
                @test dimension(_BrokenCell()) == 0
                @test atomic_positions(_BrokenCell()) == 0 
            end
            @testset "Check Method Abstractions" begin
                @test num_atoms(_BoringCell()) == 1
                @test summarize(_BoringCell()) == "\n  `_BoringCell` containing 1 atoms in 0 spatial dimensions.\n\n    I'm boring at 0"

                @test num_atoms(_BrokenCell()) == 1
                @test_throws "AbstractVecOrTuple" summarize(_BrokenCell())
            end
        end
    end
    
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
            @test (atomic_positions(cell) .== positions) |> all
        end
    end
end