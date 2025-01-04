using LatticeGeometry
using Test
using StaticArrays

@testset "UnitCells.jl" begin

    @testset "Helper Functions" begin
        @testset "_convert_vectors_into_SVectors" begin
            vec_vectors = [ [1.0, 2.0, 3.0], [4.0, 5.0, 6.0] ]
            tup_vectors = ( [1.0, 2.0, 3.0], [4.0, 5.0, 6.0] )
            svec_vectors = LatticeGeometry._convert_vectors_into_SVectors(vec_vectors)
            svec_tup = LatticeGeometry._convert_vectors_into_SVectors(tup_vectors)
            
            @test typeof(svec_vectors) === SVector{2, SVector{3, Float64}}
            @test typeof(svec_tup) === SVector{2, SVector{3, Float64}}
            @test svec_vectors[begin] == vec_vectors[begin]
            @test svec_vectors[end] == vec_vectors[end]
            @test svec_tup[begin] == vec_vectors[begin]
            @test svec_tup[end] == vec_vectors[end]
        end

    end

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
            @test LatticeGeometry.atomic_position_type(cell) === Float64
        end

        @testset "Display Functionality" begin
            warren_buffer = IOBuffer()
            cell = UnitCell( ([1.0]) )

            show(warren_buffer, cell)
            warren_str = String(take!(warren_buffer))

            @test warren_str == summarize(cell)
        end
    end
end