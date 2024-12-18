using Test
using LatticeGeometry
using LinearAlgebra

@testset "CrystalGeometer.jl" begin
    
    @testset "AbstractCrystalGeometer" begin
        
        @testset "Enforced Interface Definitions" begin
            struct _DumbGeometer <: AbstractCrystalGeometer end

            @test_throws LatticeGeometry.AbstractionError dimension(_DumbGeometer())
            @test_throws LatticeGeometry.AbstractionError basis_vectors(_DumbGeometer())
            @test_throws LatticeGeometry.AbstractionError basis_change(_DumbGeometer(), ReducedToCrystal())
            @test_throws LatticeGeometry.AbstractionError basis_change(_DumbGeometer(), CrystalToReduced())
            @test_throws LatticeGeometry.AbstractionError metric(_DumbGeometer(), ReducedBasis())
            @test_throws LatticeGeometry.AbstractionError cell_volume(_DumbGeometer())
        end

        @testset "Method Abstractions" begin
            struct _BoringGeometer <: AbstractCrystalGeometer end

            @test metric(_BoringGeometer(), CrystalBasis()) == LinearAlgebra.I
        end

        @testset "DefaultCrystalGeometer tests" begin
            
            ortho_vecs = [ [1, 0], [0, 2] ]
            ortho_geo = DefaultCrystalGeometer{Float32}(ortho_vecs)
            
            hexa_vecs = [ Float32[1.0, 0], Float32[0.5, 0.5sqrt(3)] ]
            hexa_geo = DefaultCrystalGeometer{Float32}(hexa_vecs)

            @testset "Method Abstractions (2D)" begin
                
                @test dimension(ortho_geo) == 2
                @test basis_vectors(ortho_geo)[1] == ortho_vecs[1]
                @test basis_vectors(ortho_geo)[2] == ortho_vecs[2]
                @test hcat(ortho_vecs...) ≈ basis_change(ortho_geo, ReducedToCrystal())
                @test inv(hcat(ortho_vecs...)) ≈ basis_change(ortho_geo, CrystalToReduced())
                @test cell_volume(ortho_geo) == 2.0
                
                @test dimension(hexa_geo) == 2
                @test basis_vectors(hexa_geo)[1] == hexa_vecs[1]
                @test basis_vectors(hexa_geo)[2] == hexa_vecs[2]
                @test hcat(hexa_vecs...) ≈ basis_change(hexa_geo, ReducedToCrystal())
                @test inv(hcat(hexa_vecs...)) ≈ basis_change(hexa_geo, CrystalToReduced())
                @test cell_volume(hexa_geo) ≈ basis_change(hexa_geo, ReducedToCrystal()) |> det |> abs
            end

            @testset "DefaultCrystalGeometer Constructor tests" begin
                
            end

        end

    end



end