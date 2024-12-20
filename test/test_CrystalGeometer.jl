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
                ortho_vecs = [ [1, 0], [0, 2] ]
                @test_throws AssertionError DefaultCrystalGeometer{2, Float32, Float64, 3}(ortho_vecs)
                @test_throws AssertionError DefaultCrystalGeometer{3, Float32}(ortho_vecs)
                
                ortho_vecs = [ [1, 0, 0], [0, 0, 1], [1, 1]]
                @test_throws AssertionError DefaultCrystalGeometer{3, Float32}(ortho_vecs)
            end

            @testset "DefaultCrystalGeometer Interface tests" begin
                ftype = Float32
                cube_vecs = [ [1, 0, 0], [0, 1, 0], [0, 0, 1] ]
                cube_geo = DefaultCrystalGeometer{ftype}(cube_vecs)

                Pvec = [ 1, 2, 3 ]
                Qvec = [0.5, 1, 1.5]
                @test norm(Pvec, cube_geo, CrystalBasis()) ≈ sqrt( sum( x->x^2, Pvec ))
                @test norm(Pvec, cube_geo, ReducedBasis()) ≈ sqrt( sum( x->x^2, Pvec ))
                @test dot(Pvec, Qvec, cube_geo, CrystalBasis()) ≈ sum( Pvec .* Qvec )
                @test dot(Pvec, Qvec, cube_geo, ReducedBasis()) ≈ sum( Pvec .* Qvec )
                @test isapprox(cosine_angle_between( Pvec, Qvec, cube_geo, CrystalBasis() ), sum( Pvec .* Qvec ) / (norm(Pvec) * norm(Qvec)); rtol = eps(ftype))
                @test isapprox(cosine_angle_between( Pvec, Qvec, cube_geo, ReducedBasis() ), sum( Pvec .* Qvec ) / (norm(Pvec) * norm(Qvec)); rtol = eps(ftype))
                @test isapprox(angle_between( Pvec, Qvec, cube_geo, CrystalBasis() ), acos(sum( Pvec .* Qvec ) / (norm(Pvec) * norm(Qvec))); rtol = eps(ftype))
                @test isapprox(angle_between( Pvec, Qvec, cube_geo, ReducedBasis() ), acos(sum( Pvec .* Qvec ) / (norm(Pvec) * norm(Qvec))); rtol = eps(ftype))
                
                tetra_vecs = [ [1, 0, 0], [0, 1, 0], [0, 0, 2] ]
                tetra_geo = DefaultCrystalGeometer{Float32}(tetra_vecs)

                # @test norm()
            end

        end

    end



end