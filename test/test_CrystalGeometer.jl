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
                Pvec = [ 1, 2, 3 ]
                Qvec = [0.5, 1, 1.5]

                @testset "Cubic Lattice" begin
                    cube_vecs = [ [1, 0, 0], [0, 1, 0], [0, 0, 1] ]
                    cube_geo = DefaultCrystalGeometer{ftype}(cube_vecs)

                    @test LatticeGeometry.basis_vector_type(cube_geo) == ftype
                    @test norm(Pvec, cube_geo, CrystalBasis()) ≈ sqrt( sum( x->x^2, Pvec ))
                    @test norm(Pvec, cube_geo, ReducedBasis()) ≈ sqrt( sum( x->x^2, Pvec ))
                    @test dot(Pvec, Qvec, cube_geo, CrystalBasis()) ≈ sum( Pvec .* Qvec )
                    @test dot(Pvec, Qvec, cube_geo, ReducedBasis()) ≈ sum( Pvec .* Qvec )
                    @test isapprox(cosine_angle_between( Pvec, Qvec, cube_geo, CrystalBasis() ), sum( Pvec .* Qvec ) / (norm(Pvec) * norm(Qvec)); rtol = eps(ftype))
                    @test isapprox(cosine_angle_between( Pvec, Qvec, cube_geo, ReducedBasis() ), sum( Pvec .* Qvec ) / (norm(Pvec) * norm(Qvec)); rtol = eps(ftype))
                    @test isapprox(angle_between( Pvec, Qvec, cube_geo, CrystalBasis() ), acos(sum( Pvec .* Qvec ) / (norm(Pvec) * norm(Qvec))); rtol = eps(ftype))
                    @test isapprox(angle_between( Pvec, Qvec, cube_geo, ReducedBasis() ), acos(sum( Pvec .* Qvec ) / (norm(Pvec) * norm(Qvec))); rtol = eps(ftype))
                    @test change_coordinates( [2.5, 0, 0], cube_geo, ReducedToCrystal() ) ≈ [2.5, 0.0, 0.0]
                    @test change_coordinates( [2.5, 0, 1.1], cube_geo, CrystalToReduced() ) ≈ [2.5, 0.0, 1.1]
                end

                @testset "Orthorhombic Lattice" begin
                    av = [1, 0, 0]
                    bv = [0, 2, 0]
                    cv = [0, 0, 3]
                    ortho_vecs = [ av, bv, cv ]
                    ortho_geo = DefaultCrystalGeometer{Float32}(ortho_vecs)
                    
                    @test LatticeGeometry.basis_vector_type(ortho_geo) == ftype
                    @test norm(Pvec, ortho_geo, CrystalBasis()) ≈ sqrt(sum(x -> x^2, Pvec))
                    @test norm(Pvec, ortho_geo, ReducedBasis()) ≈ sqrt(sum(x -> x^2, Pvec .* av) + sum(x -> x^2, Pvec .* bv) + sum(x -> x^2, Pvec .* cv))
                    @test dot(Pvec, Qvec, ortho_geo, CrystalBasis()) ≈ sum( Pvec .* Qvec )
                    @test dot(Pvec, Qvec, ortho_geo, ReducedBasis()) ≈ sum([ ortho_vecs[idx][idx]^2 * Pvec[idx] * Qvec[idx] for idx ∈ 1:3 ])
                    @test angle_between(Pvec, Qvec, ortho_geo, CrystalBasis()) ≈ acos( sum( Pvec .* Qvec ) / (norm(Pvec) * norm(Qvec)) )
                    @test angle_between(Pvec, Qvec, ortho_geo, ReducedBasis()) ≈ acos(sum([ ortho_vecs[idx][idx]^2 * Pvec[idx] * Qvec[idx] for idx ∈ 1:3 ]) / (norm(Pvec, ortho_geo, ReducedBasis()) * norm(Qvec, ortho_geo, ReducedBasis())))
                    @test change_coordinates([0, 0, 2.3], ortho_geo, ReducedToCrystal()) ≈ [0, 0, 3 * 2.3]
                    @test change_coordinates([1, 1, 1], ortho_geo, CrystalToReduced()) ≈ [1, 1/2, 1/3]
                end

                @testset "Rhombohedral Lattice" begin
                    rha1 = [1, 0, 0]
                    rha2 = 0.5 * [-1, sqrt(3), 0]
                    rhc  = [0, 0, 4]
                    rhomb_vecs = [ rha1, rha2, rhc ]
                    rhomb_geo = DefaultCrystalGeometer{ftype}(rhomb_vecs)

                    Rvec = [1, 1, 0]
                    Svec = Rvec + [0,0,1]
                    Tvec = [1, 0, 0]
                    
                    @test LatticeGeometry.basis_vector_type(rhomb_geo) == ftype
                    @test angle_between( Tvec, [0, 1, 0], rhomb_geo, ReducedBasis() ) ≈ 2π/3  # hard-coded because this result is important
                    @test angle_between( Tvec, [0, 0, 1], rhomb_geo, ReducedBasis() ) ≈ π/2  # hard-coded because this result is important
                    @test angle_between( Tvec, Rvec, rhomb_geo, ReducedBasis() ) ≈ π/3  # hard-coded because this result is important
                    @test angle_between( Rvec, Svec, rhomb_geo, ReducedBasis() ) ≈ atan(rhc[3] / norm(rha2))  # hard-coded because this result is important
                    @test change_coordinates( [1,1,0], rhomb_geo, ReducedToCrystal() ) ≈ (rha1 + rha2)
                    chcoords = change_coordinates( -0.1rha1 + rha2 + 3.67rhc, rhomb_geo, CrystalToReduced() )
                    excoords = [-0.1, 1, 3.67]
                    @test isapprox(chcoords, excoords; rtol = 2eps(ftype))
                end
            end

        end

    end



end