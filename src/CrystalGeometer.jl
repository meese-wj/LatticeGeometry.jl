import LinearAlgebra as LA
import LinearAlgebra: norm
import StaticArrays as SA

abstract type AbstractVectorBasis end
struct CrystalBasis <: AbstractVectorBasis end
struct ReducedBasis <: AbstractVectorBasis end

abstract type AbstractBasisTransformation end
struct ReducedToCrystal <: AbstractBasisTransformation end
struct CrystalToReduced <: AbstractBasisTransformation end

"""
    abstract type AbstractCrystalGeometer end

A generic family of structures/functors that will control
the access to metrics and distances based on basis vectors.
"""
abstract type AbstractCrystalGeometer end
metric(::AbstractCrystalGeometer, ::CrystalBasis) = LA.I

"""
    DefaultCrystalGeometer{D, T}

The go-to [`AbstractCrystalGeometer`](@ref) one can use for systems
of spatial dimension `D` and with basis vectors of with `eltype` `T <: AbstractFloat`.
"""
struct DefaultCrystalGeometer{D, T <: AbstractFloat, D2} <: AbstractCrystalGeometer
    basischange_RtoC::SA.SMatrix{D, D, T, D2}  # C = [a1 | a2 | a3 | ... ]
    basischange_CtoR::SA.SMatrix{D, D, T, D2}  # C^-1
    metric_Rbasis::SA.SMatrix{D, D, T, D2}     # C^T * C

    function DefaultCrystalGeometer{D, T, D2}(latt_vectors) where {D, T, D2}
        @assert D2 == D * D "The last template parameter **must** be the first squared. Got $D2 when it should be $(D*D) == ($D)^2."
        @assert length(latt_vectors) == D "The number of basis vectors must match the dimension of the space. Got $(length(latt_vectors)) vectors instead of $D of them."
        Tvecs = SA.SVector{D, T}[]
        for (idx, vec) ∈ enumerate(latt_vectors)
            @assert length(vec) == D "Vector $idx, $vec, is not of the appropriate length. It has $(length(vec)) elements instead of $D."
            push!(Tvecs, SA.SVector{D, T}(vec))
        end

        basischange_RtoC = hcat(Tvecs...)
        basischange_CtoR = LA.inv(basischange_RtoC)
        metric = LA.transpose(basischange_RtoC) * basischange_RtoC

        return new{D, T, D2}( basischange_RtoC, basischange_CtoR, metric )
    end
end

DefaultCrystalGeometer{D, T}(args...) where {D, T} = DefaultCrystalGeometer{D, T, D*D}(args...)
DefaultCrystalGeometer{T}( vectors::Vector ) where T = DefaultCrystalGeometer{length(vectors), T}(vectors)

dimension(::DefaultCrystalGeometer{D, T, D2}) where {D, T, D2} = D
basis_change( geo::DefaultCrystalGeometer, ::ReducedToCrystal ) = geo.basischange_RtoC
basis_change( geo::DefaultCrystalGeometer, ::CrystalToReduced ) = geo.basischange_CtoR
metric(geo::DefaultCrystalGeometer, ::ReducedBasis) = geo.metric_Rbasis
basis_vectors(geo::DefaultCrystalGeometer) = SA.SVector{dimension(geo)}( [ col for col ∈ eachcol( basis_change( geo, ReducedToCrystal() ) ) ] )

cell_volume(geo::DefaultCrystalGeometer) = LA.det( basis_change(geo, ReducedToCrystal()) ) |> abs

