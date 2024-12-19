import LinearAlgebra as LA
import LinearAlgebra: norm, dot
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

"""
    metric(::AbstractCrystalGeometer, ::CrystalBasis)

In the [`CrystalBasis`](@ref), the `metric` is the identity, 
`LinearAlgebra.I`. 
"""
metric(::AbstractCrystalGeometer, ::CrystalBasis) = LA.I
@doc raw"""
    metric(::AbstractCrystalGeometer, ::ReducedBasis)

In the [`ReducedBasis`](@ref), the `metric` is used to
compute physical _crystal_ angles between vectors measured
in the [`ReducedBasis`](@ref). For example, if ``\boldsymbol{P}``
and ``\boldsymbol{Q}`` denote points within the crystal in the 
[`ReducedBasis`](@ref), which are the 
points ``\vec{P}`` and ``\vec{Q}`` in the [`CrystalBasis`](@ref),
then the scalar product between them is

```math
\vec{P} \cdot \vec{Q} = g_{ij}P_i Q_j,
```

where ``g_{ij}`` is the metric tensor. Quantitatively, it
follows as ``g_{ij} = C_{ki}C_{kj}`` where ``C_{ij} = a_{j,i}`` is the matrix
of basis vectors ``\{\vec{a}_i}\}. Thus, the components of the metric 
tensor are 

```math
g_{ij} = a_{i,k}a_{j,k} = \vec{a}_i \cdot \vec{a}_j.
```
"""
metric(geo::AbstractCrystalGeometer, ::ReducedBasis) = throw(AbstractionError(metric, geo, ReducedBasis()))

"""
    LinearAlgebra.dot(Pvec, Qvec, ::AbstractCrystalGeometer, ::AbstractVectorBasis)

Compute the scalar `dot` product between two vectors `Pvec` and `Qvec`, both measured in 
a particular [`AbstractVectorBasis`](@ref). Use the [`AbstractCrystalGeometer`](@ref) to
determine the [`metric`](@ref) tensor, ``g_{ij}``, such that 

```math
\vec{P} \cdot \vec{Q} = g_{ij}P_i Q_j.
```
"""
LA.dot(Pvec, Qvec, geo::AbstractCrystalGeometer, basis::AbstractVectorBasis) = LA.dot(Pvec, metric(geo, basis), Qvec)

"""
    LinearAlgebra.norm(Pvec, ::AbstractCrystalGeometer, ::AbstractVectorBasis)

Compute the vector `norm` of `Pvec`, measured in the given [`AbstractVectorBasis`](@ref),
using the [`metric`](@ref) tensor. This [`norm`] is simply the square-root of the [`dot`](@ref)
product of `Pvec` with itself.
"""
LA.norm(Pvec, geo::AbstractCrystalGeometer, basis::AbstractVectorBasis) = sqrt( LA.dot(Pvec, Pvec, geo, basis) )

"""
    cosine_angle_between(Pvec, Qvec, ::AbstractCrystalGeometer, ::AbstractVectorBasis)

Compute the `cos`ine of the angle between `Pvec` and `Qvec`, both measured in a given
[`AbstractVectorBasis`](@ref) using a [`metric`](@ref) from the [`AbstractCrystalGeometer`](@ref).
This function relies on the scalar [`dot`](@ref) product and the vector [`norm`](@ref).
"""
cosine_angle_between(Pvec, Qvec, geo::AbstractCrystalGeometer, basis::AbstractVectorBasis) = LA.dot(Pvec, Qvec, geo, basis) / ( LA.norm(Pvec, geo, basis) * LA.norm(Qvec, geo, basis) )

"""
    angle_between(Pvec, Qvec, ::AbstractCrystalGeometer, ::AbstractVectorBasis)

Compute the angle between `Pvec` and `Qvec`, both measured in a given
[`AbstractVectorBasis`](@ref) using a [`metric`](@ref) from the [`AbstractCrystalGeometer`](@ref).
This function relies on the scalar [`dot`](@ref) product and the vector [`norm`](@ref).
"""
angle_between(Pvec, Qvec, geo::AbstractCrystalGeometer, basis::AbstractVectorBasis) = acos( cosine_angle_between(Pvec, Qvec, geo, basis) )

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