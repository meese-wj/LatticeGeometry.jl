

@doc raw"""
    construct_metric(lattice_vectors)

Construct the crystallographic metric tensor ``g_{ij}``

This metric tensor maps distances in the mathematical [`ReducedBasis`](@ref) to 
distances in the physical [`CrystalBasis`](@ref) using the `lattice_vectors`.

Let ``\vec{P}`` and ``\vec{Q}`` be points in space. In the [`CrystalBasis`](@ref),
these vectors are parameterized by 

```math
\begin{aligned}
    \vec{P} &= \sum_{i=1}^d p_i \vec{a}_i,
    \\
    \vec{Q} &= \sum_{i=1}^d q_i \vec{a}_i.
\end{aligned}
```

The collection of coefficients ``\{p_i\}`` and ``\{q_i\}`` are the [`ReducedBasis`](@ref)
coordinates. We denote these vectors with boldface font and write them as 

```math
\begin{aligned}
    \boldsymbol{P} &= \sum_{i=1}^d p_i \hat{e}_i,
    \\
    \boldsymbol{Q} &= \sum_{i=1}^d q_i \hat{e}_i.
\end{aligned}
```

Then the scalar product has 

```math
\vec{P}\cdot\vec{Q} = \sum_{i,j} p_iq_j \vec{a}_i\cdot\vec{a}_j \equiv \boldsymbol{P} \cdot \mathsf{g} \cdot \boldsymbol{Q},
```

where the metric tensor ``\mathsf{g}`` has the components:

```math
g_{ij} = \vec{a}_i \cdot \vec{a}_j.
```
"""
function construct_metric(lattice_vectors)
    dim = length(lattice_vectors)
    @assert all( map( v -> length(v) == dim, lattice_vectors ) ) "We must work with a complete basis."
    gij = zeros(eltype(lattice_vectors[begin]), dim, dim)
    for idx ∈ eachindex(lattice_vectors), jdx ∈ eachindex(lattice_vectors)
        gij[idx, jdx] = LA.dot(lattice_vectors[idx], lattice_vectors[jdx])
    end
    @assert LA.det(gij) ≠ zero(eltype(gij)) "The supplied basis does not span the space and the metric is noninvertible."
    return SA.SMatrix{dim, dim}(gij)
end

abstract type AbstractLatticeStructure end
LA.dot(::AbstractVectorBasis, latt::AbstractLatticeStructure, P, Q) = throw(AbstractionError(LA.norm, latt, P, Q))
LA.dot(::ReducedBasis, P, Q) = LA.dot(P, Q)
LA.dot(::ReducedBasis, latt::AbstractLatticeStructure, P, Q) = dot(ReducedBasis(), P, Q)
LA.dot(::CrystalBasis, g::AbstractMatrix, P, Q) = LA.dot( P, g, Q )
LA.dot(::CrystalBasis, latt::AbstractLatticeStructure, P, Q) = dot(CrystalBasis(), metric(CrystalBasis(), latt), P, Q)

LA.norm(basis::AbstractVectorBasis, args...) = sqrt( dot(basis, args...) )

struct LatticeStructure{D, C <: AbstractUnitCell, T} <: AbstractLatticeStructure
    lattice_constants::SA.SVector{D, T}
    lattice_vectors::SA.SVector{D, SA.SVector{D, T}}
    g_RedToCrys::SA.SMatrix{D, D, T}
    g_CrysToRed::SA.SMatrix{D, D, T}
    unitcell::C
    
    function LatticeStructure(vecs::SA.SVector{D, SA.SVector{D, T}}, unitcell::C) where {D, C, T}
        @assert D == dimension(unitcell) "The $(D)-dimensional lattice vectors must live in the $(dimension(unitcell))-dimensional space as the unit cell of type $(typeof(unitcell))."

        latt_consts = LA.norm.(vecs)
        g_RedToCrys = construct_metric(vecs)
        g_CrysToRed = LA.inv(g_RedToCrys)
        return new{D,C,T}(SA.SVector{D}(latt_consts), vecs, g_RedToCrys, g_CrysToRed, unitcell)
    end
end

metric(::CrystalBasis, latt::LatticeStructure) = latt.g_RedToCrys
metric(::ReducedBasis, latt::LatticeStructure) = latt.g_CrysToRed

