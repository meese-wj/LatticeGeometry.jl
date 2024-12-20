

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
# LA.dot(::AbstractVectorBasis, latt::AbstractLatticeStructure, P, Q) = throw(AbstractionError(LA.norm, latt, P, Q))
# LA.dot(::ReducedBasis, P, Q) = LA.dot(P, Q)
# LA.dot(::ReducedBasis, latt::AbstractLatticeStructure, P, Q) = dot(ReducedBasis(), P, Q)
# LA.dot(::CrystalBasis, g::AbstractMatrix, P, Q) = LA.dot( P, g, Q )
# LA.dot(::CrystalBasis, latt::AbstractLatticeStructure, P, Q) = dot(CrystalBasis(), metric(CrystalBasis(), latt), P, Q)

# LA.norm(basis::AbstractVectorBasis, args...) = sqrt( dot(basis, args...) )

# Necessary Interface Functions
axis_cells(latt::AbstractLatticeStructure) = throw(AbstractionError(axis_cells, latt))
num_cells(latt::AbstractLatticeStructure) = throw(AbstractionError(num_cells, latt))
dimension(latt::AbstractLatticeStructure) = throw(AbstractionError(dimension, latt))

unitcell(latt::AbstractLatticeStructure) = throw(AbstractionError(unitcell, latt))
unitcell_type(latt::AbstractLatticeStructure) = throw(AbstractionError(unitcell_type, latt))

geometer(latt::AbstractLatticeStructure) = throw(AbstractionError(geometer, latt))
geometer_type(latt::AbstractLatticeStructure) = throw(AbstractionError(geometer, latt))

# UnitCell Overloads
num_atoms(latt::AbstractLatticeStructure) = num_cells(latt) * ((num_atoms ∘ unitcell)(latt))

# Geometer Overloads
LA.dot(Pvec, Qvec, latt::AbstractLatticeStructure, basis::AbstractVectorBasis) = LA.dot(Pvec, Qvec, geometer(latt), basis)
LA.norm(Pvec, latt::AbstractLatticeStructure, basis::AbstractVectorBasis) = LA.norm(Pvec, geometer(latt), basis)
cosine_angle_between(Pvec, Qvec, latt::AbstractLatticeStructure, basis::AbstractVectorBasis) = cosine_angle_between(Pvec, Qvec, geometer(latt), basis)
angle_between(Pvec, Qvec, latt::AbstractLatticeStructure, basis::AbstractVectorBasis) = angle_between(Pvec, Qvec, geometer(latt), basis)

basis_vectors(latt::AbstractLatticeStructure) = (basis_vectors ∘ geometer)(latt)
metric(latt::AbstractLatticeStructure) = (metric ∘ geometer)(latt)
cell_volume(latt::AbstractLatticeStructure) = (cell_volume ∘ geometer)(latt)

# Julian gifts for free
function atomic_positions(latt::AbstractLatticeStructure)
    dim = dimension(latt)
    anames = names(unitcell(latt))
    apos_Rcoords = atomic_positions(unitcell(latt))
    pos_type = promote_rule( (atomic_position_type ∘ unitcell)(latt), (basis_vector_type ∘ geometer)(latt) )

    lattice_atoms = Dict{eltype(anames), Vector{SA.SVector{dim, pos_type}}}()
    for (atom, pos) ∈ zip(anames, apos_Rcoords)
        positions = SA.SVector{D, pos_type}[]

        for red_coords ∈ Iterators.product( map( L -> 1:1:L, axis_cells(latt) ) )
            Pvec = pos + SA.SVector{D}(red_coords)
            push!(positions, change_coordinates(Pvec, geometer(latt), ReducedToCrystal()) )
        end

        lattice_atoms[atom] = positions
    end

    return lattice_atoms
end

struct LatticeStructure{D, C <: AbstractUnitCell, G <: AbstractCrystalGeometer} <: AbstractLatticeStructure 
    cells_per_axis::SA.SVector{D, Int}
    unitcell::C
    geometer::G

    function LatticeStructure{D, C, G}(Lsizes, cell::C, geometer::G) where {D, C, G}
        @assert length(Lsizes) == D "The number of cells per axis does not have the right dimension. Got $(Lsizes) which is $(length(Lsizes))-dimensional. Should be $D."
        @assert D == dimension(cell) == dimension(geometer) "The unit cell lives in a different spatial dimension than the lattice? $D ≠ $(dimension(cell)) ≠ $(dimension(geometer))"
        for (idx, L) ∈ enumerate(Lsizes)
            @assert L > zero(L) "The number of cells must be nonnegative. Is it possible for axis $idx to have $L cells?"
        end

        return new{D, C, G}( SA.SVector{D, Int}(Lsizes), cell, geometer)
    end
end

LatticeStructure{C, G}(Lsizes, cell::C, geometer::C) where {C, G} = LatticeStructure{length(Lsizes), C, G}(Lsizes, cell, geometer)

axis_cells(latt::LatticeStructure) = latt.cells_per_axis
dimension(::LatticeStructure{D, C, G}) where {D, C, G} = D
num_cells(latt::LatticeStructure) = prod(latt.cells_per_axis)

unitcell(latt::LatticeStructure) = latt.cell
unitcell_type(::LatticeStructure{D, C, G}) where {D, C, G} = C

geometer(latt::LatticeStructure) = latt.geometer
geometer_type(::LatticeStructure{D, C, G}) where {D, C, G} = G

struct LatticeStructure2{D, C <: AbstractUnitCell, T} <: AbstractLatticeStructure
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

