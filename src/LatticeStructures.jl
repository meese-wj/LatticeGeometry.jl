
abstract type AbstractLatticeStructure end

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

change_coordinates(Pvec, latt::AbstractLatticeStructure, basis::AbstractBasisTransformation) = change_coordinates(Pvec, geometer(latt), basis)

# Julian gifts for free
function atomic_positions(latt::AbstractLatticeStructure)
    D = dimension(latt)
    anames = names(unitcell(latt))
    apos_Rcoords = atomic_positions(unitcell(latt))
    pos_type = promote_rule( (atomic_position_type ∘ unitcell)(latt), (basis_vector_type ∘ geometer)(latt) )

    lattice_atoms = Dict{eltype(anames), Vector{SA.SVector{D, pos_type}}}()
    for (atom, pos) ∈ zip(anames, apos_Rcoords)
        positions = SA.SVector{D, pos_type}[]

        # TODO: Maybe may the 1:1:L work for different ranges?
        for cell_idx ∈ range(1, num_cells(latt))
            Pvec = pos + coordinates( CellIndices( cell_indexer(latt), cell_idx ) )
            push!(positions, change_coordinates(Pvec, latt, ReducedToCrystal()) )
        end

        lattice_atoms[atom] = positions
    end

    return lattice_atoms
end

function unitcell_vertices(latt::AbstractLatticeStructure, ::ReducedBasis)
    # TODO: Extend to higher dimensions.
    @assert dimension(latt) == 2 "I haven't extended this to 3D or higher yet."
    vertices = [ [0, 0], [1, 0], [1, 1], [0, 1] ]
    vertices = map(v -> SA.SVector{2}(v), vertices)
    return SA.SVector{length(vertices)}(vertices)
end
unitcell_vertices(latt::AbstractLatticeStructure, ::CrystalBasis) = map( v -> change_coordinates(v, geometer(latt), ReducedToCrystal()), unitcell_vertices(latt, ReducedBasis()) )

struct LatticeStructure{D, C <: AbstractUnitCell, G <: AbstractCrystalGeometer} <: AbstractLatticeStructure 
    cell_indexer::DefaultCellIndexer{D}
    unitcell::C
    geometer::G

    function LatticeStructure{D, C, G}(Lsizes, cell::C, geometer::G) where {D, C, G}
        @assert length(Lsizes) == D "The number of cells per axis does not have the right dimension. Got $(Lsizes) which is $(length(Lsizes))-dimensional. Should be $D."
        @assert D == dimension(cell) == dimension(geometer) "The unit cell lives in a different spatial dimension than the lattice? $D ≠ $(dimension(cell)) ≠ $(dimension(geometer))"
        for (idx, L) ∈ enumerate(Lsizes)
            @assert L > zero(L) "The number of cells must be nonnegative. Is it possible for axis $idx to have $L cells?"
        end

        return new{D, C, G}( DefaultCellIndexer(Lsizes), cell, geometer)
    end
end

LatticeStructure(Lsizes, cell::C, geometer::G) where {C, G} = LatticeStructure{length(Lsizes), C, G}(Lsizes, cell, geometer)

cell_indexer(latt::LatticeStructure) = latt.cell_indexer
axis_cells(latt::LatticeStructure) = SA.SVector(  size( cell_indexer(latt) )  )
dimension(::LatticeStructure{D, C, G}) where {D, C, G} = D
num_cells(latt::LatticeStructure) = length( cell_indexer(latt) )

unitcell(latt::LatticeStructure) = latt.unitcell
unitcell_type(::LatticeStructure{D, C, G}) where {D, C, G} = C

geometer(latt::LatticeStructure) = latt.geometer
geometer_type(::LatticeStructure{D, C, G}) where {D, C, G} = G