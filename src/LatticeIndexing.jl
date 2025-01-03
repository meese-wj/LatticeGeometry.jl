
"""
    abstract type AbstractIndexingScheme 

A generic type for indexing the elements of a lattice -- that is,
the [`UnitCell`](@ref)s and the [`atom`](@ref)s within them.
"""
abstract type AbstractIndexingScheme end

@doc raw"""
    CellOrdered <: AbstractIndexingScheme

In this type of [`AbstractIndexingScheme`](@ref) one iterates
over the same [`atom`](@ref) species throughout all [`UnitCell`](@ref)s
before proceeding.

The [`index`](@ref) ordering function is given by 

```math
\mathtt{index}(\mathtt{atom}, \mathtt{cell}) = \mathtt{cell} + (\mathtt{atom} - 1) * N_{\mathrm{cells}}
```

* Note: this is the [`DefaultIndexing`](@ref) scheme.
"""
struct CellOrdered <: AbstractIndexingScheme end

@doc raw"""
    AtomOrdered <: AbstractIndexingScheme

In this type of [`AbstractIndexingScheme`](@ref) one iterates
over the all the [`atom`](@ref) species in a [`UnitCell`](@ref)s
before proceeding.

The [`index`](@ref) ordering function is given by 

```math
\mathtt{index}(\mathtt{atom}, \mathtt{cell}) = \mathtt{atom} + (\mathtt{cell} - 1) * N_{\mathrm{atoms}}
```
"""
struct AtomOrdered <: AbstractIndexingScheme end

const DefaultIndexing = CellOrdered

"""
    LatticeIndexer(; Natoms, Ncells)

This struct is almost a functor which is used 
to [`index`](@ref) [`atom`](@ref) and [`UnitCell`](@ref)s 
in the [`LatticeStructure`](@ref). Its two fields are the
[`atoms_per_cell`](@ref) (`Natoms`) and the number of cells
(`Ncells`) in the [`LatticeStructure`](@ref).
"""
struct LatticeIndexer 
    atoms_per_cell::Int
    num_cells::Int
end
LatticeIndexer(; Natoms, Ncells) = LatticeIndexer(Natoms, Ncells)

atoms_per_cell(indexer::LatticeIndexer) = indexer.atoms_per_cell
num_cells(indexer::LatticeIndexer) = indexer.num_cells
num_atoms(indexer::LatticeIndexer) = atoms_per_cell(indexer) * num_cells(indexer)

"""
    LatticeIndices(; atom, cell)

This is essentially a `NamedTuple` that encloses two `Ints` which
denote the `atom_index` (`atom`) and the `cell_index` (`cell`).
"""
struct LatticeIndices
    atom_index::Int
    cell_index::Int
end

LatticeIndices(; atom, cell) = LatticeIndices(atom, cell)

atom(index::LatticeIndices) = index.atom_index
cell(index::LatticeIndices) = index.cell_index

_forward_index(idx, jdx, N) = idx + (jdx - oneunit(jdx)) * N
function _backward_index(number, N)
    idx = oneunit(number) + (number - oneunit(number)) % N
    jdx = oneunit(number) + (number - idx) ÷ N
    return (idx, jdx)
end

"""
    index(::LatticeIndexer, ::LatticeIndices, ::CellOrdered)

Return the flattened `index` associated with the given [`LatticeIndices`](@ref)
in the [`CellOrdered`](@ref) scheme.

```jldoctest
julia> indexer = LatticeIndexer(; Natoms = 2, Ncells = 5)
LatticeIndexer(2, 5)

julia> index(indexer, LatticeIndices(; atom = 1, cell = 1), CellOrdered())
1

julia> index(indexer, LatticeIndices(; atom = 2, cell = 4), CellOrdered())
9
```
"""
function index( indexer::LatticeIndexer, inds::LatticeIndices, ::CellOrdered)
    # Ncells = num_cells(indexer)
    # return cell_idx + (atom_idx - 1) * Ncells
    return _forward_index(cell(inds), atom(inds), num_cells(indexer))
end

"""
    index(::LatticeIndexer, ::LatticeIndices, ::AtomOrdered)

Return the flattened `index` associated with the given [`LatticeIndices`](@ref)
in the [`AtomOrdered`](@ref) scheme.

```jldoctest
julia> indexer = LatticeIndexer(; Natoms = 2, Ncells = 5)
LatticeIndexer(2, 5)

julia> index(indexer, LatticeIndices(; atom = 1, cell = 1), AtomOrdered())
1

julia> index(indexer, LatticeIndices(; atom = 2, cell = 4), AtomOrdered())
8
```
"""
function index( indexer::LatticeIndexer, inds::LatticeIndices, ::AtomOrdered)
    # Natoms = atoms_per_cell(indexer)
    # return atom_idx + (cell_idx - 1) * Natoms
    return _forward_index(atom(inds), cell(inds), atoms_per_cell(indexer))
end

"""
    index(::LatticeIndexer, ::LatticeIndices [, ::AbstractIndexingScheme = DefaultIndexing()])

Default implementation of the [`index`](@ref)ing scheme for [`LatticeStructure`](@ref)s.

```jldoctest
julia> indexer = LatticeIndexer(; Natoms = 2, Ncells = 5)
LatticeIndexer(2, 5)

julia> index(indexer, LatticeIndices(; atom = 1, cell = 2)) == index(indexer, LatticeIndices(; atom = 1, cell = 2), DefaultIndexing())
true
```
"""
index(indexer::LatticeIndexer, inds::LatticeIndices) = index(indexer, inds, DefaultIndexing())

"""
    LatticeIndices(::LatticeIndexer, idx::Number, ::CellOrdered)

Create a [`LatticeIndices`](@ref) `struct` from the index `idx` in the 
[`CellOrdered`](@ref) scheme.

```jldoctest
julia> indexer = LatticeIndexer(; Natoms = 2, Ncells = 5)
LatticeIndexer(2, 5)

julia> indices = LatticeIndices(indexer, 3, CellOrdered())
LatticeIndices(1, 3)

julia> atom(indices)
1

julia> cell(indices)
3
```
"""
function LatticeIndices(indexer::LatticeIndexer, idx::Number, ::CellOrdered)
    (cell, atom) = _backward_index(idx, num_cells(indexer))
    return LatticeIndices(; atom = atom, cell = cell)
end

"""
    LatticeIndices(::LatticeIndexer, idx::Number, ::AtomOrdered)

Create a [`LatticeIndices`](@ref) `struct` from the index `idx` in the 
[`AtomOrdered`](@ref) scheme.

```jldoctest
julia> indexer = LatticeIndexer(; Natoms = 2, Ncells = 5)
LatticeIndexer(2, 5)

julia> indices = LatticeIndices(indexer, 3, AtomOrdered())
LatticeIndices(1, 2)

julia> atom(indices)
1

julia> cell(indices)
2
```
"""
function LatticeIndices(indexer::LatticeIndexer, idx::Number, ::AtomOrdered)
    (atom, cell) = _backward_index(idx, atoms_per_cell(indexer))
    return LatticeIndices(; atom = atom, cell = cell)
end

"""
    LatticeIndices(::LatticeIndexer, idx::Number)

Default implementation for [`LatticeIndices`](@ref).

```jldoctest
julia> indexer = LatticeIndexer(; Natoms = 2, Ncells = 5)
LatticeIndexer(2, 5)

julia> LatticeIndices(indexer, 3) == LatticeIndices(indexer, 3, DefaultIndexing())
true
```
"""
LatticeIndices(indexer::LatticeIndexer, idx::Number) = LatticeIndices(indexer, idx, DefaultIndexing())