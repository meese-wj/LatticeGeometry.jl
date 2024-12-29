
abstract type AbstractIndexingScheme end

struct CellOrdered <: AbstractIndexingScheme end
struct AtomOrdered <: AbstractIndexingScheme end

const DefaultIndexing = CellOrdered

struct LatticeIndexer 
    atoms_per_cell::Int
    num_cells::Int
end
LatticeIndexer(; Natoms, Ncells) = LatticeIndexer(Natoms, Ncells)

atoms_per_cell(indexer::LatticeIndexer) = indexer.atoms_per_cell
num_cells(indexer::LatticeIndexer) = indexer.num_cells
num_atoms(indexer::LatticeIndexer) = atoms_per_cell(indexer) * num_cells(indexer)

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
    jdx = oneunit(number) + (number - idx) // N
    return (idx, jdx)
end

function index( indexer::LatticeIndexer, inds::LatticeIndices, ::CellOrdered)
    # Ncells = num_cells(indexer)
    # return cell_idx + (atom_idx - 1) * Ncells
    return _forward_index(cell(inds), atom(inds), num_cells(indexer))
end

function index( indexer::LatticeIndexer, inds::LatticeIndices, ::AtomOrdered)
    # Natoms = atoms_per_cell(indexer)
    # return atom_idx + (cell_idx - 1) * Natoms
    return _forward_index(atom(inds), cell(inds), atoms_per_cell(indexer))
end

index(indexer::LatticeIndexer, inds::LatticeIndices) = index(indexer, inds, DefaultIndexing())

function LatticeIndices(indexer::LatticeIndexer, idx::Number, ::CellOrdered)
    (cell, atom) = _backward_index(idx, num_cells(indexer))
    return LatticeIndices(; atom = atom, cell = cell)
end

function LatticeIndices(indexer::LatticeIndexer, idx::Number, ::AtomOrdered)
    (atom, cell) = _backward_index(idx, atoms_per_cell(indexer))
    return LatticeIndices(; atom = atom, cell = cell)
end

LatticeIndices(indexer::LatticeIndexer, idx::Number) = LatticeIndices(indexer, idx, DefaultIndexing())