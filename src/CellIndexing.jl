import StaticArrays as SA
import Base: size, length

abstract type AbstractCellIndexer end
dimension(indexer::AbstractCellIndexer) = throw(AbstractionError(dimension, indexer))
Base.size(indexer::AbstractCellIndexer) = throw(AbstractionError(Base.size, indexer))

Base.length(indexer::AbstractCellIndexer) = prod(size(indexer))

struct DefaultCellIndexer{D} <: AbstractCellIndexer
    cells_per_axis::SA.SVector{D, Int}
    cell_products::SA.SVector{D, Int}

    function DefaultCellIndexer{D}(cells_per_axis::SA.SVector{D, Int}) where D
        prod_vec = SA.@MVector fill(zero(Int), D)
        for idx ∈ range(1, D)
            if idx == 1
                prod_vec[idx] = oneunit(Int)
            else
                prod_vec[idx] = cells_per_axis[idx - 1] * prod_vec[idx - 1]
            end
        end

        return new{D}(cells_per_axis, SA.SVector{D, Int}( prod_vec ))
    end
end
DefaultCellIndexer(dims) = DefaultCellIndexer{length(dims)}(SA.SVector{length(dims), Int}( Int.(dims) ))
const CellIndexer{D} = DefaultCellIndexer{D}

dimension(::DefaultCellIndexer{D}) where D = D
Base.size(indexer::DefaultCellIndexer) = Tuple(indexer.cells_per_axis)
products(indexer::DefaultCellIndexer) = indexer.cell_products

# Keep this abstract because maybe someday I'll add something
abstract type AbstractCellIndices end

struct CellIndices{D} <: AbstractCellIndices
    reduced_coords::SA.SVector{D, Int}
end
CellIndices(coords) = CellIndices{length(coords)}(SA.SVector{length(coords), Int}( Int.(coords) ))
dimension(::CellIndices{D}) where D = D
basis(::CellIndices) = ReducedBasis()
coordinates(indices::CellIndices) = indices.reduced_coords

function index(indexer::DefaultCellIndexer{D1}, indices::CellIndices{D2}) where {D1, D2}
    @assert D1 == D2 "Incorrect dimensions. Got $D1 for the indexer and $D2 for the indices."
    coords = coordinates(indices)
    cell_idx = coords[begin]
    for dim_k ∈ range(1, D1 - 1)
        cell_idx = _forward_index(cell_idx, coords[dim_k + 1], products(indexer)[dim_k + 1])
    end
    return cell_idx
end

function CellIndices(indexer::DefaultCellIndexer{D}, cell_idx) where D
    #= 
        TODO: This is a great spot for a @generated function on the value D 
              ** Update during testing. Apparently with MVectors, julia can 
                 figure out that it can unroll the for loop or something to
                 eliminate allocations. Good until at least D = 20. 
                 
                 ❤️julia
                 (Joe, Jan. 2, 2025)
    =#
    indices = SA.@MVector fill(zero(Int), D)
    leftover_cell = cell_idx
    for dim_k ∈ range(1, D) |> reverse
        (leftover_cell, ind_k) = _backward_index(leftover_cell, products(indexer)[dim_k])
        indices[dim_k] = ind_k
    end
    return CellIndices( SA.SVector(indices) )
end


