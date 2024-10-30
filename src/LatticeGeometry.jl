module LatticeGeometry

import LinearAlgebra
import StaticArrays
include("util.jl")

abstract type AbstractUnitCell end
names(cell::AbstractUnitCell) = throw(AbstractionError(names, cell))
num_atoms(cell::AbstractUnitCell) = length(names(cell))


struct UnitCell{D,T} <: AbstractUnitCell
    atom_names::StaticArrays.SVector{String}
    atom_positions::StaticArrays.SVector{D, T}

end

abstract type AbstractLattice end
struct Lattice <: AbstractLattice

end

end
