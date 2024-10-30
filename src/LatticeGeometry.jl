module LatticeGeometry


include("util.jl")
export showerror

include("UnitCells.jl")
export AbstractUnitCell, UnitCell, 
       names, dimension, atomic_positions, 
       num_atoms, summarize, show



abstract type AbstractLattice end
struct Lattice <: AbstractLattice

end

end
