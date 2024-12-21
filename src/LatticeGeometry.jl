module LatticeGeometry


include("util.jl")
export showerror

include("UnitCells.jl")
export AbstractUnitCell, UnitCell, 
       names, dimension, atomic_positions, 
       num_atoms, summarize, show

include("CrystalGeometer.jl")
export AbstractVectorBasis, CrystalBasis, ReducedBasis,
       AbstractCrystalGeometer, ReducedToCrystal, CrystalToReduced,
       metric, dot, norm, cosine_angle_between, angle_between,
       basis_vectors, basis_change, cell_volume, change_coordinates,
       DefaultCrystalGeometer

include("LatticeStructures.jl")
export AbstractLatticeStructure,
       LatticeStructure, unitcell_vertices

end
