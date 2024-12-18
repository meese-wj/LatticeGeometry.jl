module LatticeGeometry


include("util.jl")
export showerror

include("UnitCells.jl")
export AbstractUnitCell, UnitCell, 
       names, dimension, atomic_positions, 
       num_atoms, summarize, show

include("CrystalGeometer.jl")

include("LatticeStructures.jl")
export AbstractVectorBasis, CrystalBasis, ReducedBasis,
       construct_metric, 
       AbstractLatticeStructure, dot, norm
       LatticeStructure, metric

end
