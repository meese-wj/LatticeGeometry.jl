var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = LatticeGeometry","category":"page"},{"location":"#LatticeGeometry","page":"Home","title":"LatticeGeometry","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for LatticeGeometry.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [LatticeGeometry]","category":"page"},{"location":"#LatticeGeometry.AbstractCrystalGeometer","page":"Home","title":"LatticeGeometry.AbstractCrystalGeometer","text":"abstract type AbstractCrystalGeometer end\n\nA generic family of structures/functors that will control the access to metrics and distances based on basis vectors.\n\n\n\n\n\n","category":"type"},{"location":"#LatticeGeometry.AbstractIndexingScheme","page":"Home","title":"LatticeGeometry.AbstractIndexingScheme","text":"abstract type AbstractIndexingScheme\n\nA generic type for indexing the elements of a lattice – that is, the UnitCells and the atoms within them.\n\n\n\n\n\n","category":"type"},{"location":"#LatticeGeometry.AbstractionError","page":"Home","title":"LatticeGeometry.AbstractionError","text":"AbstractionError <: Exception\n\nA pretty generic catch all for abtract typing on the fly.\n\n\n\n\n\n","category":"type"},{"location":"#LatticeGeometry.AtomOrdered","page":"Home","title":"LatticeGeometry.AtomOrdered","text":"AtomOrdered <: AbstractIndexingScheme\n\nIn this type of AbstractIndexingScheme one iterates over the all the atom species in a UnitCells before proceeding.\n\nThe index ordering function is given by \n\nmathttindex(mathttatom mathttcell) = mathttatom + (mathttcell - 1) * N_mathrmatoms\n\n\n\n\n\n","category":"type"},{"location":"#LatticeGeometry.CellOrdered","page":"Home","title":"LatticeGeometry.CellOrdered","text":"CellOrdered <: AbstractIndexingScheme\n\nIn this type of AbstractIndexingScheme one iterates over the same atom species throughout all UnitCells before proceeding.\n\nThe index ordering function is given by \n\nmathttindex(mathttatom mathttcell) = mathttcell + (mathttatom - 1) * N_mathrmcells\n\nNote: this is the DefaultIndexing scheme.\n\n\n\n\n\n","category":"type"},{"location":"#LatticeGeometry.DefaultCrystalGeometer","page":"Home","title":"LatticeGeometry.DefaultCrystalGeometer","text":"DefaultCrystalGeometer{D, T}\n\nThe go-to AbstractCrystalGeometer one can use for systems of spatial dimension D and with basis vectors of with eltype T <: AbstractFloat.\n\n** Note     The metric will be stored in a secondary default type of MT = Float64 for      higher precision calculations, particularly of angles between vectors. If, for     some crazy reason, one wants to store the basis vectors as a type T with more     precision than Float64, then MT will be be promoted to T.\n\n\n\n\n\n","category":"type"},{"location":"#LatticeGeometry.LatticeIndexer","page":"Home","title":"LatticeGeometry.LatticeIndexer","text":"LatticeIndexer(; Natoms, Ncells)\n\nThis struct is almost a functor which is used  to index atom and UnitCells  in the LatticeStructure. Its two fields are the atoms_per_cell (Natoms) and the number of cells (Ncells) in the LatticeStructure.\n\n\n\n\n\n","category":"type"},{"location":"#LatticeGeometry.LatticeIndices","page":"Home","title":"LatticeGeometry.LatticeIndices","text":"LatticeIndices(; atom, cell)\n\nThis is essentially a NamedTuple that encloses two Ints which denote the atom_index (atom) and the cell_index (cell).\n\n\n\n\n\n","category":"type"},{"location":"#LatticeGeometry.LatticeIndices-Tuple{LatticeIndexer, Number, AtomOrdered}","page":"Home","title":"LatticeGeometry.LatticeIndices","text":"LatticeIndices(::LatticeIndexer, idx::Number, ::AtomOrdered)\n\nCreate a LatticeIndices struct from the index idx in the  AtomOrdered scheme.\n\n```jldoctest julia> indexer = LatticeIndexer(; Natoms = 2, Ncells = 5) LatticeIndexer(2, 5)\n\njulia> indices = LatticeIndices(indexer, 3, AtomOrdered()) LatticeIndices(1, 2)\n\njulia> atom(indices) 1\n\njulia> cell(indices) 2\n\n\n\n\n\n","category":"method"},{"location":"#LatticeGeometry.LatticeIndices-Tuple{LatticeIndexer, Number, CellOrdered}","page":"Home","title":"LatticeGeometry.LatticeIndices","text":"LatticeIndices(::LatticeIndexer, idx::Number, ::CellOrdered)\n\nCreate a LatticeIndices struct from the index idx in the  CellOrdered scheme.\n\n```jldoctest julia> indexer = LatticeIndexer(; Natoms = 2, Ncells = 5) LatticeIndexer(2, 5)\n\njulia> indices = LatticeIndices(indexer, 3, CellOrdered()) LatticeIndices(1, 3)\n\njulia> atom(indices) 1\n\njulia> cell(indices) 3\n\n\n\n\n\n","category":"method"},{"location":"#LatticeGeometry.LatticeIndices-Tuple{LatticeIndexer, Number}","page":"Home","title":"LatticeGeometry.LatticeIndices","text":"LatticeIndices(::LatticeIndexer, idx::Number)\n\nDefault implementation for LatticeIndices.\n\njulia> indexer = LatticeIndexer(; Natoms = 2, Ncells = 5)\nLatticeIndexer(2, 5)\n\njulia> LatticeIndices(indexer, 3) == LatticeIndices(indexer, 3, DefaultIndexing())\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#LatticeGeometry.UnitCell","page":"Home","title":"LatticeGeometry.UnitCell","text":"UnitCell{Natoms, Dimension, Type} <: AbstractUnitCell\n\nGeneric UnitCell structure with Natoms in a given spatial Dimension. Type signifies the numerical type used to store the positions. We use the standard convention of writing the positions in the cell coordinates, where each component of the positions must lie on [0,1].\n\n\n\n\n\n","category":"type"},{"location":"#LatticeGeometry.UnitCell-Union{Tuple{StaticArraysCore.SArray{Tuple{NA}, StaticArraysCore.SVector{D, T}, 1, NA}}, Tuple{T}, Tuple{D}, Tuple{NA}} where {NA, D, T}","page":"Home","title":"LatticeGeometry.UnitCell","text":"UnitCell([atom_names], atomic_positions)\n\nConstructor for the UnitCell. The atomic_positions must be supplied as a StaticArrays.SVector of StaticArrays.SVectors to optimize memory. The atom_names are optional, though if they are supplied then length(atom_names) == length(atomic_positions). If they are not supplied, default atom_names will be used.\n\njulia> using StaticArrays\n\njulia> positions = @SVector [ SVector(0.0, 0.0, 0.0), SVector(0.5, 0.5, 0.5) ];\n\njulia> labels = SVector( \"Fe\", \"As\" );\n\njulia> UnitCell(labels, positions)\n\n  `UnitCell` containing 2 atoms in 3 spatial dimensions.\n\n    Fe at [0.0, 0.0, 0.0]\n    As at [0.5, 0.5, 0.5]\n\njulia> UnitCell(positions; warn_on = false)\n\n  `UnitCell` containing 2 atoms in 3 spatial dimensions.\n\n    Atom_1 at [0.0, 0.0, 0.0]\n    Atom_2 at [0.5, 0.5, 0.5]\n\n\n\n\n\n","category":"method"},{"location":"#LatticeGeometry._positions_within_bounds-Tuple{Any}","page":"Home","title":"LatticeGeometry._positions_within_bounds","text":"_positions_within_bounds(positions; [interval = [0, 1]])\n\nInternal function to check that each position coordinate is on the supplied interval.\n\n\n\n\n\n","category":"method"},{"location":"#LatticeGeometry.angle_between-Tuple{Any, Any, AbstractCrystalGeometer, AbstractVectorBasis}","page":"Home","title":"LatticeGeometry.angle_between","text":"angle_between(Pvec, Qvec, ::AbstractCrystalGeometer, ::AbstractVectorBasis)\n\nCompute the angle between Pvec and Qvec, both measured in a given AbstractVectorBasis using a metric from the AbstractCrystalGeometer. This function relies on the scalar dot product and the vector norm.\n\n\n\n\n\n","category":"method"},{"location":"#LatticeGeometry.basis_change-Tuple{AbstractCrystalGeometer, LatticeGeometry.AbstractBasisTransformation}","page":"Home","title":"LatticeGeometry.basis_change","text":"basis_change(::AbstractCrystalGeometer, ::AbstractBasisTransformation)\n\nReturn the transformation matrix from the AbstractCrystalGeometer for the appropriate AbstractBasisTransformation.\n\n\n\n\n\n","category":"method"},{"location":"#LatticeGeometry.basis_vectors-Tuple{AbstractCrystalGeometer}","page":"Home","title":"LatticeGeometry.basis_vectors","text":"basis_vectors(::AbstractCrystalGeometer)\n\nReturn the lattice basis vectors (in the CrystalBasis).\n\n\n\n\n\n","category":"method"},{"location":"#LatticeGeometry.cell_volume-Tuple{AbstractCrystalGeometer}","page":"Home","title":"LatticeGeometry.cell_volume","text":"cell_volume(::AbstractCrystalGeometer)\n\nReturn the unit-cell volume for a given set of basis_vectors. This is just the determinant of the basis_change(::AbstractCrystalGeometer, ReducedToCrystal()).\n\n\n\n\n\n","category":"method"},{"location":"#LatticeGeometry.change_coordinates-Tuple{Any, AbstractCrystalGeometer, LatticeGeometry.AbstractBasisTransformation}","page":"Home","title":"LatticeGeometry.change_coordinates","text":"change_coordinates(Pvec, ::AbstractCrystalGeometer, ::AbstractBasisTransformation)\n\nConvert the Pvec vector into a different set of coordinates using the  given AbstractCrystalGeometer and AbstractBasisTransformation. For example, to go from reduced coordinates to the physical crystal coordinates, use\n\nchange_coordinates(Pvec, geo, ReducedToCrystal())  # geo <: AbstractCrystalGeometer\n\nMathematically, the coordinate change in the example above amounts to\n\nvecP = veca_i P_i = C cdot boldsymbolP\n\nwhere boldsymbolP is the vector of reduced coordinates and  C = veca_1  veca_2  dots is the matrix of lattice vectors.\n\n\n\n\n\n","category":"method"},{"location":"#LatticeGeometry.construct_metric-Tuple{Any}","page":"Home","title":"LatticeGeometry.construct_metric","text":"construct_metric(lattice_vectors)\n\nConstruct the crystallographic metric tensor g_ij\n\nThis metric tensor maps distances in the mathematical ReducedBasis to  distances in the physical CrystalBasis using the lattice_vectors.\n\nLet vecP and vecQ be points in space. In the CrystalBasis, these vectors are parameterized by \n\nbeginaligned\n    vecP = sum_i=1^d p_i veca_i\n    \n    vecQ = sum_i=1^d q_i veca_i\nendaligned\n\nThe collection of coefficients p_i and q_i are the ReducedBasis coordinates. We denote these vectors with boldface font and write them as \n\nbeginaligned\n    boldsymbolP = sum_i=1^d p_i hate_i\n    \n    boldsymbolQ = sum_i=1^d q_i hate_i\nendaligned\n\nThen the scalar product has \n\nvecPcdotvecQ = sum_ij p_iq_j veca_icdotveca_j equiv boldsymbolP cdot mathsfg cdot boldsymbolQ\n\nwhere the metric tensor mathsfg has the components:\n\ng_ij = veca_i cdot veca_j\n\n\n\n\n\n","category":"method"},{"location":"#LatticeGeometry.cosine_angle_between-Tuple{Any, Any, AbstractCrystalGeometer, AbstractVectorBasis}","page":"Home","title":"LatticeGeometry.cosine_angle_between","text":"cosine_angle_between(Pvec, Qvec, ::AbstractCrystalGeometer, ::AbstractVectorBasis)\n\nCompute the cosine of the angle between Pvec and Qvec, both measured in a given AbstractVectorBasis using a metric from the AbstractCrystalGeometer. This function relies on the scalar dot product and the vector norm.\n\n\n\n\n\n","category":"method"},{"location":"#LatticeGeometry.dimension-Tuple{AbstractCrystalGeometer}","page":"Home","title":"LatticeGeometry.dimension","text":"dimension(::AbstractCrystalGeometer)\n\nReturn the spatial dimension of the AbstractCrystalGeometer.\n\n\n\n\n\n","category":"method"},{"location":"#LatticeGeometry.index-Tuple{LatticeIndexer, LatticeIndices, AtomOrdered}","page":"Home","title":"LatticeGeometry.index","text":"index(::LatticeIndexer, ::LatticeIndices, ::AtomOrdered)\n\nReturn the flattened index associated with the given LatticeIndices in the AtomOrdered scheme.\n\njulia> indexer = LatticeIndexer(; Natoms = 2, Ncells = 5)\nLatticeIndexer(2, 5)\n\njulia> index(indexer, LatticeIndices(; atom = 1, cell = 1), AtomOrdered())\n1\n\njulia> index(indexer, LatticeIndices(; atom = 2, cell = 4), AtomOrdered())\n8\n\n\n\n\n\n","category":"method"},{"location":"#LatticeGeometry.index-Tuple{LatticeIndexer, LatticeIndices, CellOrdered}","page":"Home","title":"LatticeGeometry.index","text":"index(::LatticeIndexer, ::LatticeIndices, ::CellOrdered)\n\nReturn the flattened index associated with the given LatticeIndices in the CellOrdered scheme.\n\njulia> indexer = LatticeIndexer(; Natoms = 2, Ncells = 5)\nLatticeIndexer(2, 5)\n\njulia> index(indexer, LatticeIndices(; atom = 1, cell = 1), CellOrdered())\n1\n\njulia> index(indexer, LatticeIndices(; atom = 2, cell = 4), CellOrdered())\n9\n\n\n\n\n\n","category":"method"},{"location":"#LatticeGeometry.index-Tuple{LatticeIndexer, LatticeIndices}","page":"Home","title":"LatticeGeometry.index","text":"index(::LatticeIndexer, ::LatticeIndices [, ::AbstractIndexingScheme = DefaultIndexing()])\n\nDefault implementation of the indexing scheme for LatticeStructures.\n\njulia> indexer = LatticeIndexer(; Natoms = 2, Ncells = 5)\nLatticeIndexer(2, 5)\n\njulia> index(indexer, LatticeIndices(; atom = 1, cell = 2)) == index(indexer, LatticeIndices(; atom = 1, cell = 2), DefaultIndexing())\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#LatticeGeometry.metric-Tuple{AbstractCrystalGeometer, CrystalBasis}","page":"Home","title":"LatticeGeometry.metric","text":"metric(::AbstractCrystalGeometer, ::CrystalBasis)\n\nIn the CrystalBasis, the metric is the identity,  LinearAlgebra.I. \n\n\n\n\n\n","category":"method"},{"location":"#LatticeGeometry.metric-Tuple{AbstractCrystalGeometer, ReducedBasis}","page":"Home","title":"LatticeGeometry.metric","text":"metric(::AbstractCrystalGeometer, ::ReducedBasis)\n\nIn the ReducedBasis, the metric is used to compute physical crystal angles between vectors measured in the ReducedBasis. For example, if boldsymbolP and boldsymbolQ denote points within the crystal in the  ReducedBasis, which are the  points vecP and vecQ in the CrystalBasis, then the scalar product between them is\n\nvecP cdot vecQ = g_ijP_i Q_j\n\nwhere g_ij is the metric tensor. Quantitatively, it follows as g_ij = C_kiC_kj where C_ij = a_ji is the matrix of basis vectors ``{\\vec{a}_i}}. Thus, the components of the metric  tensor are \n\ng_ij = a_ika_jk = veca_i cdot veca_j\n\n\n\n\n\n","category":"method"},{"location":"#LinearAlgebra.dot-Tuple{Any, Any, AbstractCrystalGeometer, AbstractVectorBasis}","page":"Home","title":"LinearAlgebra.dot","text":"LinearAlgebra.dot(Pvec, Qvec, ::AbstractCrystalGeometer, ::AbstractVectorBasis)\n\nCompute the scalar dot product between two vectors Pvec and Qvec, both measured in  a particular AbstractVectorBasis. Use the AbstractCrystalGeometer to determine the metric tensor, g_ij, such that \n\nvecP cdot vecQ = g_ijP_i Q_j\n\n\n\n\n\n","category":"method"},{"location":"#LinearAlgebra.norm-Tuple{Any, AbstractCrystalGeometer, AbstractVectorBasis}","page":"Home","title":"LinearAlgebra.norm","text":"LinearAlgebra.norm(Pvec, ::AbstractCrystalGeometer, ::AbstractVectorBasis)\n\nCompute the vector norm of Pvec, measured in the given AbstractVectorBasis, using the metric tensor. This [norm] is simply the square-root of the dot product of Pvec with itself.\n\n\n\n\n\n","category":"method"}]
}
