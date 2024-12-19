var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = LatticeGeometry","category":"page"},{"location":"#LatticeGeometry","page":"Home","title":"LatticeGeometry","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for LatticeGeometry.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [LatticeGeometry]","category":"page"},{"location":"#LatticeGeometry.AbstractCrystalGeometer","page":"Home","title":"LatticeGeometry.AbstractCrystalGeometer","text":"abstract type AbstractCrystalGeometer end\n\nA generic family of structures/functors that will control the access to metrics and distances based on basis vectors.\n\n\n\n\n\n","category":"type"},{"location":"#LatticeGeometry.AbstractionError","page":"Home","title":"LatticeGeometry.AbstractionError","text":"AbstractionError <: Exception\n\nA pretty generic catch all for abtract typing on the fly.\n\n\n\n\n\n","category":"type"},{"location":"#LatticeGeometry.DefaultCrystalGeometer","page":"Home","title":"LatticeGeometry.DefaultCrystalGeometer","text":"DefaultCrystalGeometer{D, T}\n\nThe go-to AbstractCrystalGeometer one can use for systems of spatial dimension D and with basis vectors of with eltype T <: AbstractFloat.\n\n\n\n\n\n","category":"type"},{"location":"#LatticeGeometry.UnitCell","page":"Home","title":"LatticeGeometry.UnitCell","text":"UnitCell{Natoms, Dimension, Type} <: AbstractUnitCell\n\nGeneric UnitCell structure with Natoms in a given spatial Dimension. Type signifies the numerical type used to store the positions. We use the standard convention of writing the positions in the cell coordinates, where each component of the positions must lie on [0,1].\n\n\n\n\n\n","category":"type"},{"location":"#LatticeGeometry.UnitCell-Union{Tuple{StaticArraysCore.SArray{Tuple{NA}, StaticArraysCore.SVector{D, T}, 1, NA}}, Tuple{T}, Tuple{D}, Tuple{NA}} where {NA, D, T}","page":"Home","title":"LatticeGeometry.UnitCell","text":"UnitCell([atom_names], atomic_positions)\n\nConstructor for the UnitCell. The atomic_positions must be supplied as a StaticArrays.SVector of StaticArrays.SVectors to optimize memory. The atom_names are optional, though if they are supplied then length(atom_names) == length(atomic_positions). If they are not supplied, default atom_names will be used.\n\njulia> using StaticArrays\n\njulia> positions = @SVector [ SVector(0.0, 0.0, 0.0), SVector(0.5, 0.5, 0.5) ];\n\njulia> labels = SVector( \"Fe\", \"As\" );\n\njulia> UnitCell(labels, positions)\n\n  `UnitCell` containing 2 atoms in 3 spatial dimensions.\n\n    Fe at [0.0, 0.0, 0.0]\n    As at [0.5, 0.5, 0.5]\n\njulia> UnitCell(positions; warn_on = false)\n\n  `UnitCell` containing 2 atoms in 3 spatial dimensions.\n\n    Atom_1 at [0.0, 0.0, 0.0]\n    Atom_2 at [0.5, 0.5, 0.5]\n\n\n\n\n\n","category":"method"},{"location":"#LatticeGeometry._positions_within_bounds-Tuple{Any}","page":"Home","title":"LatticeGeometry._positions_within_bounds","text":"_positions_within_bounds(positions; [interval = [0, 1]])\n\nInternal function to check that each position coordinate is on the supplied interval.\n\n\n\n\n\n","category":"method"},{"location":"#LatticeGeometry.angle_between-Tuple{Any, Any, AbstractCrystalGeometer, AbstractVectorBasis}","page":"Home","title":"LatticeGeometry.angle_between","text":"angle_between(Pvec, Qvec, ::AbstractCrystalGeometer, ::AbstractVectorBasis)\n\nCompute the angle between Pvec and Qvec, both measured in a given AbstractVectorBasis using a metric from the AbstractCrystalGeometer. This function relies on the scalar dot product and the vector norm.\n\n\n\n\n\n","category":"method"},{"location":"#LatticeGeometry.basis_change-Tuple{AbstractCrystalGeometer, LatticeGeometry.AbstractBasisTransformation}","page":"Home","title":"LatticeGeometry.basis_change","text":"basis_change(::AbstractCrystalGeometer, ::AbstractBasisTransformation)\n\nReturn the transformation matrix from the AbstractCrystalGeometer for the appropriate AbstractBasisTransformation.\n\n\n\n\n\n","category":"method"},{"location":"#LatticeGeometry.basis_vectors-Tuple{AbstractCrystalGeometer}","page":"Home","title":"LatticeGeometry.basis_vectors","text":"basis_vectors(::AbstractCrystalGeometer)\n\nReturn the lattice basis vectors (in the CrystalBasis).\n\n\n\n\n\n","category":"method"},{"location":"#LatticeGeometry.cell_volume-Tuple{AbstractCrystalGeometer}","page":"Home","title":"LatticeGeometry.cell_volume","text":"cell_volume(::AbstractCrystalGeometer)\n\nReturn the unit-cell volume for a given set of basis_vectors. This is just the determinant of the basis_change(::AbstractCrystalGeometer, ReducedToCrystal()).\n\n\n\n\n\n","category":"method"},{"location":"#LatticeGeometry.construct_metric-Tuple{Any}","page":"Home","title":"LatticeGeometry.construct_metric","text":"construct_metric(lattice_vectors)\n\nConstruct the crystallographic metric tensor g_ij\n\nThis metric tensor maps distances in the mathematical ReducedBasis to  distances in the physical CrystalBasis using the lattice_vectors.\n\nLet vecP and vecQ be points in space. In the CrystalBasis, these vectors are parameterized by \n\nbeginaligned\n    vecP = sum_i=1^d p_i veca_i\n    \n    vecQ = sum_i=1^d q_i veca_i\nendaligned\n\nThe collection of coefficients p_i and q_i are the ReducedBasis coordinates. We denote these vectors with boldface font and write them as \n\nbeginaligned\n    boldsymbolP = sum_i=1^d p_i hate_i\n    \n    boldsymbolQ = sum_i=1^d q_i hate_i\nendaligned\n\nThen the scalar product has \n\nvecPcdotvecQ = sum_ij p_iq_j veca_icdotveca_j equiv boldsymbolP cdot mathsfg cdot boldsymbolQ\n\nwhere the metric tensor mathsfg has the components:\n\ng_ij = veca_i cdot veca_j\n\n\n\n\n\n","category":"method"},{"location":"#LatticeGeometry.cosine_angle_between-Tuple{Any, Any, AbstractCrystalGeometer, AbstractVectorBasis}","page":"Home","title":"LatticeGeometry.cosine_angle_between","text":"cosine_angle_between(Pvec, Qvec, ::AbstractCrystalGeometer, ::AbstractVectorBasis)\n\nCompute the cosine of the angle between Pvec and Qvec, both measured in a given AbstractVectorBasis using a metric from the AbstractCrystalGeometer. This function relies on the scalar dot product and the vector norm.\n\n\n\n\n\n","category":"method"},{"location":"#LatticeGeometry.dimension-Tuple{AbstractCrystalGeometer}","page":"Home","title":"LatticeGeometry.dimension","text":"dimension(::AbstractCrystalGeometer)\n\nReturn the spatial dimension of the AbstractCrystalGeometer.\n\n\n\n\n\n","category":"method"},{"location":"#LatticeGeometry.metric-Tuple{AbstractCrystalGeometer, CrystalBasis}","page":"Home","title":"LatticeGeometry.metric","text":"metric(::AbstractCrystalGeometer, ::CrystalBasis)\n\nIn the CrystalBasis, the metric is the identity,  LinearAlgebra.I. \n\n\n\n\n\n","category":"method"},{"location":"#LatticeGeometry.metric-Tuple{AbstractCrystalGeometer, ReducedBasis}","page":"Home","title":"LatticeGeometry.metric","text":"metric(::AbstractCrystalGeometer, ::ReducedBasis)\n\nIn the ReducedBasis, the metric is used to compute physical crystal angles between vectors measured in the ReducedBasis. For example, if boldsymbolP and boldsymbolQ denote points within the crystal in the  ReducedBasis, which are the  points vecP and vecQ in the CrystalBasis, then the scalar product between them is\n\nvecP cdot vecQ = g_ijP_i Q_j\n\nwhere g_ij is the metric tensor. Quantitatively, it follows as g_ij = C_kiC_kj where C_ij = a_ji is the matrix of basis vectors ``{\\vec{a}_i}}. Thus, the components of the metric  tensor are \n\ng_ij = a_ika_jk = veca_i cdot veca_j\n\n\n\n\n\n","category":"method"},{"location":"#LinearAlgebra.dot-Tuple{Any, Any, AbstractCrystalGeometer, AbstractVectorBasis}","page":"Home","title":"LinearAlgebra.dot","text":"LinearAlgebra.dot(Pvec, Qvec, ::AbstractCrystalGeometer, ::AbstractVectorBasis)\n\nCompute the scalar dot product between two vectors Pvec and Qvec, both measured in  a particular AbstractVectorBasis. Use the AbstractCrystalGeometer to determine the metric tensor, g_ij, such that \n\nvecP cdot vecQ = g_ijP_i Q_j\n\n\n\n\n\n","category":"method"},{"location":"#LinearAlgebra.norm-Tuple{Any, AbstractCrystalGeometer, AbstractVectorBasis}","page":"Home","title":"LinearAlgebra.norm","text":"LinearAlgebra.norm(Pvec, ::AbstractCrystalGeometer, ::AbstractVectorBasis)\n\nCompute the vector norm of Pvec, measured in the given AbstractVectorBasis, using the metric tensor. This [norm] is simply the square-root of the dot product of Pvec with itself.\n\n\n\n\n\n","category":"method"}]
}
