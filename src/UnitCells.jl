import Base: show, names
import LinearAlgebra
import StaticArrays as SA

abstract type AbstractUnitCell end

##############################################
# Enforced Interface Definitions
##############################################
Base.names(cell::AbstractUnitCell) = throw(AbstractionError(names, cell))
dimension(cell::AbstractUnitCell) = throw(AbstractionError(dimension, cell))
atomic_positions(cell::AbstractUnitCell) = throw(AbstractionError(atomic_positions, cell))

##############################################
# Method Abstractions
##############################################
num_atoms(cell::AbstractUnitCell) = length(names(cell))
function summarize(cell::AbstractUnitCell)
    names(cell) isa Base.AbstractVecOrTuple ? nothing : throw(ArgumentError("The names and atomic positions must be `<: AbstractVecOrTuple`. Check the `names` definition.")) 
    atomic_positions(cell) isa Base.AbstractVecOrTuple ? nothing : throw(ArgumentError("The names and atomic positions must be `<: AbstractVecOrTuple`. Check the `atomic_positions` definition.")) 

    output = "\n  `$(typeof(cell).name.wrapper)` containing $(num_atoms(cell)) atoms in $(dimension(cell)) spatial dimensions.\n"
    output *= _tabulate_label_positions(names(cell), atomic_positions(cell))
    return output
end
Base.show(io::IO, cell::AbstractUnitCell) = print(io, summarize(cell))

"""
    UnitCell{Natoms, Dimension, Type} <: AbstractUnitCell

Generic `UnitCell` structure with `Natoms` in a given spatial `Dimension`.
`Type` signifies the numerical type used to store the positions. We use the
standard convention of writing the positions in the `cell coordinates`, where
each component of the positions must lie on [0,1].
"""
struct UnitCell{NA, D, T} <: AbstractUnitCell
    atom_names::SA.SVector{NA, String}
    atom_positions::SA.SVector{NA, SA.SVector{D, T}}

    function UnitCell(anames, apos::SA.SVector{NA, SA.SVector{D, T}}) where {NA, D, T}
        length(anames) == NA ? nothing : throw(ArgumentError("There must be a name and position for each atom in the UnitCell. Got $(anames) of length $(length(anames)) and $(apos) of length $(length(apos))."))
        _positions_within_bounds(apos) ? nothing : throw(ArgumentError("Each atomic position must be listed in the 'cell' coordinates such that each component is in [0,1]. Got $(_tabulate_label_positions(anames, apos))."))
        return new{NA, D, T}(anames, apos)
    end
end

"""
    UnitCell([atom_names], atomic_positions)

Constructor for the `UnitCell`. The `atomic_positions` must be supplied
as a `StaticArrays.SVector` of `StaticArrays.SVector`s to optimize memory.
The `atom_names` are optional, though if they are supplied then `length(atom_names) == length(atomic_positions)`.
If they are not supplied, default `atom_names` will be used.

```jldoctest
julia> using StaticArrays

julia> positions = @SVector [ SVector(0.0, 0.0, 0.0), SVector(0.5, 0.5, 0.5) ];

julia> labels = SVector( "Fe", "As" );

julia> UnitCell(labels, positions)

  `UnitCell` containing 2 atoms in 2 spatial dimensions.

    Fe at [0.0, 0.0, 0.0]
    As at [0.5, 0.5, 0.5]

julia> UnitCell(positions; warn_on = false)

  `UnitCell` containing 2 atoms in 2 spatial dimensions.

    Atom_1 at [0.0, 0.0, 0.0]
    Atom_2 at [0.5, 0.5, 0.5]
```
"""
function UnitCell(apos::SA.SVector{NA, SA.SVector{D, T}}; name_base = "Atom", warn_on = true) where {NA, D, T}
    warn_on ? (@warn "No atomic names were supplied. Creating our own.") : nothing
    anames = [ "$(name_base)_$(atom)" for atom ∈ eachindex(apos) ]
    anames = SA.SVector{NA}(anames)
    return UnitCell(anames, apos)
end

Base.names(cell::UnitCell) = cell.atom_names
dimension(::UnitCell{N, D, T}) where {N, D,T} = D
atomic_positions(cell::UnitCell) = cell.atom_positions