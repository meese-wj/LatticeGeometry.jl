
import Base: Exception, showerror
import StaticArrays

"""
    AbstractionError <: Exception

A pretty generic catch all for abtract typing on the fly.
"""
struct AbstractionError <: Base.Exception 
    func::Function
    args::Vector{Type}

    function AbstractionError(func, args...)
        arg_vec = map( a -> typeof(a), [args...] )
        return new(func, arg_vec)
    end
end
function message(err::AbstractionError)
    signature = "($(err.args[begin])"
    for arg ∈ err.args[begin+1:end-1]
        signature *= ", $(arg)"
    end
    if length(err.args) ≥ 2
        signature *= ", $(err.args[end]))"
    end
    return "The function `$(err.func)` has no method definition with the signature: `$(signature)`.\n\n\tAre you sure you implemented the API correctly? 🤨\n"
end

Base.showerror(io::IO, err::AbstractionError) = print(io, "AbstractionError: $(message(err))")

"""
    _positions_within_bounds(positions; [interval = [0, 1]])

Internal function to check that each position coordinate is on the supplied `interval`.
"""
function _positions_within_bounds(positions; interval = StaticArrays.@SVector [0.0, 1.0])
    output = true
    for pos ∈ positions, comp ∈ pos
        output *= interval[begin] ≤ comp ≤ interval[end]
    end
    return output  
end

function _tabulate_label_positions(labels, positions)
    output = ""
    for (label, pos) ∈ zip(labels, positions)
        output *= "\n    $(label) at $(pos)"
    end
    return output
end