
import Base: Exception, showerror
export showerror

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