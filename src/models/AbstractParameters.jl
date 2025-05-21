abstract type AbstractParameters end


function Base.show(io::IO, p::T) where {T <: AbstractParameters}
    fnames = fieldnames(typeof(p))   
    # Create a named tuple-style string, avoiding module names
    fields_str = join(["$f=$(getfield(p, f))" for f ∈ fnames], ", ")
    print(io, "($fields_str)")
end