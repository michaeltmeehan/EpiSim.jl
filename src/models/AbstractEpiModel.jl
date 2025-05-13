abstract type AbstractEpiModel end

isagentic(model::AbstractEpiModel) = typeof(model.initial_state) <: AgenticState


function Base.show(io::IO, model::T) where {T<:AbstractEpiModel}
    println(io, nameof(T))
    names = fieldnames(T)
    maxlen = maximum(length.(string.(names)))
    for name in names
        val = getfield(model, name)
        padded = rpad(string(name), maxlen)
        println(io, "  $padded: ", val)
    end
end



abstract type AbstractEpiParameters end