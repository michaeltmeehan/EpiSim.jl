abstract type AbstractEpiModel end


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




"""
    update_event_rates!(event_rates, model, ...)

Update the event rates for the given model.
Must be implemented for each subtype of `AbstractEpiModel`.
"""