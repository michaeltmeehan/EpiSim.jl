abstract type AbstractEpiState end


abstract type AbstractEpiStateSlice end


function Base.convert(::Type{DataFrame}, slices::Vector{<:AbstractEpiStateSlice})
    T = typeof(slices[1])
    names = fieldnames(T)
    cols = NamedTuple{names}(Tuple(getfield.(slices, name) for name in names))
    return DataFrame(cols)
end


function Base.convert(::Type{DataFrame}, slices::Vector{MultiTypeBirthDeathStateSlice})
    n = length(slices[1].I)
    dct = Dict(:t => getfield.(slices, :t))
    for i in 1:n                                                                                                                                                   
        dct[Symbol("I_$i")] = getindex.(getfield.(slices, :I), i)                                                                                                      
        end
    return DataFrame(dct)
end



function Base.show(io::IO, ::MIME"text/plain", slices::Vector{<:AbstractEpiStateSlice})
    df = convert(DataFrame, slices)
    show(io, MIME"text/plain"(), df)
end


function Base.copy(state::AbstractEpiState)
    T = typeof(state)
    T(getfield.(Ref(state), fieldnames(T))...)
end