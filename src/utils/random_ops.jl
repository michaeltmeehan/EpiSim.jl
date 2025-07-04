@inline function pop_random!(rng::AbstractRNG, v::Vector)
    isempty(v) && throw(ArgumentError("Cannot pop from an empty vector."))
    i = rand(rng, 1:length(v))
    v[i], v[end] = v[end], v[i]   # swap with last element
    return pop!(v)
end