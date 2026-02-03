struct State
    time::Float64
    S::Int64
    E::Int64
    I::Int64
    R::Int64
end


function get_prevalence(states::Vector{State})
    return [state.I for state in states]
end


function get_prevalence(states::Vector{State}, times::Vector{Float64})
    preva = Int[]
    state_index = 1
    for t in times
        while state_index < length(states) && states[state_index + 1].time <= t
            state_index += 1
        end
        push!(preva, states[state_index].I)
    end
    return preva
end