abstract type AbstractState end

abstract type AgenticState <: AbstractState end
abstract type AggregateState <: AbstractState end

# TODO: It may be sensible to define AbstractEpiState as a subtype of AbstractState (but then what would I do with AgenticState and AggregateState?)
capture(state::AbstractState) = (; t=state.t, I=state.I)


function initialize_event_log(state::AgenticState)::Vector{AbstractEpiEvent}
    event_log = Vector{AbstractEpiEvent}()
    for i in 1:state.n_cumulative
        push!(event_log, Seed(i, 0.0))
    end
    return event_log
end


function initialize_event_log(state::AggregateState)::Vector{AbstractEpiEvent}
    event_log = Vector{AbstractEpiEvent}()
    for i in 1:sum(state.I)
        push!(event_log, Seed(nothing, 0.0))
    end
    return event_log
end