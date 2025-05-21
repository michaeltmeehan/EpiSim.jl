abstract type AbstractModel end


isagentic(model::AbstractModel) = typeof(model.initial_state) <: AgenticState

function strip_module_prefixes(s::String)
    replace(s, r"EpiSim\.(EpiEvent|Models)\." => "")
end


struct Model{P<:AbstractParameters, S<:AbstractState} <: AbstractModel
    par::P
    event_types::Vector{DataType}
    initial_state::S
end


function Base.show(io::IO, model::M) where {M <: Model}
    prefix = isagentic(model) ? "Agentic " : ""
    println(io, prefix * string(nameof(M)))

    # Parameters block (indented, multi-line)
    par_str = sprint(show, model.par) |> strip_module_prefixes
    println(io, "  parameters:      $par_str")

    # Event types (compact, short names)
    shortnames = [string(Base.typename(t).name) for t in model.event_types]
    println(io, "  event_types:     [", join(shortnames, ", "), "]")

    # Initial state (compact, single line)
    state_str = sprint(show, model.initial_state) |> strip_module_prefixes
    println(io, "  initial_state:   $state_str")
end


get_default_stop_condition(model::Model{<:AbstractParameters, <:AgenticState}) = s -> s.n_sampled >= 100

get_default_stop_condition(model::Model{<:AbstractParameters, <:AggregateState}) = s -> s.t >= 100.0
