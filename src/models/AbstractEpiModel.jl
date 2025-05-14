abstract type AbstractEpiModel end

isagentic(model::AbstractEpiModel) = typeof(model.initial_state) <: AgenticState


function strip_module_prefixes(s::String)
    replace(s, r"EpiSim\.(EpiEvent|Models)\." => "")
end



function Base.show(io::IO, model::T) where {T <: AbstractEpiModel}
    prefix = isagentic(model) ? "Agentic " : ""
    println(io, prefix * string(nameof(T)))

    # Parameters block (indented, multi-line)
    par_str = sprint(show, model.parameters) |> strip_module_prefixes
    println(io, "  parameters:      $par_str")

    # Event types (compact, short names)
    shortnames = [string(Base.typename(t).name) for t in model.event_types]
    println(io, "  event_types:     [", join(shortnames, ", "), "]")

    # Initial state (compact, single line)
    state_str = sprint(show, model.initial_state) |> strip_module_prefixes
    println(io, "  initial_state:   $state_str")
end


abstract type AbstractEpiParameters end


function Base.show(io::IO, p::T) where {T <: AbstractEpiParameters}
    fnames = fieldnames(typeof(p))   
    # Create a named tuple-style string, avoiding module names
    fields_str = join(["$f=$(getfield(p, f))" for f ∈ fnames], ", ")
    print(io, "($fields_str)")
end
