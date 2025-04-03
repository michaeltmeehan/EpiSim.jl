abstract type AbstractEpiModel end


"""
    update_event_rates!(event_rates, model, ...)

Update the event rates for the given model.
Must be implemented for each subtype of `AbstractEpiModel`.
"""


"""
    simulate_chain(model::AbstractEpiModel; kwargs...) -> TransmissionChain

Simulate a stochastic epidemic outbreak under the specified epidemiological model.

This is an abstract interface function. Concrete subtypes of `AbstractEpiModel` must implement
a specific method for `simulate_chain`, defining how the outbreak proceeds under the dynamics
of that model (e.g., SIR, SEIR, birth-death, etc.).

The function typically returns a `TransmissionChain` object, which records the chain of infections,
their timing, and the identities of sampled individuals.

# Arguments
- `model`: an instance of a subtype of `AbstractEpiModel`, defining the parameters and structure
           of the epidemic model.
- `kwargs...`: optional keyword arguments passed to the specific simulation method, such as:
    - `S_init`: initial number of susceptibles
    - `I_init`: initial number of infectives
    - `S_max`: maximum number of samples to collect
    - `N`: total population size
    - etc.

# Returns
- A `TransmissionChain` recording who infected whom, when, and who was sampled.

# Notes
- This function should be specialized for each model type.
- The base method will throw a `MethodError` if not implemented for the given `model` type.

# Example
```julia
model = SIRModel(β, γ, ψ, N)
chain = simulate_chain(model; S_init=999, I_init=1, S_max=100)
```
"""

