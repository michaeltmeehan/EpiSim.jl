function SuperSpreaderModel(; R0::Float64=3.,
                              relative_transmissibility::Float64=0.2,
                              superspreader_fraction::Float64=0.1,
                              recovery_rate::Vector{Float64}=[0.9, 0.9], 
                              sampling_rate::Vector{Float64}=[0.1, 0.1], 
                              I::Vector{Int}=[1, 0], 
                              agentic::Bool=true)
    ρ = relative_transmissibility
    c = superspreader_fraction
    δ = recovery_rate .+ sampling_rate
    λ = [c * R0 * δ[1]         ρ * c * R0 * δ[2];
                  (1. - c) * R0 * δ[1]  ρ * (1. - c) * R0 * δ[2]]
    λ ./= 1. - (1. - ρ) * (1. - c)
    par = MTBDParameters(λ, δ, sampling_rate)
    state = agentic ?
        AgenticMTBDState(; I=I) :
        AggregateMTBDState(; I=I)
    return Model(par, MTBD_EVENT_TYPES, state)
end