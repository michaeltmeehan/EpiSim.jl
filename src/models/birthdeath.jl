# =========================================
# Homogeneous Birth–Death model
# =========================================

using Distributions

export BirthDeathModel

struct BirthDeathModel <: AbstractEpidemicModel
    dβ::Distribution
    dτi::Distribution
    dτs::Distribution
    r::Float64
end


function BirthDeathModel(β, τi, τs, r::Float64)

    dβ  = β isa Distribution ? β : DiscreteNonParametric([β], [1.0])
    dτi = τi isa Distribution ? τi : Exponential(τi)
    dτs = τs isa Distribution ? τs : Exponential(τs)

    BirthDeathModel(dβ, dτi, dτs, r)
end

# None of this is necessary for simulation but is required for dispatch safety and likelihood calculations
# draw_beta(model::BirthDeathModel) = rand(model.dβ)
# draw_tau_i(model::BirthDeathModel) = rand(model.dτi)
# draw_tau_s(model::BirthDeathModel) = rand(model.dτs)

# draw_beta(rng::AbstractRNG, model::BirthDeathModel) = rand(rng, model.dβ)
# draw_tau_i(rng::AbstractRNG, model::BirthDeathModel) = rand(rng, model.dτi)
# draw_tau_s(rng::AbstractRNG, model::BirthDeathModel) = rand(rng, model.dτs)

# # Not used but required for dispatch safety
# draw_tau_e(::BirthDeathModel) = 0.0
# draw_tau_e(rng::AbstractRNG, ::BirthDeathModel) = 0.0


has_latent_stage(::BirthDeathModel) = false


default_stopping(::BirthDeathModel) =
    StopWhenCumulativeInfected(100)