# =========================================
# Homogeneous SIR model (no latent stage)
# =========================================

using Distributions

export SIRModel

struct SIRModel <: AbstractEpidemicModel
    dβ::Distribution
    dτi::Distribution
    dτs::Distribution
    r::Float64
end

function SIRModel(β, τi, τs, r::Float64)

    dβ  = β  isa Distribution ? β  : DiscreteNonParametric([β],  [1.0])
    dτi = τi isa Distribution ? τi : Exponential(τi)
    dτs = τs isa Distribution ? τs : Exponential(τs)

    SIRModel(dβ, dτi, dτs, r)
end

infection_slope(model::SIRModel, actives, N) =
    actives.total_β / N

draw_beta(model::SIRModel) = rand(model.dβ)
draw_tau_i(model::SIRModel) = rand(model.dτi)
draw_tau_s(model::SIRModel) = rand(model.dτs)

draw_beta(rng::AbstractRNG, model::SIRModel) = rand(rng, model.dβ)
draw_tau_i(rng::AbstractRNG, model::SIRModel) = rand(rng, model.dτi)
draw_tau_s(rng::AbstractRNG, model::SIRModel) = rand(rng, model.dτs)

has_latent_stage(::SIRModel) = false

# Not used but must exist for dispatch safety
draw_tau_e(::SIRModel) = 0.0
draw_tau_e(rng::AbstractRNG, ::SIRModel) = 0.0


default_stopping(::SIRModel) = NoStopping()