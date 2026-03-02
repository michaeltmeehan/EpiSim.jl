# =========================================
# Homogeneous SEIR model
# =========================================

using Distributions

export SEIRModel

struct SEIRModel <: AbstractEpidemicModel
    dβ::Distribution
    dτe::Distribution
    dτi::Distribution
    dτs::Distribution
    r::Float64
end

function SEIRModel(β, τe, τi, τs, r::Float64)

    dβ  = β  isa Distribution ? β  : DiscreteNonParametric([β],  [1.0])
    dτe = τe isa Distribution ? τe : Exponential(τe)
    dτi = τi isa Distribution ? τi : Exponential(τi)
    dτs = τs isa Distribution ? τs : Exponential(τs)

    SEIRModel(dβ, dτe, dτi, dτs, r)
end

infection_slope(model::SEIRModel, actives, N) =
    actives.total_β / N

draw_beta(model::SEIRModel) = rand(model.dβ)
draw_tau_e(model::SEIRModel) = rand(model.dτe)
draw_tau_i(model::SEIRModel) = rand(model.dτi)
draw_tau_s(model::SEIRModel) = rand(model.dτs)

draw_beta(rng::AbstractRNG, model::SEIRModel) = rand(rng, model.dβ)
draw_tau_e(rng::AbstractRNG, model::SEIRModel) = rand(rng, model.dτe)
draw_tau_i(rng::AbstractRNG, model::SEIRModel) = rand(rng, model.dτi)
draw_tau_s(rng::AbstractRNG, model::SEIRModel) = rand(rng, model.dτs)

has_latent_stage(::SEIRModel) = true