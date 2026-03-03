# =========================================
# Abstract epidemic model interface
# =========================================

export AbstractEpidemicModel,
       infection_slope,
       draw_beta,
       draw_tau_e,
       draw_tau_i,
       draw_tau_s,
       has_latent_stage

abstract type AbstractEpidemicModel end

"""
Return instantaneous infection slope (dΛ/dt).
Engine supplies actives and population size.
"""
function infection_slope end

"""
Draw transmission rate β.
"""
function draw_beta end

"""
Draw incubation period.
Only used if has_latent_stage(model) == true.
"""
function draw_tau_e end

"""
Draw infectious duration.
"""
function draw_tau_i end

"""
Draw sampling interval.
"""
function draw_tau_s end

"""
Return true if model has a latent stage (SEIR),
false for SIR-type models.
"""
function has_latent_stage end


default_stopping(::AbstractEpidemicModel) = NoStopping()