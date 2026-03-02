# =========================================
# Engine interface
# =========================================

export AbstractSimulationEngine,
       SellkeEngine,
       GillespieEngine,
       simulate

abstract type AbstractSimulationEngine end

struct SellkeEngine <: AbstractSimulationEngine end
struct GillespieEngine <: AbstractSimulationEngine end

"""
Primary simulation entry point.

Dispatches on engine type.
"""
function simulate end