using Test
using EpiSim

@test isdefined(EpiSim, :sellke)
@test isdefined(EpiSim, :gillespie)
@test isdefined(EpiSim, :EventLog)
@test isdefined(EpiSim, :EventKind)