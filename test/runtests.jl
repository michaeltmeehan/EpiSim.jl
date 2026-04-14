using EpiSim
using Test

const ACTIVE_EVENT_KINDS = Set(instances(EventKind))

function assert_eventlog_invariants(el::EventLog)
    n = length(el)
    @test length(el.time) == n
    @test length(el.host) == n
    @test length(el.infector) == n
    @test length(el.kind) == n
    @test all(k -> k in ACTIVE_EVENT_KINDS, el.kind)
end

@testset "EventLog" begin
    el = EventLog(3)

    @test hasproperty(el, :time)
    @test hasproperty(el, :host)
    @test hasproperty(el, :infector)
    @test hasproperty(el, :kind)
    assert_eventlog_invariants(el)
    @test el.time == zeros(3)
    @test el.host == [1, 2, 3]
    @test el.infector == zeros(Int, 3)
    @test el.kind == fill(EK_Seeding, 3)

    EpiSim.update_event_log!(el, 1.5, 4, 2, EK_Transmission)

    assert_eventlog_invariants(el)
    @test length(el) == 4
    @test el.time[end] == 1.5
    @test el.host[end] == 4
    @test el.infector[end] == 2
    @test el.kind[end] == EK_Transmission
end

@testset "sellke smoke" begin
    el = sellke(4, 0, 1, 0.8, 0.25, 0.75, 0.2, 0.5)

    @test el isa EventLog
    assert_eventlog_invariants(el)
    @test !isempty(el.time)
    @test all(isfinite, el.time)
    @test all(el.time .>= 0.0)
end

@testset "gillespie smoke" begin
    el = gillespie(4, 0, 1, 0.8, 1.0, 1.0, 0.5, 0.5)

    @test el isa EventLog
    assert_eventlog_invariants(el)
    @test !isempty(el.time)
    @test all(isfinite, el.time)
    @test all(el.time .>= 0.0)
end

@testset "tree extraction smoke" begin
    el = EventLog(
        [0.0, 1.0],
        [1, 1],
        [0, 0],
        [EK_Seeding, EK_SerialSampling],
    )

    trees = extract_sampled_trees(el, 1)

    @test trees isa Vector{Tree}
    @test length(trees) == 1
    @test length(trees[1]) == 2
    @test trees[1].kind == [NK_Root, NK_SampledLeaf]
end
