using EpiSim
using Test
using Random
using Distributions: DiscreteNonParametric

const ACTIVE_EVENT_KINDS = Set(instances(EventKind))

function assert_eventlog_invariants(el::EventLog)
    n = length(el)
    @test length(el.time) == n
    @test length(el.host) == n
    @test length(el.infector) == n
    @test length(el.kind) == n
    @test all(k -> k in ACTIVE_EVENT_KINDS, el.kind)
    @test validate_event_log(el; throw=false)
end

function mean_summary(simulator, nrep::Int; seed::Int)
    Random.seed!(seed)
    final_sizes = Vector{Int}(undef, nrep)
    transmissions = Vector{Int}(undef, nrep)

    for i in 1:nrep
        el = simulator()
        validate_event_log(el)
        final_sizes[i] = final_size(el)
        transmissions[i] = event_count(el, EK_Transmission)
    end

    return (; mean_final_size=sum(final_sizes) / nrep, mean_transmissions=sum(transmissions) / nrep)
end

function seir_state_at(el::EventLog, t::Float64, S0::Int, E0::Int, I0::Int)
    S = S0
    E = E0
    I = I0

    for i in eachindex(el.time)
        el.time[i] > t && break
        kind = el.kind[i]
        if kind == EK_Transmission
            S -= 1
            E += 1
        elseif kind == EK_Activation
            E -= 1
            I += 1
        elseif kind == EK_Removal || kind == EK_SerialSampling
            I -= 1
        end
    end

    return (; S, E, I)
end

function exact_seir_distribution_uniformization(population_size::Int, S0::Int, E0::Int, I0::Int,
                                                β::Float64, α::Float64, γ::Float64, t::Float64;
                                                atol::Float64=1e-12)
    states = Tuple{Int,Int,Int}[]
    index = Dict{Tuple{Int,Int,Int},Int}()

    for S in 0:population_size, E in 0:(population_size - S), I in 0:(population_size - S - E)
        state = (S, E, I)
        index[state] = length(states) + 1
        push!(states, state)
    end

    rates = [Tuple{Int,Float64}[] for _ in states]
    max_exit_rate = 0.0
    for (j, (S, E, I)) in pairs(states)
        exit_rate = 0.0
        if S > 0 && I > 0
            rate = β * S * I / population_size
            push!(rates[j], (index[(S - 1, E + 1, I)], rate))
            exit_rate += rate
        end
        if E > 0
            rate = α * E
            push!(rates[j], (index[(S, E - 1, I + 1)], rate))
            exit_rate += rate
        end
        if I > 0
            rate = γ * I
            push!(rates[j], (index[(S, E, I - 1)], rate))
            exit_rate += rate
        end
        max_exit_rate = max(max_exit_rate, exit_rate)
    end

    p = zeros(Float64, length(states))
    p[index[(S0, E0, I0)]] = 1.0
    max_exit_rate == 0.0 && return states, p

    result = zeros(Float64, length(states))
    poisson_weight = exp(-max_exit_rate * t)
    poisson_cdf = poisson_weight
    result .+= poisson_weight .* p
    k = 0

    while 1.0 - poisson_cdf > atol
        next_p = similar(p)
        fill!(next_p, 0.0)

        for (j, state_rates) in pairs(rates)
            stay_probability = 1.0
            for (dest, rate) in state_rates
                transition_probability = rate / max_exit_rate
                next_p[dest] += p[j] * transition_probability
                stay_probability -= transition_probability
            end
            next_p[j] += p[j] * stay_probability
        end

        k += 1
        poisson_weight *= max_exit_rate * t / k
        poisson_cdf += poisson_weight
        result .+= poisson_weight .* next_p
        p = next_p
    end

    result ./= sum(result)
    return states, result
end

function infectious_count_distribution(states, probabilities, population_size::Int)
    distribution = zeros(Float64, population_size + 1)
    for i in eachindex(states)
        infectious_count = states[i][3]
        distribution[infectious_count + 1] += probabilities[i]
    end
    return distribution
end

function empirical_infectious_count_distribution(simulator, nrep::Int, t::Float64,
                                                 population_size::Int, S0::Int, E0::Int, I0::Int;
                                                 seed::Int)
    Random.seed!(seed)
    counts = zeros(Int, population_size + 1)

    for _ in 1:nrep
        el = simulator()
        validate_event_log(el; population_size=population_size)
        state = seir_state_at(el, t, S0, E0, I0)
        counts[state.I + 1] += 1
    end

    return counts ./ nrep
end

total_variation_distance(p, q) = 0.5 * sum(abs.(p .- q))

@testset "public export surface" begin
    exported = Set(names(EpiSim))

    for name in (
        :EventLog,
        :gillespie,
        :sellke,
        :event_times,
        :event_kind,
        :host_event_summary,
        :event_time_state_counts,
        :transmission_tree,
        :transmission_edges,
        :transmission_chain,
        :run_ensemble,
        :ensemble_aggregate_summary,
    )
        @test name in exported
    end

    @test :EK_None ∉ exported
    @test :popr! ∉ exported
    @test :wsample ∉ exported
    @test :wsampleindex ∉ exported
    @test :wsampleindex_cols ∉ exported

    @test EpiSim.EK_None isa EventKind
    @test EpiSim.wsample([:a, :b], [1.0, 0.0]) == :a
end

@testset "polish wording consistency" begin
    root = dirname(@__DIR__)
    read_package_file(parts...) = read(joinpath(root, parts...), String)

    readme = read_package_file("README.md")
    events_src = read_package_file("src", "core", "events.jl")
    transmission_src = read_package_file("src", "analysis", "transmission_views.jl")
    ensemble_src = read_package_file("src", "analysis", "ensemble.jl")

    @test occursin("Simulation functions return an `EventLog`", readme)
    @test occursin("derived views over that record", readme)
    @test occursin("canonical event record", readme)
    @test occursin("tree-native", readme)
    @test !occursin("tree extraction", readme)
    @test !occursin("complete event semantics", readme)

    @test occursin("tree-native representation", events_src)
    @test !occursin("transmission tree API", events_src)

    @test occursin("earliest host with no recorded infector", transmission_src)
    @test occursin("canonical event record", transmission_src)
    @test !occursin("seeded ancestor", transmission_src)
    @test !occursin("complete event semantics", transmission_src)

    @test occursin("derived views", ensemble_src)
    @test !occursin("recovery layers", ensemble_src)
end

@testset "public collection ergonomics" begin
    log = EventLog(
        [0.0, 0.4, 0.4, 1.2],
        [1, 2, 2, 1],
        [0, 1, 0, 0],
        [EK_Seeding, EK_Transmission, EK_Activation, EK_Removal],
    )

    @test length(log) == 4
    @test !isempty(log)
    @test firstindex(log) == 1
    @test lastindex(log) == 4
    @test eltype(EventLog) === EventRecord
    @test log[1] == EventRecord(0.0, 1, 0, EK_Seeding)
    @test log[end] == EventRecord(1.2, 1, 0, EK_Removal)
    @test first(log) == log[1]
    @test last(log) == log[end]
    @test collect(log) == [
        EventRecord(0.0, 1, 0, EK_Seeding),
        EventRecord(0.4, 2, 1, EK_Transmission),
        EventRecord(0.4, 2, 0, EK_Activation),
        EventRecord(1.2, 1, 0, EK_Removal),
    ]
    @test_throws BoundsError log[0]
    @test sprint(show, log[1]) == "EventRecord(time=0.0, host=1, kind=Seeding)"
    @test sprint(show, EventRecord(1.1926, 1, 0, EK_Activation)) == "EventRecord(time=1.193, host=1, kind=Activation)"
    @test sprint(show, EventRecord(1.3121, 1, 0, EK_SerialSampling)) == "EventRecord(time=1.312, host=1, kind=SerialSampling)"
    @test sprint(show, EventRecord(2.4812, 4, 1, EK_Transmission)) == "EventRecord(time=2.481, host=4, infector=1, kind=Transmission)"

    empty_log = EventLog(Float64[], Int[], Int[], EventKind[])
    @test isempty(empty_log)
    @test collect(empty_log) == EventRecord[]

    traj = StateCountTrajectory(
        [0.0, 0.4, 0.4, 1.2],
        [9, 8, 8, 8],
        [0, 1, 0, 0],
        [1, 1, 2, 1],
        [0, 0, 0, 1],
    )
    @test length(traj) == 4
    @test firstindex(traj) == 1
    @test lastindex(traj) == 4
    @test eltype(StateCountTrajectory) === StateCountPoint
    @test traj[2] == StateCountPoint(0.4, 8, 1, 1, 0)
    @test traj[3].time == 0.4
    @test first(traj) == StateCountPoint(0.0, 9, 0, 1, 0)
    @test last(traj) == StateCountPoint(1.2, 8, 0, 1, 1)
    @test collect(traj)[2:3] == [
        StateCountPoint(0.4, 8, 1, 1, 0),
        StateCountPoint(0.4, 8, 0, 2, 0),
    ]
    @test occursin("StateCountPoint(time=0.4, S=8, E=1, I=1, R=0)", sprint(show, traj[2]))

    host_summary = HostEventSummary([2, 5], [1, 3], [0, 2], [1, 0], [1, 1])
    @test length(host_summary) == 2
    @test eltype(HostEventSummary) === HostEventRecord
    @test host_summary[1] == HostEventRecord(2, 1, 0, 1, 1)
    @test host_summary[2] == HostEventRecord(5, 3, 2, 0, 1)
    @test host_summary[1].host_id == 2
    @test host_summary[2].host_id == 5
    @test collect(host_summary) == [host_summary[1], host_summary[2]]
    @test_throws BoundsError host_summary[5]
    @test occursin("HostEventRecord(host_id=5, transmissions_caused=3, samples=2, removals=0, activations=1)", sprint(show, host_summary[2]))

    view = transmission_tree(log)
    @test length(view) == 1
    @test !isempty(view)
    @test eltype(TransmissionTreeView) === TransmissionEdge
    @test view[1] == TransmissionEdge(1, 2, 0.4)
    @test first(view) == view[1]
    @test last(view) == view[1]
    @test collect(view) == [TransmissionEdge(1, 2, 0.4)]
    @test transmission_edges(view) == [(infector=1, infectee=2, time=0.4)]
    @test occursin("TransmissionEdge(1 -> 2 at time 0.4)", sprint(show, view[1]))

    empty_view = TransmissionTreeView(Int[], Int[], Float64[])
    @test isempty(empty_view)
    @test collect(empty_view) == TransmissionEdge[]

    chain = transmission_chain(view, 2)
    @test length(chain) == 2
    @test eltype(TransmissionChain) === TransmissionChainStep
    @test chain[1] == TransmissionChainStep(1, nothing)
    @test chain[end] == TransmissionChainStep(2, 0.4)
    @test collect(chain) == [TransmissionChainStep(1, nothing), TransmissionChainStep(2, 0.4)]
    @test occursin("TransmissionChainStep(host_id=2, infection_time=0.4)", sprint(show, chain[end]))

    ensemble = EnsembleSummary(
        2,
        [1, 2],
        [1, 4],
        [0.0, 1.2],
        [0, 1],
        [0, 1],
        [1, 1],
        [0, 0],
        [0, 0],
        nothing,
    )
    @test length(ensemble) == 2
    @test eltype(EnsembleSummary) === EnsembleReplicateSummary
    @test ensemble[1] == EnsembleReplicateSummary(1, 1, 1, 0.0, 0, 0, 1, 0, 0, nothing)
    @test ensemble[end].log === nothing
    @test first(ensemble) == ensemble[1]
    @test last(ensemble) == ensemble[2]
    @test collect(ensemble) == [ensemble[1], ensemble[2]]
    @test occursin("log not retained", sprint(show, ensemble[1]))

    retained = EnsembleSummary(
        1,
        [2],
        [4],
        [1.2],
        [1],
        [1],
        [1],
        [0],
        [0],
        [log],
    )
    @test retained[1].log === log
    @test occursin("log retained", sprint(show, retained[1]))
end

@testset "aggregate summary helpers" begin
    ensemble = EnsembleSummary(
        3,
        [2, 4, 6],
        [5, 7, 9],
        [1.0, 2.0, 4.0],
        [1, 3, 5],
        [0, 1, 2],
        [1, 1, 3],
        [2, 0, 1],
        [0, 1, 1],
        nothing,
    )
    ensemble_before = (
        copy(ensemble.final_size),
        copy(ensemble.total_events),
        copy(ensemble.final_time),
        copy(ensemble.transmissions),
        copy(ensemble.activations),
        copy(ensemble.removals),
        copy(ensemble.fossilized_samples),
        copy(ensemble.serial_samples),
    )

    ensemble_stats = ensemble_aggregate_summary(ensemble)
    @test ensemble_stats isa EnsembleAggregateSummary
    @test ensemble_stats.final_size == ScalarSummary(4.0, 2.0, 6.0)
    @test ensemble_stats.final_time == ScalarSummary(7 / 3, 1.0, 4.0)
    @test ensemble_stats.transmissions == ScalarSummary(3.0, 1.0, 5.0)
    @test ensemble_stats.activations == ScalarSummary(1.0, 0.0, 2.0)
    @test ensemble_stats.removals == ScalarSummary(5 / 3, 1.0, 3.0)
    @test ensemble_stats.samples == ScalarSummary(5 / 3, 1.0, 2.0)
    @test ensemble_stats.total_samples == 5
    @test occursin("EnsembleAggregateSummary", sprint(show, ensemble_stats))
    @test occursin("final outbreak size", sprint(show, ensemble_stats))
    @test occursin("total samples: 5", sprint(show, ensemble_stats))
    @test !occursin("Per-replicate quantities", sprint(show, ensemble_stats))
    @test !occursin("Ensemble-wide totals", sprint(show, ensemble_stats))

    ensemble_show = sprint(show, ensemble)
    @test occursin("EnsembleSummary(3 replicates", ensemble_show)
    @test occursin("logs not retained", ensemble_show)

    trajs = [
        StateCountTrajectory([0.0, 1.0, 2.0], [5, 4, 4], [0, 1, 0], [1, 2, 1], [0, 0, 2]),
        StateCountTrajectory([0.0, 0.5], [3, 2], [0, 0], [1, 3], [0, 1]),
    ]
    trajs_before = [(copy(t.time), copy(t.S), copy(t.E), copy(t.I), copy(t.R)) for t in trajs]

    trajectory_stats = trajectory_aggregate_summary(trajs)
    @test trajectory_stats isa TrajectoryAggregateSummary
    @test trajectory_stats.peak_infectious == ScalarSummary(2.5, 2.0, 3.0)
    @test trajectory_stats.peak_infectious_time == ScalarSummary(0.75, 0.5, 1.0)
    @test trajectory_stats.final_removed == ScalarSummary(1.5, 1.0, 2.0)
    @test trajectory_stats.final_time == ScalarSummary(1.25, 0.5, 2.0)
    @test occursin("TrajectoryAggregateSummary", sprint(show, trajectory_stats))
    @test occursin("first time of peak infectious count", sprint(show, trajectory_stats))
    @test !occursin("no common time axis", sprint(show, trajectory_stats))

    traj_show = sprint(show, trajs[1])
    @test occursin("StateCountTrajectory(3 points", traj_show)
    @test occursin("time range 0.0 to 2.0", traj_show)

    host_summaries = [
        HostEventSummary([1, 2], [2, 0], [1, 3], [0, 1], [1, 1]),
        HostEventSummary([4], [3], [0], [2], [1]),
    ]
    host_before = [
        (copy(s.host_id), copy(s.transmissions_caused), copy(s.samples), copy(s.removals), copy(s.activations))
        for s in host_summaries
    ]

    host_stats = host_aggregate_summary(host_summaries)
    @test host_stats isa HostAggregateSummary
    @test host_stats.observed_hosts == ScalarSummary(1.5, 1.0, 2.0)
    @test host_stats.mean_transmissions_per_host == ScalarSummary(2.0, 1.0, 3.0)
    @test host_stats.max_transmissions_per_host == ScalarSummary(2.5, 2.0, 3.0)
    @test host_stats.mean_samples_per_host == ScalarSummary(1.0, 0.0, 2.0)
    @test host_stats.mean_removals_per_host == ScalarSummary(1.25, 0.5, 2.0)
    @test host_stats.mean_activations_per_host == ScalarSummary(1.0, 1.0, 1.0)
    @test occursin("HostAggregateSummary", sprint(show, host_stats))
    @test occursin("observed hosts per replicate", sprint(show, host_stats))
    @test !occursin("host IDs are not matched", sprint(show, host_stats))

    host_show = sprint(show, host_summaries[1])
    @test occursin("HostEventSummary(2 observed hosts", host_show)
    @test occursin("host id range 1 to 2", host_show)

    @test ensemble_before == (
        ensemble.final_size,
        ensemble.total_events,
        ensemble.final_time,
        ensemble.transmissions,
        ensemble.activations,
        ensemble.removals,
        ensemble.fossilized_samples,
        ensemble.serial_samples,
    )
    @test [(t.time, t.S, t.E, t.I, t.R) for t in trajs] == trajs_before
    @test [(s.host_id, s.transmissions_caused, s.samples, s.removals, s.activations) for s in host_summaries] == host_before

    @test_throws ArgumentError EpiSim.scalar_summary(Int[])
    @test_throws ArgumentError ensemble_aggregate_summary(EnsembleSummary(0, Int[], Int[], Float64[], Int[], Int[], Int[], Int[], Int[], nothing))
    @test_throws ArgumentError trajectory_aggregate_summary(StateCountTrajectory[])
    @test_throws ArgumentError host_aggregate_summary(HostEventSummary[])
end

@testset "transmission tree and chain views" begin
    el = EventLog(
        [0.0, 0.0, 0.5, 0.8, 1.0, 1.3],
        [1, 2, 3, 4, 3, 4],
        [0, 0, 1, 3, 0, 0],
        [EK_Seeding, EK_Seeding, EK_Transmission, EK_Transmission, EK_Activation, EK_Removal],
    )
    original = (copy(el.time), copy(el.host), copy(el.infector), copy(el.kind))

    view = transmission_tree(el)
    @test view isa TransmissionTreeView
    @test length(view) == 2
    @test view.infector == [1, 3]
    @test view.infectee == [3, 4]
    @test view.time == [0.5, 0.8]
    @test transmission_edges(view) == [
        (infector=1, infectee=3, time=0.5),
        (infector=3, infectee=4, time=0.8),
    ]
    @test transmission_edges(el) == transmission_edges(view)

    view_show = sprint(show, view)
    @test occursin("TransmissionTreeView", view_show)
    @test occursin("2 edges", view_show)
    @test !occursin("not a TreeSim tree", view_show)

    chain = transmission_chain(view, 4)
    @test chain isa TransmissionChain
    @test length(chain) == 3
    @test chain.host_id == 4
    @test chain.host_path == [1, 3, 4]
    @test chain.infection_time == Union{Nothing,Float64}[nothing, 0.5, 0.8]
    @test transmission_chain(el, 4).host_path == chain.host_path

    seed_chain = transmission_chain(view, 2)
    @test seed_chain.host_path == [2]
    @test seed_chain.infection_time == Union{Nothing,Float64}[nothing]

    chain_show = sprint(show, chain)
    @test occursin("TransmissionChain", chain_show)
    @test occursin("host 4", chain_show)
    @test occursin("source 1 -> host 4", chain_show)

    empty_view = transmission_tree(EventLog([0.0], [1], [0], [EK_Seeding]))
    @test length(empty_view) == 0
    @test transmission_edges(empty_view) == NamedTuple{(:infector, :infectee, :time),Tuple{Int,Int,Float64}}[]
    @test occursin("no transmission events", sprint(show, empty_view))

    @test_throws ArgumentError transmission_chain(view, 0)
    cyclic = TransmissionTreeView([2, 1], [1, 2], [0.1, 0.2])
    @test_throws ArgumentError transmission_chain(cyclic, 1)
    @test (el.time, el.host, el.infector, el.kind) == original
end

@testset "visualization support summaries" begin
    traj1 = StateCountTrajectory(
        [0.0, 0.5, 0.5, 1.0],
        [5, 4, 4, 4],
        [0, 1, 0, 0],
        [1, 1, 2, 1],
        [0, 0, 0, 1],
    )
    traj2 = StateCountTrajectory(
        [0.0, 0.25, 0.75],
        [3, 2, 2],
        [0, 1, 0],
        [1, 1, 2],
        [0, 0, 0],
    )
    traj1_before = (copy(traj1.time), copy(traj1.S), copy(traj1.E), copy(traj1.I), copy(traj1.R))

    time, infectious = trajectory_series(traj1, :I)
    @test time === traj1.time
    @test infectious === traj1.I
    @test time == [0.0, 0.5, 0.5, 1.0]
    @test infectious == [1, 1, 2, 1]

    all_series = trajectory_series(traj1)
    @test all_series.S == (traj1.time, traj1.S)
    @test all_series.E == (traj1.time, traj1.E)
    @test all_series.I == (traj1.time, traj1.I)
    @test all_series.R == (traj1.time, traj1.R)
    @test_throws ArgumentError trajectory_series(traj1, :X)

    trajs = [traj1, traj2]
    @test trajectory_final_sizes(trajs) == [2, 2]
    @test trajectory_final_times(trajs) == [1.0, 0.75]
    @test peak_infectious(trajs) == [2, 2]
    @test peak_infectious_time(trajs) == [0.5, 0.75]

    host_summary = HostEventSummary([1, 3], [2, 0], [1, 2], [0, 1], [1, 0])
    host_summary_before = (
        copy(host_summary.host_id),
        copy(host_summary.transmissions_caused),
        copy(host_summary.samples),
        copy(host_summary.removals),
        copy(host_summary.activations),
    )
    hosts, transmissions = host_series(host_summary, :transmissions)
    @test hosts === host_summary.host_id
    @test transmissions === host_summary.transmissions_caused
    @test host_series(host_summary, :transmissions_caused) == (host_summary.host_id, host_summary.transmissions_caused)
    @test host_series(host_summary, :samples) == (host_summary.host_id, host_summary.samples)
    @test host_series(host_summary, :removals) == (host_summary.host_id, host_summary.removals)
    @test host_series(host_summary, :activations) == (host_summary.host_id, host_summary.activations)
    @test_throws ArgumentError host_series(host_summary, :secondary_cases)

    @test (traj1.time, traj1.S, traj1.E, traj1.I, traj1.R) == traj1_before
    @test (host_summary.host_id, host_summary.transmissions_caused, host_summary.samples,
           host_summary.removals, host_summary.activations) == host_summary_before
end

@testset "event-log processing helpers" begin
    el = EventLog(
        [0.0, 0.0, 0.5, 0.7, 1.0, 1.2, 1.2, 1.5],
        [1, 2, 3, 1, 3, 2, 3, 1],
        [0, 0, 1, 0, 0, 0, 0, 0],
        [EK_Seeding, EK_Seeding, EK_Transmission, EK_FossilizedSampling,
         EK_Activation, EK_SerialSampling, EK_Removal, EK_Activation],
    )
    original = (copy(el.time), copy(el.host), copy(el.infector), copy(el.kind))

    @test total_events(el) == 8
    @test event_count(el, EK_Seeding) == 2
    @test event_count(el, :seeding) == 2
    @test event_count(el, :seedings) == 2
    @test event_count(el, :serial_sampling) == 1
    @test event_count(el, "serial_sampling") == 1
    @test event_count(el, EK_Activation) == 2
    @test event_indices(el, EK_Activation) == [5, 8]
    @test event_indices(el, :activation) == [5, 8]
    @test event_times(el, EK_Activation) == [1.0, 1.5]
    @test event_times(el, :activation) == [1.0, 1.5]
    @test event_times(el, "activation") == [1.0, 1.5]
    @test event_indices(el, EpiSim.EK_None) == Int[]
    @test event_times(el, EpiSim.EK_None) == Float64[]

    @test first_event_time(el, EK_Activation) == 1.0
    @test first_event_time(el, :activation) == 1.0
    @test last_event_time(el, EK_Activation) == 1.5
    @test last_event_time(el, :activation) == 1.5
    @test first_event_time(el, EpiSim.EK_None) === nothing
    @test last_event_time(el, EpiSim.EK_None) === nothing
    @test has_event_kind(el, EK_SerialSampling)
    @test has_event_kind(el, :serial_sampling)
    @test !has_event_kind(el, EpiSim.EK_None)
    @test event_kind(:transmission) == EK_Transmission
    @test event_kind(:transmissions) == EK_Transmission
    @test event_kind(:fossilized_sampling) == EK_FossilizedSampling
    @test event_kind(:fossilised_sampling) == EK_FossilizedSampling
    @test event_kind(" fossilised_sampling ") == EK_FossilizedSampling
    @test event_kind(EK_Removal) == EK_Removal
    @test_throws ArgumentError event_kind(:seed)
    @test_throws ArgumentError event_kind(:none)
    @test_throws ArgumentError event_times(el, :secondary_case)
    invalid_message = try
        event_kind(:secondary_case)
    catch err
        sprint(showerror, err)
    end
    @test occursin("unknown event kind :secondary_case", invalid_message)
    @test occursin(":transmission", invalid_message)
    @test occursin(":fossilised_sampling", invalid_message)
    @test occursin("plural forms are also accepted", invalid_message)

    @test observed_hosts(el) == [1, 2, 3]
    @test observed_hosts(EventLog([0.5], [3], [9], [EK_Transmission]); include_sources=false) == [3]
    @test distinct_host_count(el) == 3

    summary = host_event_summary(el)
    @test summary isa HostEventSummary
    @test length(summary) == 3
    @test summary.host_id == [1, 2, 3]
    @test summary.transmissions_caused == [1, 0, 0]
    @test summary.samples == [1, 1, 0]
    @test summary.removals == [0, 1, 1]
    @test summary.activations == [1, 0, 1]

    empty = EventLog(Float64[], Int[], Int[], EventKind[])
    empty_summary = host_event_summary(empty)
    @test total_events(empty) == 0
    @test observed_hosts(empty) == Int[]
    @test distinct_host_count(empty) == 0
    @test event_indices(empty, EK_Transmission) == Int[]
    @test event_times(empty, EK_Transmission) == Float64[]
    @test first_event_time(empty, EK_Transmission) === nothing
    @test last_event_time(empty, EK_Transmission) === nothing
    @test empty_summary.host_id == Int[]
    @test empty_summary.transmissions_caused == Int[]
    @test empty_summary.samples == Int[]
    @test empty_summary.removals == Int[]
    @test empty_summary.activations == Int[]

    @test (el.time, el.host, el.infector, el.kind) == original
end

@testset "event-time state-count recovery" begin
    el = EventLog(
        [0.0, 0.0, 0.5, 1.0, 1.0, 1.5],
        [1, 2, 3, 1, 3, 2],
        [0, 0, 1, 0, 0, 0],
        [EK_Seeding, EK_Seeding, EK_Transmission, EK_Activation, EK_Activation, EK_Removal],
    )
    original = (copy(el.time), copy(el.host), copy(el.infector), copy(el.kind))

    traj = event_time_state_counts(el; S0=3, E0=1, I0=1)

    @test traj isa StateCountTrajectory
    @test length(traj) == length(el) + 1
    @test traj.time == [0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.5]
    @test traj.S == [3, 3, 3, 2, 2, 2, 2]
    @test traj.E == [1, 1, 1, 2, 1, 0, 0]
    @test traj.I == [1, 1, 1, 1, 2, 3, 2]
    @test traj.R == [0, 0, 0, 0, 0, 0, 1]
    @test (el.time, el.host, el.infector, el.kind) == original

    empty = EventLog(Float64[], Int[], Int[], EventKind[])
    empty_traj = event_time_state_counts(empty; S0=5, E0=0, I0=0, R0=0)
    @test length(empty_traj) == 1
    @test empty_traj.time == [0.0]
    @test empty_traj.S == [5]
    @test empty_traj.E == [0]
    @test empty_traj.I == [0]
    @test empty_traj.R == [0]

    sampled = EventLog(
        [0.0, 0.2, 0.2, 0.3],
        [1, 1, 1, 1],
        [0, 0, 0, 0],
        [EK_Seeding, EK_FossilizedSampling, EK_FossilizedSampling, EK_SerialSampling],
    )
    sampled_traj = event_time_state_counts(sampled; S0=0, E0=0, I0=1)
    @test sampled_traj.time == [0.0, 0.0, 0.2, 0.2, 0.3]
    @test sampled_traj.I == [1, 1, 1, 1, 0]
    @test sampled_traj.R == [0, 0, 0, 0, 1]

    simulated = gillespie(Random.MersenneTwister(901), 5, 0, 1, 0.7, 1.0, 1.0, 0.2, 0.5)
    simulated_traj = event_time_state_counts(simulated; S0=5, E0=0, I0=1)
    @test simulated_traj.S[end] == 5 - event_count(simulated, EK_Transmission)
    @test simulated_traj.E[end] == event_count(simulated, EK_Transmission) - event_count(simulated, EK_Activation)
    @test simulated_traj.I[end] == 1 + event_count(simulated, EK_Activation) -
                                   event_count(simulated, EK_Removal) -
                                   event_count(simulated, EK_SerialSampling)
    @test simulated_traj.R[end] == event_count(simulated, EK_Removal) +
                                   event_count(simulated, EK_SerialSampling)

    @test_throws ArgumentError event_time_state_counts(el; S0=-1, E0=1, I0=1)
    @test_throws ArgumentError event_time_state_counts(el; S0=0, E0=0, I0=0)
end

@testset "ensemble derived helpers over retained logs" begin
    retained = run_ensemble(
        rng -> gillespie(rng, 5, 0, 1, 0.4, 1.0, 1.0, 0.0, 0.0),
        3;
        rng=Random.MersenneTwister(7701),
        retain_logs=true,
    )
    logs_before = [(copy(log.time), copy(log.host), copy(log.infector), copy(log.kind)) for log in retained.logs]
    summary_before = (
        copy(retained.final_size),
        copy(retained.total_events),
        copy(retained.final_time),
        copy(retained.transmissions),
        copy(retained.activations),
        copy(retained.removals),
        copy(retained.fossilized_samples),
        copy(retained.serial_samples),
    )

    trajectories = ensemble_state_trajectories(retained; S0=5, E0=0, I0=1)
    @test length(trajectories) == 3
    @test all(traj -> traj isa StateCountTrajectory, trajectories)
    @test trajectories[1].time == event_time_state_counts(retained.logs[1]; S0=5, E0=0, I0=1).time
    @test trajectories[1].S == event_time_state_counts(retained.logs[1]; S0=5, E0=0, I0=1).S
    @test [last(traj.R) for traj in trajectories] == retained.removals .+ retained.serial_samples

    host_summaries = ensemble_host_event_summaries(retained)
    @test length(host_summaries) == 3
    @test all(summary -> summary isa HostEventSummary, host_summaries)
    @test host_summaries[1].host_id == host_event_summary(retained.logs[1]).host_id
    @test host_summaries[1].removals == host_event_summary(retained.logs[1]).removals

    @test [(log.time, log.host, log.infector, log.kind) for log in retained.logs] == logs_before
    @test summary_before == (
        retained.final_size,
        retained.total_events,
        retained.final_time,
        retained.transmissions,
        retained.activations,
        retained.removals,
        retained.fossilized_samples,
        retained.serial_samples,
    )

    not_retained = run_ensemble(
        rng -> gillespie(rng, 5, 0, 1, 0.4, 1.0, 1.0, 0.0, 0.0),
        1;
        rng=Random.MersenneTwister(7702),
    )
    @test_throws ArgumentError ensemble_state_trajectories(not_retained; S0=5, E0=0, I0=1)
    @test_throws ArgumentError ensemble_host_event_summaries(not_retained)

    empty = run_ensemble(
        rng -> EventLog(Float64[], Int[], Int[], EventKind[]),
        0;
        retain_logs=true,
    )
    @test ensemble_state_trajectories(empty; S0=5, E0=0, I0=0) == StateCountTrajectory[]
    @test ensemble_host_event_summaries(empty) == HostEventSummary[]
end

@testset "compact ensemble summaries" begin
    summary = run_ensemble(
        rng -> gillespie(rng, 5, 0, 1, 0.0, 1.0, 1.0, 0.0, 0.0),
        4;
        rng=Random.MersenneTwister(3301),
    )

    @test summary isa EnsembleSummary
    @test length(summary) == 4
    @test summary.nrep == 4
    @test summary.logs === nothing
    @test summary.final_size == fill(1, 4)
    @test summary.transmissions == fill(0, 4)
    @test summary.removals == fill(1, 4)
    @test summary.serial_samples == fill(0, 4)
    @test summary.fossilized_samples == fill(0, 4)
    @test summary.total_events == fill(2, 4)
    @test all(>(0.0), summary.final_time)
    @test mean_final_size(summary) == 1.0
    @test attack_rate(summary, 6) == 1 / 6
    summary_show = sprint(show, summary)
    @test occursin("EnsembleSummary(4 replicates", summary_show)
    @test occursin("logs not retained", summary_show)

    retained = run_ensemble(
        rng -> gillespie(rng, 5, 0, 1, 0.0, 1.0, 1.0, 0.0, 0.0),
        2;
        rng=Random.MersenneTwister(3302),
        retain_logs=true,
    )

    @test retained.logs !== nothing
    @test length(retained.logs) == 2
    @test all(log -> log isa EventLog, retained.logs)
    @test retained.final_size == final_size.(retained.logs)
    @test retained.final_time == final_time.(retained.logs)

    again = run_ensemble(
        rng -> gillespie(rng, 5, 0, 1, 0.0, 1.0, 1.0, 0.0, 0.0),
        2;
        rng=Random.MersenneTwister(3302),
        retain_logs=true,
    )
    @test retained.final_time == again.final_time
    @test retained.logs[1].time == again.logs[1].time
    @test occursin("logs retained", sprint(show, retained))

    Random.seed!(4401)
    expected_global_draw = rand()
    Random.seed!(4401)
    run_ensemble(
        rng -> gillespie(rng, 5, 0, 1, 0.0, 1.0, 1.0, 0.0, 0.0),
        2;
        rng=Random.MersenneTwister(4402),
    )
    @test rand() == expected_global_draw

    @test_throws ArgumentError run_ensemble(rng -> nothing, 1)
    @test_throws ArgumentError run_ensemble(rng -> EventLog(1), -1)
    @test_throws ArgumentError attack_rate(summary, 0)
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
    eventlog_show = sprint(show, el)
    @test occursin("EventLog(3 events", eventlog_show)
    @test occursin("time range 0.0 to 0.0", eventlog_show)
    @test !occursin("canonical simulation event record", eventlog_show)

    EpiSim.update_event_log!(el, 1.5, 4, 2, EK_Transmission)

    assert_eventlog_invariants(el)
    @test length(el) == 4
    @test el.time[end] == 1.5
    @test el.host[end] == 4
    @test el.infector[end] == 2
    @test el.kind[end] == EK_Transmission
    @test occursin("time range 0.0 to 1.5", sprint(show, el))
end

@testset "EventLog semantic validation" begin
    valid = EventLog(
        [0.0, 1.0, 1.5, 2.0, 2.5],
        [1, 2, 2, 2, 1],
        [0, 1, 0, 0, 0],
        [EK_Seeding, EK_Transmission, EK_Activation, EK_FossilizedSampling, EK_Removal],
    )
    @test validate_event_log(valid; population_size=2)

    @test !validate_event_log(EventLog([0.0, 0.5, 0.4], [1, 2, 2], [0, 1, 0],
                                       [EK_Seeding, EK_Transmission, EK_Activation]); throw=false)
    @test !validate_event_log(EventLog([0.0], [0], [0], [EK_Seeding]); throw=false)
    @test !validate_event_log(EventLog([0.0, 1.0], [1, 2], [0, 0],
                                       [EK_Seeding, EK_Transmission]); throw=false)
    @test !validate_event_log(EventLog([0.0, 1.0], [1, 2], [0, 2],
                                       [EK_Seeding, EK_Activation]); throw=false)
    @test !validate_event_log(EventLog([0.0, 1.0], [1, 1], [0, 2],
                                       [EK_Seeding, EK_Transmission]); throw=false)
    @test !validate_event_log(EventLog([0.0], [1], [0], [EpiSim.EK_None]); throw=false)
    @test_throws ErrorException validate_event_log(EventLog([Inf], [1], [0], [EK_Seeding]))
end

@testset "sellke smoke" begin
    Random.seed!(1234)
    el = sellke(4, 0, 1, 0.8, 0.25, 0.75, 0.2, 0.5)

    @test el isa EventLog
    assert_eventlog_invariants(el)
    @test !isempty(el.time)
    @test all(isfinite, el.time)
    @test all(el.time .>= 0.0)
end

@testset "gillespie smoke" begin
    Random.seed!(1234)
    el = gillespie(4, 0, 1, 0.8, 1.0, 1.0, 0.5, 0.5)

    @test el isa EventLog
    assert_eventlog_invariants(el)
    @test !isempty(el.time)
    @test all(isfinite, el.time)
    @test all(el.time .>= 0.0)
end

@testset "fixed-seed reproducibility" begin
    Random.seed!(202402)
    s1 = sellke(8, 1, 2, 0.5, 0.25, 0.75, 0.2, 0.4)
    Random.seed!(202402)
    s2 = sellke(8, 1, 2, 0.5, 0.25, 0.75, 0.2, 0.4)
    @test s1.time == s2.time
    @test s1.host == s2.host
    @test s1.infector == s2.infector
    @test s1.kind == s2.kind

    Random.seed!(202403)
    g1 = gillespie(8, 1, 2, 0.5, 1.2, 0.7, 0.3, 0.4)
    Random.seed!(202403)
    g2 = gillespie(8, 1, 2, 0.5, 1.2, 0.7, 0.3, 0.4)
    @test g1.time == g2.time
    @test g1.host == g2.host
    @test g1.infector == g2.infector
    @test g1.kind == g2.kind
end

@testset "simple finite-population boundary behaviour" begin
    Random.seed!(77)
    el = gillespie(5, 0, 1, 0.0, 1.0, 1.0, 0.0, 0.0)

    @test validate_event_log(el; population_size=6)
    @test count(==(EK_Seeding), el.kind) == 1
    @test count(==(EK_Transmission), el.kind) == 0
    @test count(==(EK_Removal), el.kind) == 1
    @test el.kind[end] == EK_Removal
end

@testset "shared zero-transmission limiting behaviour" begin
    Random.seed!(1701)
    sellke_log = sellke(5, 1, 2, 0.0, 1.0, 1.0, 1000.0, 0.0)
    Random.seed!(1701)
    gillespie_log = gillespie(5, 1, 2, 0.0, 1.0, 1.0, 0.0, 0.0)

    for el in (sellke_log, gillespie_log)
        @test validate_event_log(el; population_size=8)
        @test event_count(el, EK_Transmission) == 0
        @test final_size(el) == 3
        @test event_count(el, EK_Activation) == 1
        @test event_count(el, EK_Removal) == 3
    end
end

@testset "cross-engine ensemble agreement" begin
    nrep = 250
    sellke_summary = mean_summary(() -> sellke(10, 0, 1, 1.1, 1.0, 1.0, 1000.0, 0.0), nrep; seed=8101)
    gillespie_summary = mean_summary(() -> gillespie(10, 0, 1, 1.1, 1.0, 1.0, 0.0, 0.0), nrep; seed=8102)

    @test abs(sellke_summary.mean_final_size - gillespie_summary.mean_final_size) ≤ 1.0
    @test abs(sellke_summary.mean_transmissions - gillespie_summary.mean_transmissions) ≤ 1.0
end

@testset "ensemble transmission monotonicity" begin
    nrep = 250
    low = mean_summary(() -> gillespie(15, 0, 1, 0.25, 1.0, 1.0, 0.0, 0.0), nrep; seed=9101)
    high = mean_summary(() -> gillespie(15, 0, 1, 1.5, 1.0, 1.0, 0.0, 0.0), nrep; seed=9102)

    @test high.mean_final_size > low.mean_final_size + 2.0
    @test high.mean_transmissions > low.mean_transmissions + 2.0
end

@testset "exact finite-state SEIR transient validation" begin
    population_size = 4
    S0, E0, I0 = 3, 0, 1
    β, α, γ = 1.0, 1.3, 0.8
    t_observe = 1.2

    states, probabilities = exact_seir_distribution_uniformization(population_size, S0, E0, I0, β, α, γ, t_observe)
    exact_infectious_distribution = infectious_count_distribution(states, probabilities, population_size)
    exact_mean_infected = sum(probabilities[i] * states[i][3] for i in eachindex(states))
    exact_extinction_probability = sum(probabilities[i] for i in eachindex(states) if states[i][3] == 0)

    nrep = 6000
    no_sampling = DiscreteNonParametric([Inf], [1.0])
    engine_distributions = (
        gillespie = empirical_infectious_count_distribution(
            () -> gillespie(S0, E0, I0, β, α, γ, 0.0, 0.0),
            nrep, t_observe, population_size, S0, E0, I0; seed=12031),
        sellke = empirical_infectious_count_distribution(
            () -> sellke(S0, E0, I0, β, 1 / α, 1 / γ, no_sampling, 0.0),
            nrep, t_observe, population_size, S0, E0, I0; seed=12032),
    )

    for estimated_distribution in engine_distributions
        estimated_mean_infected = sum((k - 1) * estimated_distribution[k] for k in eachindex(estimated_distribution))
        estimated_extinction_probability = estimated_distribution[1]

        @test total_variation_distance(estimated_distribution, exact_infectious_distribution) ≤ 0.04
        @test maximum(abs.(estimated_distribution .- exact_infectious_distribution)) ≤ 0.04
        @test abs(estimated_mean_infected - exact_mean_infected) ≤ 0.08
        @test abs(estimated_extinction_probability - exact_extinction_probability) ≤ 0.04
    end
end
