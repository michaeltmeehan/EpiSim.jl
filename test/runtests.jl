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
    @test event_count(el, EK_Activation) == 2
    @test event_indices(el, EK_Activation) == [5, 8]
    @test event_times(el, EK_Activation) == [1.0, 1.5]
    @test event_indices(el, EK_None) == Int[]
    @test event_times(el, EK_None) == Float64[]

    @test first_event_time(el, EK_Activation) == 1.0
    @test last_event_time(el, EK_Activation) == 1.5
    @test first_event_time(el, EK_None) === nothing
    @test last_event_time(el, EK_None) === nothing
    @test has_event_kind(el, EK_SerialSampling)
    @test !has_event_kind(el, EK_None)

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

    EpiSim.update_event_log!(el, 1.5, 4, 2, EK_Transmission)

    assert_eventlog_invariants(el)
    @test length(el) == 4
    @test el.time[end] == 1.5
    @test el.host[end] == 4
    @test el.infector[end] == 2
    @test el.kind[end] == EK_Transmission
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
    @test !validate_event_log(EventLog([0.0], [1], [0], [EK_None]); throw=false)
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
