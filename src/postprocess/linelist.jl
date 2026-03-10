using DataFrames

struct LineList{H<:AbstractDataFrame,S<:AbstractDataFrame}
    hosts::H
    samples::S
end


function LineList(log::EventLog)

    N = maximum(log.host)

    infection_time  = Vector{Union{Missing,Time}}(missing, N)
    activation_time = Vector{Union{Missing,Time}}(missing, N)
    removal_time    = Vector{Union{Missing,Time}}(missing, N)
    infector        = Vector{Union{Missing,HostID}}(missing, N)

    n_samples        = zeros(Int, N)
    first_sample     = Vector{Union{Missing,Time}}(missing, N)
    last_sample      = Vector{Union{Missing,Time}}(missing, N)

    sample_host = HostID[]
    sample_time = Time[]
    sample_kind = EventKind[]

    for ev in log

        id = ev.host
        t  = ev.t

        if ev.kind == Transmission
            infection_time[id] = t
            infector[id] = ev.src

        elseif ev.kind == Seeding
            infection_time[id] = t
            infector[id] = 0

        elseif ev.kind == Activation
            activation_time[id] = t

        elseif ev.kind == Removal
            removal_time[id] = t

        # TODO: Need to allow SerialSampling to also remove host
        elseif ev.kind == SerialSampling || ev.kind == FossilisedSampling

            push!(sample_host, id)
            push!(sample_time, t)
            push!(sample_kind, ev.kind)

            n_samples[id] += 1

            if ismissing(first_sample[id])
                first_sample[id] = t
                last_sample[id]  = t
            else
                first_sample[id] = min(first_sample[id], t)
                last_sample[id]  = max(last_sample[id], t)
            end

            if ev.kind == SerialSampling
                # Assume serial sampling also removes host
                removal_time[id] = t
            end
        end
    end

    hosts = DataFrame(
        id = 1:N,
        infection_time = infection_time,
        activation_time = activation_time,
        removal_time = removal_time,
        infector = infector,
        n_samples = n_samples,
        first_sample_time = first_sample,
        last_sample_time = last_sample
    )

    samples = DataFrame(
        host_id = sample_host,
        sample_time = sample_time,
        sample_kind = sample_kind
    )

    return LineList(hosts, samples)
end


function secondary_cases(ll::LineList)

    df = ll.hosts
    N = nrow(df)

    counts = zeros(Int, N)

    for row in eachrow(df)

        if !ismissing(row.infector) && row.infector > 0
            counts[row.infector] += 1
        end

    end

    return counts
end


function infectious_periods(ll::LineList)
    df = ll.hosts
    periods = Float64[]

    for row in eachrow(df)
        if !ismissing(row.infection_time) &&
           !ismissing(row.removal_time)
        
           start_time = !ismissing(row.activation_time) ? row.activation_time : row.infection_time
           push!(periods, row.removal_time - start_time)
        end
    end

    return periods
end


function generation_intervals(ll::LineList)
    df = ll.hosts
    gi = Float64[]

    for row in eachrow(df)
        if !ismissing(row.infector) && row.infector > 0
            parent_inf_time = df.infection_time[row.infector]
            if !ismissing(parent_inf_time)
                push!(gi, row.infection_time - parent_inf_time)
            end
        end
    end

    return gi
end


function serial_intervals(ll::LineList)
    df = ll.hosts
    si = Float64[]

    for row in eachrow(df)
        if !ismissing(row.infector) && row.infector > 0
            parent_act = df.activation_time[row.infector]
            if !ismissing(parent_act)
                push!(si, row.activation_time - parent_act)
            end
        end
    end

    return si
end


function attack_rate(ll::LineList)
    df = ll.hosts
    infected = count(!ismissing, df.infection_time)
    return infected / nrow(df)
end


function empirical_R(ll::LineList)
    sc = secondary_cases(ll)
    return mean(sc)
end


function empirical_R_completed(ll::LineList)

    df = ll.hosts
    counts = secondary_cases(ll)

    completed = .!ismissing.(df.removal_time)

    return mean(counts[completed])
end


function empirical_sampling_proportion(ll::LineList)

    df = ll.hosts

    sampled = 0
    removed = 0

    for row in eachrow(df)

        if !ismissing(row.removal_time)

            removed += 1

            if row.n_samples > 0
                sampled += 1
            end

        end
    end

    return sampled / removed
end


function susceptible_fraction(log::EventLog, S₀::Int, E₀::Int, I₀::Int)

    S = S₀
    N = S₀ + E₀ + I₀

    frac = Dict{Int,Float64}()

    for ev in log

        if ev.kind == Transmission || ev.kind == Seeding
            frac[ev.host] = S/N
            S -= 1
        end

    end

    return frac
end


function normalized_offspring(ll::LineList, log::EventLog, S₀::Int, E₀::Int, I₀::Int)

    counts = secondary_cases(ll)
    sfrac = susceptible_fraction(log, S₀, E₀, I₀)

    N = S₀ + E₀ + I₀

    norm = Float64[]

    for (i,c) in enumerate(counts)

        if haskey(sfrac,i)
            push!(norm, c / sfrac[i])
        end

    end

    return norm
end