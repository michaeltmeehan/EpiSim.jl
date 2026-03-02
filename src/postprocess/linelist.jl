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

    counts = Dict{HostID,Int}()

    for row in eachrow(ll.hosts)
        if !ismissing(row.infection_time)
            counts[row.id] = 0
            if !ismissing(row.infector) && row.infector > 0
                counts[row.infector] = get(counts, row.infector, 0) + 1
            end
        end
    end

    return collect(values(counts))
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