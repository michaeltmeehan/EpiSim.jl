@inline function nsampled(log::EventLog)
    count(ev -> ev == SerialSampling || ev == FossilisedSampling,
          log.kind)
end


function final_size(log::EventLog)
    count(ev -> ev == Transmission || ev == Seeding,
          log.kind)
end


function extinction_time(log::EventLog)
    last(filter(i -> log.kind[i] == Removal || log.kind[i] == SerialSampling,
                eachindex(log.kind))) |> i -> log.t[i]
end