
@enum EventKind::UInt8 begin
    EK_None = 0
    EK_Seeding = 1
    EK_Transmission = 2
    EK_Recovery = 3
    EK_Activation = 4
    EK_FossilizedSampling = 5
    EK_SerialSampling = 6
end


const EVENTKIND_LABELS = Dict(
    EK_None              => "None",
    EK_Seeding           => "Seeding",
    EK_Transmission      => "Transmission",
    EK_Recovery          => "Recovery",
    EK_Activation        => "Activation",
    EK_FossilizedSampling => "FossilizedSampling",
    EK_SerialSampling    => "SerialSampling"
)


function Base.show(io::IO, x::EventKind)
    print(io, EVENTKIND_LABELS[x])
end


struct EventLog
    time::Vector{Float64}
    host::Vector{Int}
    infector::Vector{Int}  # only for Transmission events; 0 otherwise
    kind::Vector{EventKind} # 0: None, 1: Seeding, 2: Transmission, 3: Recovery, 4: Activation, 5: FossilizedSampling, 6: SerialSampling
end


Base.length(el::EventLog) = length(el.time)

# Write iterators for time, host, infector, kind
eachtime(el::EventLog) = (el.time[i] for i in 1:length(el))
eachhost(el::EventLog) = (el.host[i] for i in 1:length(el))
eachinfector(el::EventLog) = (el.infector[i] for i in 1:length(el))
eachkind(el::EventLog) = (el.kind[i] for i in 1:length(el))