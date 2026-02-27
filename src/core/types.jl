# ================================
# Core type definitions
# ================================

export HostID, Time, EventKind

const HostID = Int
const Time   = Float64

@enum EventKind::UInt8 begin
    Seeding
    Transmission
    Activation
    Removal
    SerialSampling
    FossilisedSampling
end