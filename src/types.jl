# Type tree and enums for TrackMatcher

## Define type tree of abstract types
abstract type DataSet{T<:AbstractFloat} end
abstract type MeasuredSet{T<:AbstractFloat} <: DataSet{T} end
abstract type PrimarySet{T<:AbstractFloat} <: MeasuredSet{T} end
abstract type PrimaryTrack{T<:AbstractFloat} <: PrimarySet{T} end
abstract type FlightTrack{T<:AbstractFloat} <: PrimaryTrack{T} end
abstract type CloudTrack{T<:AbstractFloat} <: PrimaryTrack{T} end
abstract type SecondarySet{T<:AbstractFloat} <: MeasuredSet{T} end
abstract type SecondaryTrack{T<:AbstractFloat} <: SecondarySet{T} end
abstract type SatTrack{T<:AbstractFloat} <: SecondaryTrack{T} end
abstract type ObservationSet{T<:AbstractFloat} <: MeasuredSet{T} end
abstract type ComputedSet{T<:AbstractFloat} <: DataSet{T} end
abstract type Intersection{T<:AbstractFloat} <: ComputedSet{T} end


## Define Enums for aerosol and cloud classification

"""
# Enum SkyCondition

- `invalid (0)` (bad or missing data)
- `clear (1)`
- `cloud (2)`
- `aerosol (3)`
- `stratospheric (4)` (polar stratospheric cloud (PSC) or stratospheric aerosol)
- `surface (5)`
- `subsurface (6)`
- `no_signal (7)` (totally attenuated)
"""
@enum SkyCondition::UInt16 begin
    invalid
    clear
    cloud
    aerosol
    stratospheric
    surface
    subsurface
    no_signal
end


"""
# Enum AerosolType

- `undetermined (0)` (subtype not determined)
- `marine (1)` (clean marine)
- `dust (2)`
- `polluted (3)` (polluted continental)
- `remote (4)` (remote continental)
- `polluted_dust (5)`
- `smoke (6)`
- `other (7)`
"""
@enum AerosolType::UInt16 begin
    undetermined
    marine
    dust
    polluted
    remote
    polluted_dust
    smoke
    other
end


"""
# Enum CloudType

- `low_trasparent (0)` (low overcast, transparent)
- `low_opaque (1)` (low overcast, opaque)
- `transition_sc (2)` (transition stratocumulus)
- `cu (3)` (low, broken cumulus)
- `ac (4)` (altocumulus, transparent)
- `as (5)` (altostratus, opaque)
- `ci (6)` (cirrus, transparent)
- `cb (7)` (cumulonimbus, opaque)
"""
@enum CloudType::UInt16 begin
    low_trasparent
    low_opaque
    transition_sc
    cu
    ac
    as
    ci
    cb
end

@enum SatProduct::UInt16 begin
    unknown
    clay
    cpro
    apro
    alay
end
