# Import Base methods for overloading
import Base: ==, !=, isapprox

TMType = Union{PrimaryTrack,PrimarySet,SatTrack,SecondarySet,MeasuredData,Intersection,ObservationSet,MeasuredSet,DataSet}

## Helper functions for fast nested equality/approximation for arrays with missing-aware scalar semantics
@inline _nested_equal(a, b) = isequal(a, b)

@inline function _nested_equal(a::AbstractArray, b::AbstractArray)
    axes(a) == axes(b) || return false
    @inbounds for i in eachindex(a, b)
        _nested_equal(a[i], b[i]) || return false
    end
    return true
end

@inline function _nested_approx(a, b; atol::Real=0.0, rtol::Real=sqrt(eps(Float64)))
    if ismissing(a) || ismissing(b)
        return a === b
    elseif a isa Number && b isa Number
        return Base.isapprox(a, b; atol=atol, rtol=rtol)
    else
        return isequal(a, b)
    end
end

@inline function _nested_approx(a::AbstractArray, b::AbstractArray; atol::Real=0.0, rtol::Real=sqrt(eps(Float64)))
    axes(a) == axes(b) || return false
    @inbounds for i in eachindex(a, b)
        _nested_approx(a[i], b[i]; atol=atol, rtol=rtol) || return false
    end
    return true
end


## Base overloads

"""
    ==(a::T, b::T) where {T<:TMType} -> Bool

where `TMType = Union{PrimaryTrack,PrimarySet,SatTrack,SecondarySet,MeasuredData,Intersection,ObservationSet,MeasuredSet,DataSet}`.

Fast nested equality for CPro and CLay types.

In contrast to the Base `==` method, this method provides NaN/missing-aware scalar semantics
and distinguishes between `0` and `-0` for floating-point types.
"""
function ==(a::T, b::T)::Bool where {T<:TMType}
    @inbounds for field in fieldnames(T)
        field == :metadata && continue  # ℹ skip metadata field for equality check
        _nested_equal(getproperty(a, field), getproperty(b, field)) || return false
    end
    return true
end


"""
    isapprox(a::T, b::T; atol::Real=0.0, rtol::Real=sqrt(eps(Float64))) where {T<:TMType} -> Bool

where `TMType = Union{PrimaryTrack,PrimarySet,SatTrack,SecondarySet,MeasuredData,Intersection,ObservationSet,MeasuredSet,DataSet}`.

Fast nested approximate equality for CPro and CLay types.

In contrast to the Base `isapprox` method, this method provides NaN/missing-aware scalar semantics.
"""
function isapprox(a::T, b::T; atol::Real=0.0, rtol::Real=sqrt(eps(Float64)))::Bool where {T<:TMType}
    @inbounds for field in fieldnames(T)
        field == :metadata && continue  # ℹ skip metadata field for approximate equality check
        _nested_approx(getproperty(a, field), getproperty(b, field); atol=atol, rtol=rtol) || return false
    end
    return true
end


"""
    !=(a::T, b::T) where {T<:TMType} -> Bool

where `TMType = Union{PrimaryTrack,PrimarySet,SatTrack,SecondarySet,MeasuredData,Intersection,ObservationSet,MeasuredSet,DataSet}`.

Fast nested inequality for CPro and CLay types.
In contrast to the Base `!=` method, this method provides NaN/missing-aware scalar semantics and distinguishes between `0` and `-0` for floating-point types.
"""
!=(a::T, b::T) where {T<:TMType} = !(a == b)
