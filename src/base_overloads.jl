# Import Base methods for overloading
import Base: ==, isapprox, isempty

TMType = Union{PrimaryTrack,PrimarySet,SatTrack,SecondarySet,MeasuredData,Intersection,ObservationSet,MeasuredSet,DataSet}

## Helper functions for fast nested equality/approximation with missing-aware scalar semantics

@inline _same_type(a::DataType, b::DataType) = a.name === b.name

@inline _nested_equal(a, b) = isequal(a, b)

@inline function _nested_equal(a::AbstractArray, b::AbstractArray)
    axes(a) == axes(b) || return false
    @inbounds for i in eachindex(a, b)
        _nested_equal(a[i], b[i]) || return false
    end
    return true
end

@inline function _nested_approx(a, b; atol::Real=0.0, rtol::Real=sqrt(eps(Float64)))
    if a isa TMType && b isa TMType
        return isapprox(a, b; atol, rtol)
    elseif a isa Number && b isa Number
        return isapprox(a, b; atol, rtol)
    elseif a isa AbstractString && b isa AbstractString
        return a == b
    else
        return isequal(a, b)
    end
end

@inline function _nested_approx(a::AbstractArray, b::AbstractArray; atol::Real=0.0, rtol::Real=sqrt(eps(Float64)))
    axes(a) == axes(b) || return false
    @inbounds for i in eachindex(a, b)
        _nested_approx(a[i], b[i]; atol, rtol) || return false
    end
    return true
end

@inline function _nested_approx(a::AbstractDataFrame, b::AbstractDataFrame; atol::Real=0.0, rtol::Real=sqrt(eps(Float64)))
    names(a) == names(b) || return false
    @inbounds for name in names(a)
        _nested_approx(a[!, name], b[!, name]; atol, rtol) || return false
    end
    return true
end


## Base overloads

"""
    ==(a::T1, b::T2) where {T1<:TMType,T2<:TMType} -> Bool

where `TMType = Union{PrimaryTrack,PrimarySet,SatTrack,SecondarySet,MeasuredData,Intersection,ObservationSet,MeasuredSet,DataSet}`.

Fast nested equality for _TrackMatcher_ types.

In contrast to the Base `==` method, this method provides NaN/missing-aware scalar semantics
and distinguishes between `0` and `-0` for floating-point types.
"""
function ==(a::T1, b::T2)::Bool where {T1<:TMType,T2<:TMType}
    _same_type(T1, T2) || return false
    @inbounds for field in fieldnames(T1)
        field == :metadata && continue  # ℹ skip metadata field for equality check
        _nested_equal(getproperty(a, field), getproperty(b, field)) || return false
    end
    return true
end


"""
    isapprox(a::T1, b::T2; atol::Real=0.0, rtol::Real=sqrt(eps(Float64))) where {T1<:TMType,T2<:TMType} -> Bool

where `TMType = Union{PrimaryTrack,PrimarySet,SatTrack,SecondarySet,MeasuredData,Intersection,ObservationSet,MeasuredSet,DataSet}`.

Fast nested approximate equality for _TrackMatcher_ types.

In contrast to the Base `isapprox` method, this method provides NaN/missing-aware scalar semantics.
"""
function isapprox(a::T1, b::T2; atol::Real=0.0, rtol::Real=sqrt(eps(Float64)))::Bool where {T1<:TMType,T2<:TMType}
    _same_type(T1, T2) || return false
    @inbounds for field in fieldnames(T1)
        field == :metadata && continue  # ℹ skip metadata field for approximate equality check
        _nested_approx(getproperty(a, field), getproperty(b, field); atol=atol, rtol=rtol) || return false
    end
    return true
end

function isempty(a::TMType)::Bool
    @inbounds for field in fieldnames(typeof(a))
        field == :metadata && continue  # ℹ skip metadata field for emptiness check
        isempty(getproperty(a, field)) || return false
    end
    return true
end
