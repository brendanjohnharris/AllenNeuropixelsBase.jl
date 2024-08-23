using SparseArrays
using DimensionalData
export AbstractSparseToolsArray, SparseToolsArray
import DimensionalData: dims, refdims, data, name, metadata, parent, rebuild

abstract type AbstractSparseToolsArray{T, N, D, A} <: AbstractToolsArray{T, N, D, A} end

struct SparseToolsArray{T, N, D <: Tuple, R <: Tuple, A <: AbstractSparseArray{T, Int64, N},
                        Na, Me} <: AbstractSparseToolsArray{T, N, D, A}
    data::A
    dims::D
    refdims::R
    name::Na
    metadata::Me
end

function DimensionalData.format(dims, x::AbstractSparseToolsArray)
    DimensionalData.format(dims, axes(x))
end
# DimensionalData.format(::Base.KeySet, ::Type, ::AbstractRange) = DimensionalData.format(DimensionalData.Dimensions.LookupArrays.Aligned, ::Type, ::AbstractRange)
# format(val::AbstractSparseArray, D::Type, axis::AbstractRange) = format(AutoLookup(), D, val, axis)

function SparseToolsArray(data::A, dims::D, refdims::R = (),
                          name::Na = DimensionalData.NoName(),
                          metadata = DimensionalData.NoMetadata()) where {D, R, A, Na}
    @info "Constructing a sparsedimarray properly"
    SparseToolsArray(data, DimensionalData.format(dims, axes(data)), refdims, name,
                     metadata)
end

function goSparseToolsArray(data::A, dims::D, refdims::R = (),
                            name::Na = DimensionalData.NoName(),
                            metadata = DimensionalData.NoMetadata()) where {D, R, A, Na}
    @info "Constructing a sparsedimarray properly"
    SparseToolsArray(data, DimensionalData.format(dims, axes(data)), refdims, name,
                     metadata)
end

# function FeatureArray(data::AbstractArray, features::Union{Tuple{Symbol}, Vector{Symbol}}, args...)
#     FeatureArray(data, (Dim{:feature}(features), fill(AnonDim(), ndims(data)-1)...), args...)
# end

# FeatureArray(D::ToolsArray) = FeatureArray(D.data, D.dims, D.refdims, D.name, D.metadata)
# # DimensionalData.ToolsArray(D::FeatureArray) = ToolsArray(D.data, D.dims, D.refdims, D.name, D.metadata)

dims(A::AbstractSparseToolsArray) = A.dims
refdims(A::AbstractSparseToolsArray) = A.refdims
data(A::AbstractSparseToolsArray) = A.data
name(A::AbstractSparseToolsArray) = A.name
metadata(A::AbstractSparseToolsArray) = A.metadata
parent(A::AbstractSparseToolsArray) = data(A)
Base.Array(A::AbstractSparseToolsArray) = parent(A)

function DimensionalData.rebuild(A::AbstractSparseToolsArray, data::AbstractSparseArray,
                                 dims::Tuple, refdims::Tuple, name, metadata)
    SparseToolsArray(data, dims, refdims, name, metadata)
end

AbstractSparseDimVector = AbstractSparseToolsArray{T, 1, D, A} where {T, D, A}
