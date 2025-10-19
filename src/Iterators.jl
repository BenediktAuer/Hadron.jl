# Implemntation of a custom partiion Iterators
using Base: HasShape, HasLength, IteratorSize

function varpartition(c,n::AbstractArray{<:Int,1})
     (any(x->x<1,n) || sum(n) != length(c)) && throw(ArgumentError("cannot create partitions of lengths $n"))
     sz = IteratorSize(typeof(c))
    sz isa Base.SizeUnknown && throw(ArgumentError("cannot create partitons with unknowen iterator"))
    sz isa Base.IsInfinite && throw(ArgumentError("cannot create partitons with infinite iterator"))
    return VarPartitionIterator(c,n)
end
struct VarPartitionIterator{T}
c::T
n::AbstractArray{<:Int,1}
end

eltype(::Type{VarPartitionIterator{T}}) where {T} = Vector{eltype(T)}
eltype(::Type{VarPartitionIterator{T}}) where {T<:AbstractArray} = AbstractVector{eltype(T)}
# But for some common implementations in Base we know the answer exactly
eltype(::Type{VarPartitionIterator{T}}) where {T<:Vector} = Base.SubArray{eltype(T), 1, T, Tuple{UnitRange{Int}}, true}
IteratorEltype(::Type{VarPartitionIterator{T}}) where {T} = IteratorEltype(T)
IteratorEltype(::Type{VarPartitionIterator{T}}) where {T<:AbstractArray} = EltypeUnknown()
IteratorEltype(::Type{VarPartitionIterator{T}}) where {T<:Vector} = IteratorEltype(T)
function Base.length(itr::VarPartitionIterator)
    return Base.length(itr.n)
end
varpartition_iteratorsize(::HasShape) = HasLength()
varpartition_iteratorsize(isz) = isz
function Base.IteratorSize(::Type{VarPartitionIterator{T}}) where {T} 
    varpartition_iteratorsize(IteratorSize(T))
end
function Base.iterate(itr::VarPartitionIterator{<:AbstractRange},state=(firstindex(itr.c),firstindex(itr.n)))
    state[1] >lastindex(itr.c) && return nothing
    r = min(state[1]+itr.n[state[2]]-1, lastindex(itr.c))
    return @inbounds itr.c[state[1]:r], (r+1,state[2]+1)
    
end
function Base.iterate(itr::VarPartitionIterator{<:AbstractArray},state=(firstindex(itr.c),firstindex(itr.n)))
    state[1] >lastindex(itr.c) && return nothing
    r = min(state[1]+itr.n[state[2]]-1, lastindex(itr.c))
    return @inbounds Base.view( itr.c,state[1]:r), (r+1,state[2]+1)
    
end