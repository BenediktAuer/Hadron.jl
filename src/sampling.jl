using Random: default_rng,Sampler

"""
    direct_sample!([rng], a::AbstractArray, x::AbstractArray)

Direct sampling: for each `j` in `1:k`, randomly pick `i` from `1:n`,
and set `x[j] = a[i]`, with `n=length(a)` and `k=length(x)`.

This algorithm consumes `k` random numbers.
"""
function direct_sample!(rng::AbstractRNG, a::AbstractVector, x::AbstractVector)
    1 == firstindex(a) == firstindex(x) ||
        throw(ArgumentError("non 1-based arrays are not supported"))
    Base.mightalias(a, x) &&
        throw(ArgumentError("output array x must not share memory with input array a"))
    s = Sampler(rng, 1:length(a))
    @inbounds for i = eachindex(x)
        x[i] = a[rand(rng, s)]
    end
    return x
end
"""
    direct_sample!([rng], a::AbstractMatrix, x::AbstractMatrix)

Direct sampling: for each `j` in `1:k`, randomly pick `i` from `1:n`,
and set `x[j] = a[i]`, with `n=cols(a)` and `k=rows(x)`.

This algorithm consumes `k*n` random numbers.
"""
function direct_sample!(rng::AbstractRNG, a::AbstractMatrix, x::AbstractMatrix)
    nrows_a, ncols_a = size(a)
    nrows_x, ncols_x = size(x)

    ncols_a == ncols_x ||
        throw(ArgumentError("number of columns of a and x must match"))
    1 == firstindex(a) == firstindex(x) ||
        throw(ArgumentError("non 1-based arrays are not supported"))
    Base.mightalias(a, x) &&
        throw(ArgumentError("output array x must not share memory with input array a"))

    s = Sampler(rng, 1:nrows_a)
idxbuf = Vector{Int}(undef, nrows_x)
for j in 1:ncols_a
    rand!(rng, idxbuf, s)
    @inbounds @simd for k in 1:nrows_x
        x[k, j] = a[idxbuf[k], j]
    end
end
    return x
end
direct_sample!(a::T, x::T) where {T} = direct_sample!(default_rng(), a, x)