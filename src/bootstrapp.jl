using StatsBase: sample!
import Tables
using Statistics: mean, std
using Random: rand, seed!

# Type mapping
_promote_type(::Type{T}) where {T<:Integer} = FloatTypeForInt(T)
_promote_type(::Type{T}) where {T<:AbstractFloat} = T
_promote_type(::Type{Complex{U}}) where {U} = Complex{_promote_type(U)}

# helper: find float type matching integer width
FloatTypeForInt(::Type{Int8})  = Float16
FloatTypeForInt(::Type{UInt8}) = Float16
FloatTypeForInt(::Type{Int16}) = Float16
FloatTypeForInt(::Type{UInt16}) = Float16
FloatTypeForInt(::Type{Int32}) = Float32
FloatTypeForInt(::Type{UInt32}) = Float32
FloatTypeForInt(::Type{Int64}) = Float64
FloatTypeForInt(::Type{UInt64}) = Float64
FloatTypeForInt(::Type{Int128}) = Float64   # no Float128 in base Julia
FloatTypeForInt(::Type{UInt128}) = Float64  # ditto

abstract type AbstractBootstrap end

struct Bootstrap{ T<:Number, V<:AbstractVector{T}}<: AbstractBootstrap
data::V
seed ::Int
    function Bootstrap(_data::V,  _seed::Int) where {T<:Number,V<:AbstractVector{T} }
    seed!(_seed)
    new{T,V}(_data,_seed)
end
# ["Blocksize","μ", "σ", "δσ","τ_int", "Bias"]
end
function Bootstrap(data::V) where{T<:Number,V<:AbstractVector{T}}
     _seed = rand(Int)
    Bootstrap(data,_seed)
end
function Bootstrap(data::D, column::Symbol) where{D}
         Tables.istable(D) ||ArgumentError( "The Provided Data must implemment the Tables.jl interface")
     _seed = rand(Int)
    Bootstrap(Tables.getcolumn(data,column),_seed)
end

"""
    Bootstrap(f,data::D, column::Symbol) where{D} => Bootstrap

    Selects `column` form `data` and applies function `f` to it.
    `data` must implement the Tables.jl interface.


"""
function Bootstrap(f,data::D, column::Symbol) where{D}
         Tables.istable(D) ||ArgumentError( "The Provided Data must implemment the Tables.jl interface")
     _seed = rand(Int)
    Bootstrap(f(Tables.getcolumn(data,column)),_seed)
end

function ts_boot(bs::T,f;R::Int =500, skip::Int=0) where {T<:AbstractBootstrap} 
    data = bs.data[skip:length(data)]



end

"""
    boot(bs::T,f;R::Int =500, skip::Int=0) where {T<:AbstractBootstrap} => BootstrapResult

Generate R bootstrap samples from the data in `bs` and apply the function `f` to each sample.
    Returns a `BootstrapResult` containing the results of applying `f` to each bootstrap sample.
"""
function boot(bs::T,f;R::Int =500, skip::Int=0) where {T<:AbstractBootstrap}::BootstrapResult
    data = bs.data[1+skip:end]
     n = length(data)
    temp = Vector{Float64}(undef, n)
    result = Vector{Float64}(undef, R)
    for b in 1:R
        resample = sample!(data, temp, replace=true)
        result[b ] = f(resample)
    end

    res = BootstrapResult(result,f)
    return res
end