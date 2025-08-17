using StatsBase: sample!
using DataFrames: AbstractDataFrame, DataFrame
using Statistics: mean, std
using Random: rand, seed!

abstract type AbstractBootstrap end
struct Bootstrap{ T<:Number, V<:AbstractVector{T},F, DF<:AbstractDataFrame}<: AbstractBootstrap
data::V
f::F
seed ::Int
R::Int
result::DF
    function Bootstrap(_data::V, _R::Int, _f::F, _seed::Int) where {T<:Number,V<:AbstractVector{T}, F }
    seed!(_seed)
    df = DataFrame(Symbol(_f) =>ones(T,_R ) )
    new{T,V,F,typeof(df)}(_data,_f,_seed,_R,df)
end
# ["Blocksize","μ", "σ", "δσ","τ_int", "Bias"]
end
function Bootstrap(data::V,R::Int, f::F) where{T<:Number,V<:AbstractVector{T},  F}
     _seed = rand(Int)
    Bootstrap(data,R,f,_seed)
end

function ts_boot(bs::T,skip::Int) where {T<:AbstractBootstrap} 
    data = bs.data[skip:length(data)]


end
function boot!(bs::T,skip::Int) where {T<:AbstractBootstrap} 
    data = bs.data[skip:end]
    df_col = Symbol(bs.f)
     n = length(data)
    temp = Vector{Float64}(undef, n)
    for b in 1:bs.R
        resample = sample!(data, temp, replace=true)
        bs.result[b,df_col ] =bs.f(resample)
    end
    return bs
end