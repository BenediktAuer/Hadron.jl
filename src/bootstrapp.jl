using StatsBase: sample!
import Tables
using Statistics: mean, std
using Random: rand, seed!

include("helperfunctions.jl")
include("blocking.jl")

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

struct Bootstrap{ T<:Number, V<:AbstractVecOrMat{T}}<: AbstractBootstrap
data::V
seed ::Int
header::Vector{Symbol}
    function Bootstrap(_data::V,  _seed::Int, _header::Vector{Symbol}) where {T<:Number,V<:AbstractVecOrMat{T} }
    seed!(_seed)
    new{T,V}(_data,_seed,_header)
end
# ["Blocksize","μ", "σ", "δσ","τ_int", "Bias"]
end
function Bootstrap(data::V,) where{T<:Number,V<:AbstractVecOrMat{T}}
     _seed = rand(Int)
    Bootstrap(data,_seed,[Symbol("Column "*string(i)) for i in 1: size(data,2)])
end
function Bootstrap(data::D, column::Symbol) where{D}
         Tables.istable(D) ||ArgumentError( "The Provided Data must implemment the Tables.jl interface")
     _seed = rand(Int)
    Bootstrap(Tables.getcolumn(data,column),_seed, [column])
end

"""
    Bootstrap(f,data::D, column::Symbol) where{D} => Bootstrap

    Selects `column` form `data` and applies function `f` to it.
    `data` must implement the Tables.jl interface.


"""
function Bootstrap(f,data::D, column::Symbol) where{D}
         Tables.istable(D) ||ArgumentError( "The Provided Data must implemment the Tables.jl interface")
     _seed = rand(Int)
    Bootstrap(f.(Tables.getcolumn(data,column)),_seed,[column])
end

function Bootstrap(data::D; skip::Int=0) where{D}
         Tables.istable(D) ||ArgumentError( "The Provided Data must implemment the Tables.jl interface")
     _seed = rand(Int)
    Bootstrap(Matrix(data)[:,skip+1:end],_seed,Tables.columnnames(data)[skip+1:end])
end
function Bootstrap(f,data::D, column::Symbol, seed::Int) where{D}
         Tables.istable(D) ||ArgumentError( "The Provided Data must implemment the Tables.jl interface")
    Bootstrap(f.(Tables.getcolumn(data,column)),seed,[column])
end

function ts_boot(bs::T;l::Int = 2,R::Int =500, skip::Int=0) where {T<:AbstractBootstrap} 
    data = bs.data[1+skip:length(bs.data)]
    naive_σ = std(data)/sqrt(length(data))
    println("Naive: $(mean(data)) ± $(naive_σ) ")
    j=0
   
        temp_res = Vector{Float64}(undef, 6)
        Output = zeros(Float64,0,6)
    while length(data)/l>20
        j+=1
        blckdData = blocking(data,l)
        res = boot(Bootstrap(blckdData),[mean,std];R=R)
        temp_res[1] = l
        temp_res[2] = mean(blckdData)
        temp_res[3] = std(res[:mean])
        temp_res[4] = std(res[:std])/sqrt(length(blckdData))
        temp_res[5] = (std(res[:mean])^2/naive_σ^2)/2
        temp_res[6] = mean(blckdData)- mean(res[:mean])
        Output = vcat(Output, temp_res')
        if l<32
            l=l*2
        else
            l=l+20
            
        end

        
    end
    res = TSBootstrapResult(Output,[:Blocksize,:μ,:σ,:δσ,:τ_int,:Bias])


end

"""
    boot(bs::T,f;R::Int =500, skip::Int=0) where {T<:AbstractBootstrap} => BootstrapResult

Generate `R` bootstrap samples from the data in `bs` and apply the function `f` to each sample.
    Returns a `BootstrapResult` containing the results of applying `f` to each bootstrap sample.
    `f` can also be a Array of Functions
"""
function boot(bs::T,f;R::Int =500, skip::Int=0)::BootstrapResult where {T<:Bootstrap{<:Number,<:AbstractVector}}
    data = bs.data[1+skip:end]
     n = size(data)
    temp = Vector{Float64}(undef, n)
    result = Vector{Float64}(undef, R)
    for b in 1:R
     sample!(data, temp, replace=true)
        result[b ] = f(temp)
    end

    res = BootstrapResult(result,f,[f(data)],bs.header)
    return res
end
function boot(bs::T,f;R::Int =2, skip::Int=0)::BootstrapResult where {T<:Bootstrap{<:Number,<:AbstractArray{<:Number,2}}}
    data =  bs.data[ skip+1:last(axes(bs.data,1)), :]
     nrows,ncols = size(data)
    temp = similar(bs.data,Float64 )
    result = similar(bs.data,Float64, R, ncols)
    for b in 1:R
        direct_sample!(data, temp)
        res_row = f(temp,dims=1)

        for k in 1:ncols
        result[b,k ] = res_row[k]
    end
end
    res = BootstrapResult(result,f,vec(f(data,dims=1)),bs.header)
    return res
end
function boot(bs::T,f::F;R::Int =500, skip::Int=0)::BootstrapResult where {T<:Bootstrap{<:Number,<:AbstractVector}, F<: AbstractArray{Function}}
    data = bs.data[1+skip:end]
     n = length(data)
    temp = Vector{Float64}(undef, n)
    result = Matrix{Float64}(undef, R, length(f))
    for b in 1:R
        resample = sample!(data, temp, replace=true)
        for i in eachindex(f)
            result[b, i] = f[i](resample) 
            
        end
    end

    res = BootstrapResult(result,f, map(x->x(data), f),header=bs.header)
    return res
end
    
