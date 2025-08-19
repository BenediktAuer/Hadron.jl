import  Tables
abstract type AbstractBootstrapResult<: Tables.AbstractColumns end
# To be a vaild AbstractBootstrapResult, the field data has to be defined
## make use of Tables interface
Tables.istable(::Type{<:AbstractBootstrapResult}) = true
names(m::AbstractBootstrapResult) = getfield(getfield(m, :data),:names)
matrix(m::AbstractBootstrapResult) = getfield(getfield(m, :data),:matrix)
lookup(m::AbstractBootstrapResult) = getfield(getfield(m, :data),:lookup)
Tables.schema(m::AbstractBootstrapResult) = Tables.schema((getfield(m, :data)))

#colum interface
Tables.columnaccess(::Type{<:AbstractBootstrapResult})=true
Tables.columns(m::AbstractBootstrapResult) = m
Tables.getcolumn(m::AbstractBootstrapResult, ::Type{T}, col::Int, nm::Symbol) where {T} = matrix(m)[:, col]
Tables.getcolumn(m::AbstractBootstrapResult, nm::Symbol) = matrix(m)[:, lookup(m)[nm]]
Tables.getcolumn(m::AbstractBootstrapResult, i::Int) = matrix(m)[:, i]
Tables.columnnames(m::AbstractBootstrapResult) = names(m)

# declare that any MatrixTable defines its own `Tables.rows` method
rowaccess(::Type{<:AbstractBootstrapResult}) = true
# just return itself, which means MatrixTable must iterate `Tables.AbstractRow`-compatible objects
rows(m::AbstractBootstrapResult) = rows(getfield(m, :data))
# the iteration interface, at a minimum, requires `eltype`, `length`, and `iterate`
# for `MatrixTable` `eltype`, we're going to provide a custom row type
Base.eltype(m::AbstractBootstrapResult)  = Base.eltype(getfield(m, :data))
Base.length(m::AbstractBootstrapResult) = length(getfield(m, :data))

Base.iterate(m::AbstractBootstrapResult, st=1) = iterate(getfield(m, :data), st)

getcolumn(m::eltype(AbstractBootstrapResult), ::Type, col::Int, nm::Symbol) =
    getfield(getfield(m, :source), :matrix)[getfield(m, :row), col]
getcolumn(m::eltype(AbstractBootstrapResult), i::Int) =
    getfield(getfield(m, :source), :matrix)[getfield(m, :row), i]
getcolumn(m::eltype(AbstractBootstrapResult), nm::Symbol) =
    getfield(getfield(m, :source), :matrix)[getfield(m, :row), getfield(getfield(m, :source), :lookup)[nm]]
columnnames(m::eltype(AbstractBootstrapResult)) = names(getfield(m, :source))

struct BootstrapResult{DF} <:AbstractBootstrapResult
    data::DF
    f::Symbol
    function BootstrapResult(data::AbstractVector{T}, func) where {T<:Number}
#use table default constructor to make Table
    df = Tables.table(data, header=[Symbol(func)])
    new{typeof(df)}(df, Symbol(func))
end
end


struct TSBootstrapResult{DF} <:AbstractBootstrapResult
    data::DF
    f::Symbol
function TSBootstrapResult(data::AbstractVector{T}, func) where {T<:Number}
#use table default constructor to make Table
    df = Tables.table( data,header=[Symbol(func)])
    new{typeof(df)}(df,Symbol(func))
end
end
