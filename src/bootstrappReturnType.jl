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

# Tables.getcolumn(m::eltype(AbstractBootstrapResult), ::Type, col::Int, nm::Symbol) =
#     getfield(getfield(m, :source), :matrix)[getfield(m, :row), col]
# Tables.getcolumn(m::eltype(AbstractBootstrapResult), i::Int) =
#     getfield(getfield(m, :source), :matrix)[getfield(m, :row), i]
# Tables.getcolumn(m::eltype(AbstractBootstrapResult), nm::Symbol) =
#     getfield(getfield(m, :source), :matrix)[getfield(m, :row), getfield(getfield(m, :source), :lookup)[nm]]
# Tables.columnnames(m::eltype(AbstractBootstrapResult)) = names(m)

function describeBoot(m::AbstractBootstrapResult)
    println("Summary of Bootstrap Result: ")
    for idx in eachindex(Tables.columnnames(m))
        colname = Tables.columnnames(m)[idx]
        coldata = Tables.getcolumn(m, colname)
        println("Column: $colname")
        println("  Oberservable: $(getfield(m,:observable)[idx])")
        println("  std: $(std(coldata))")
        
    end
end

struct BootstrapResult{DF} <:AbstractBootstrapResult
    data::DF
    f::Vector{Symbol}
    observable::Vector{Float64}
    function BootstrapResult(data::AbstractVecOrMat{T}, func, observable::Float64) where {T<:Number}
#use table default constructor to make Table
    df = Tables.table(data, header=[Symbol(func)])
    new{typeof(df)}(df, [Symbol(func)],[observable])
end
function BootstrapResult(data::AbstractVecOrMat{T}, func::Vector{Function},observable::Vector{Float64}) where {T<:Number}
#use table default constructor to make Table
    func = map(f -> Symbol(f), func)
    df = Tables.table(data, header=func)
    new{typeof(df)}(df, func,observable)
end
end


struct TSBootstrapResult{DF} <:AbstractBootstrapResult
    data::DF
    f:: Vector{Symbol}


function TSBootstrapResult(data::AbstractVecOrMat{T}, names::Vector{Symbol}) where {T<:Number}
#use table default constructor to make Table

    df = Tables.table(data, header=names)
    new{typeof(df)}(df, [:μ])
end

end
