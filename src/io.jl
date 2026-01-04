using PrettyTables: pretty_table
function Base.show(io::IO, m::AbstractBootstrapResult)
    pretty_table(m)

end

function Base.show(io::IO, m::AbstractBootstrap)
    println("$(typeof(m)) with $(size(m.data)) values")
    println("Seed: $(m.seed)")
end
"Struct to plot all Obs"
struct AllObs 
end