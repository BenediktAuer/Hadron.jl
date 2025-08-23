module  MakieExtension

if isdefined(Base, :get_extension)
    using Makie
else
    using ..Makie
end

using Hadron
import Hadron: boothist,boothist!, BootstrapResult



Makie.convert_arguments(P::Type{<:QQNorm}, m::BootstrapResult; qqline=:none, col=1)= return Makie.convert_arguments(P, m[col]; qqline=qqline)
Makie.used_attributes(::Type{<:QQNorm},m::BootstrapResult ) = (:qqline,:col)

# @recipe Bootqq begin
#     qqline=:fit
#     Makie.mixin_generic_plot_attributes()...
# end
# function Makie.plot!(plt::QQNorm{<:Tuple{BootstrapResult}}) 
# input_nodes = [:converted_1, :col, :qqline]
# output_nodes = [:data, :qqline]
# map!(plt.attributes, input_nodes,output_nodes) do data, col, qqline

#     return (data[col],qqline)
# end
#     qqnorm!(plt,plt.attributes, plt.data; qqline=plt.qqline)
# end
# Makie.convert_singe_argument(m::BootstrapResult,i::Int=1) = m[i]
@recipe BootHist  begin 

        nbins = 5        # number of bins
        color = :steelblue    # histogram color
        linecolor = :red      # observable line color
        linewidth = 2
        col = 1

        Makie.mixin_generic_plot_attributes()...
end
# Makie.convert_arguments(P::Type{<:Analyse},v::BootstrapResult) = convert_arguments(P,v[1],getfield(v,:observable))
function Makie.plot!(plt::BootHist{<:Tuple{BootstrapResult}})
    input_nodes = [:converted_1, :col, :nbins] 
    output_nodes = [:counts, :bin_edges, :obs ]
    map!(plt.attributes, input_nodes, output_nodes) do data, col,nbins
        obs = getfield(data, :observable)[col]
        lo, hi = extrema(data[col])
 halfwidth = max(obs - lo, hi - obs)
 step = halfwidth / nbins
 edges = (obs - nbins*step) : step : (obs + nbins*step)
 return (data[col],edges, obs)
    end

        hist!(plt, plt.attributes, plt.counts; bins=plt.bin_edges)
        vlines!(plt, plt.attributes , plt.obs,color=plt.linecolor)



    return plt
end
function Hadron.analyse(m::BootstrapResult, ;col=1, QQmarkercolor=Makie.wong_colors()[2], QQcolor = Makie.wong_colors()[1],kw_args...)
    fig = Figure()
    a1 = Axis(fig[1,1], title = "Histogram",xlabel=String(getfield(m,:f)[col]) ,ylabel="Hits") 
    a2 = Axis(fig[1,2], title = "Q-Q Plot - $(String(getfield(m,:f)[col]))")
    boothist!(a1,m,col=col,kw_args... )
    qqnorm!(a2,m,col=col,markercolor=QQmarkercolor, color=QQcolor, qqline=:fit)
    return fig
end

end