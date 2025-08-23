module  MakieExtension

if isdefined(Base, :get_extension)
    using Makie
else
    using ..Makie
end

using Hadron
import Hadron: boothist,boothist!, BootstrapResult
import Statistics: std
" QQNORM conversion"

Makie.convert_arguments(P::Type{<:QQNorm}, m::BootstrapResult; qqline=:none, col=1)=  Makie.convert_arguments(P, m[col]; qqline=qqline)
Makie.used_attributes(::Type{<:QQNorm},m::BootstrapResult ) = (:qqline,:col)

" Errorbars conversion"
Makie.convert_arguments(P::Type{<:Errorbars},x::Vector, m::Vector{BootstrapResult}; col=1) = Makie.convert_arguments(P,x,[getfield(k,:observable)[col] for k in m], [std(k[col]) for k in m ] )
Makie.convert_arguments(P::Type{<:Errorbars}, x::Real, m::BootstrapResult; col=1) = Makie.convert_arguments(P,[x],[getfield(m,:observable)[col]], [std(m[col])])
Makie.used_attributes(::Type{<:Errorbars},x,m::Union{BootstrapResult, Vector{BootstrapResult}}) = (:col,)

Makie.convert_arguments(P::Type{<:Errorbars}, m::Union{BootstrapResult, Vector{BootstrapResult}},x::Vector; col=1) = Makie.convert_arguments(P,[getfield(k,:observable)[col] for k in m],x, [std(k[col]) for k in m ]; direction = :x )
function Makie.convert_arguments(P::Type{<:Errorbars}, m::BootstrapResult, x::Real; col=1) 
    return Makie.convert_arguments(P,[getfield(m,:observable)[col]],[x], [std(m[col])])
end
Makie.used_attributes(::Type{<:Errorbars},m::Union{BootstrapResult, Vector{BootstrapResult}},x) = (:col, )


" Scatter conversion"
Makie.convert_arguments(P::Type{<:Scatter}, x::Vector, m::Vector{BootstrapResult}; col=1  ) = Makie.convert_arguments(P, x,[getfield(k,:observable)[col] for k in m] )
Makie.convert_arguments(P::Type{<:Scatter}, x::Real, m::BootstrapResult; col=1) = Makie.convert_arguments(P,[x],[getfield(m,:observable)[col]])
Makie.used_attributes(::Type{<:Scatter},x,m::Union{BootstrapResult, Vector{BootstrapResult}}) = (:col,)

Makie.convert_arguments(P::Type{<:Scatter}, m::Vector{BootstrapResult}, x::Vector; col=1  ) = Makie.convert_arguments(P, [getfield(k,:observable)[col] for k in m],x )
Makie.convert_arguments(P::Type{<:Scatter}, m::BootstrapResult, x::Real; col=1) = Makie.convert_arguments(P,[getfield(m,:observable)[col]],[x])
Makie.used_attributes(::Type{<:Scatter},m::Union{BootstrapResult, Vector{BootstrapResult}},x) = (:col, )
"""
    boothist(bootstrap)

Plot a histogram of `bootsrap` centered around the observable idicated by a vertical line 

"""
@recipe BootHist  begin 
"""
Number of bins left and right of the value of the observable
"""
        nbins = 5        # number of bins

        """
    Color can either be:
    * a vector of `bins` colors
    * a single color
    * `:values`, to color the bars with the values from the histogram
        """
        color = :steelblue 
        """
        Color of the vertical line
        """
        linecolor = :red     
        linewidth = 2

        """
        If `bootstrap` has multiple observables, `col` selects which observables to plot.
        Type `Int`
        """
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
    Label(fig[1, 1:2, Top()], " ", valign = :top,
    font = :bold,
    padding = (5, 5, 20, 5))
    boothist!(a1,m,col=col,kw_args... )
    qqnorm!(a2,m,col=col,markercolor=QQmarkercolor, color=QQcolor, qqline=:fit)
    return fig
end

end