module  MakieExtension

if isdefined(Base, :get_extension)
    using Makie
else
    using ..Makie
end

using Hadron
import Hadron: histboot,histboot!, BootstrapResult


used_attributes(::Type{<:QQNorm},m::BootstrapResult ) = (:qqline,)
Makie.convert_arguments(P::Type{<:QQNorm},m::BootstrapResult; qqline = :none )= convert_arguments(P,m[1],qqline= qqline)
# Makie.convert_singe_argument(m::BootstrapResult,i::Int=1) = m[i]
@recipe Histboot ( result,) begin 

        nbins = 5        # number of bins
        color = :steelblue    # histogram color
        linecolor = :red      # observable line color
        linewidth = 2
        col = 1

        Makie.mixin_generic_plot_attributes()...
end
# Makie.convert_arguments(P::Type{<:Analyse},v::BootstrapResult) = convert_arguments(P,v[1],getfield(v,:observable))
function Makie.plot!(plt::Histboot)
    result = plt[:result][]
    nbins = plt[:nbins][]


    observables = getfield(result,:observable)
 function symmetric_bins(d, obs, nbins)
        lo, hi = extrema(d)
        # span symmetric around obs
        halfwidth = max(obs - lo, hi - obs)
        # step size so we get <= nbins bins per side
        step = halfwidth / nbins
        edges = (obs - nbins*step) : step : (obs + nbins*step)
        return collect(edges)
    end
    # Handle single vs multiple observables
        obs = observables[plt.col[]]
        d = result[plt.col[]]
 edges = symmetric_bins(d, obs, nbins)
#    grid =plt.plots[1,1] =GridLayout()
#     ax_hist = Axis(grid[1, 1])
#     ax_qq   = Axis(grid[1, 2])
        hist!(plt, plt.attributes, d; bins=edges)
        vlines!(plt, plt.attributes , [obs],color=plt.linecolor)



    return plt
end
function Hadron.analyse(m::BootstrapResult, ;col=1, QQmarkercolor=Makie.wong_colors()[2], QQcolor = Makie.wong_colors()[1],kw_args...)
    fig = Figure()
    a1 = Axis(fig[1,1], title = "Histogram",xlabel=String(getfield(m,:f)[col]) ,ylabel="Hits") 
    a2 = Axis(fig[1,2], title = "Q-Q Plot - $(String(getfield(m,:f)[col]))")
    histboot!(a1,m,col=col,kw_args... )
    qqnorm!(a2,m[col],markercolor=QQmarkercolor, color=QQcolor, qqline=:fit)
    return fig
end

end