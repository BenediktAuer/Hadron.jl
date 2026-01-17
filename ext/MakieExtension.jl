module  MakieExtension

if isdefined(Base, :get_extension)
    using Makie
else
    using ..Makie
end

using Hadron
import Hadron: boothist,boothist!, BootstrapResult, TSBootstrapResult, GammaErrorReturnTye,AllObs
import Statistics: std
import Tables: columnnames, getcolumn, matrix


" QQNORM conversion"

Makie.convert_arguments(P::Type{<:QQNorm}, m::BootstrapResult; qqline=:none, col=1)=  Makie.convert_arguments(P, m[col]; qqline=qqline)
Makie.used_attributes(::Type{<:QQNorm},m::BootstrapResult ) = (:qqline,:col)

" Errorbars conversion"
Makie.convert_arguments(P::Type{<:Errorbars},x::AbstractVector, m::AbstractVector{<:BootstrapResult}; col=1) = Makie.convert_arguments(P,x,[getfield(k,:observable)[col] for k in m], [std(k[col]) for k in m ] )
Makie.convert_arguments(P::Type{<:Errorbars},x::AbstractVector, m::BootstrapResult, all::AllObs) = Makie.convert_arguments(P,x,getfield(m,:observable), vec(std(matrix(m),dims=1)) )
Makie.convert_arguments(P::Type{<:Errorbars}, x::Real, m::BootstrapResult; col=1) = Makie.convert_arguments(P,[x],[getfield(m,:observable)[col]], [std(m[col])])
Makie.used_attributes(::Type{<:Errorbars},x,m::Union{BootstrapResult, AbstractVector{<:BootstrapResult}}) = (:col,)


Makie.convert_arguments(P::Type{<:Errorbars}, m::Union{BootstrapResult, AbstractVector{<:BootstrapResult}},x::AbstractVector; col=1) = Makie.convert_arguments(P,[getfield(k,:observable)[col] for k in m],x, [std(k[col]) for k in m ]; direction = :x )
function Makie.convert_arguments(P::Type{<:Errorbars}, m::BootstrapResult, x::Real; col=1) 
    return Makie.convert_arguments(P,[getfield(m,:observable)[col]],[x], [std(m[col])])
end
Makie.used_attributes(::Type{<:Errorbars},m::Union{BootstrapResult, AbstractVector{<:BootstrapResult}},x) = (:col, )


" Scatter conversion"
Makie.convert_arguments(P::Type{<:Scatter}, x::AbstractVector, m::AbstractVector{<:BootstrapResult}; col=1  ) = Makie.convert_arguments(P, x,[getfield(k,:observable)[col] for k in m] )
Makie.convert_arguments(P::Type{<:Scatter}, x::Real, m::BootstrapResult; col=1) = Makie.convert_arguments(P,[x],[getfield(m,:observable)[col]])
Makie.used_attributes(::Type{<:Scatter},x,m::Union{BootstrapResult, AbstractVector{BootstrapResult}}) = (:col,)

Makie.convert_arguments(P::Type{<:Scatter}, m::AbstractVector{<:BootstrapResult}, x::AbstractVector; col=1  ) = Makie.convert_arguments(P, [getfield(k,:observable)[col] for k in m],x )
Makie.convert_arguments(P::Type{<:Scatter}, m::BootstrapResult, x::Real; col=1) = Makie.convert_arguments(P,[getfield(m,:observable)[col]],[x])
Makie.used_attributes(::Type{<:Scatter},m::Union{BootstrapResult, AbstractVector{<:BootstrapResult}},x) = (:col, )
Makie.convert_arguments(P::Type{<:Scatter}, m::BootstrapResult, all::AllObs) = Makie.convert_arguments(P,getfield(m,:observable))
Makie.convert_arguments(P::Type{<:Scatter}, x::AbstractVector,m::BootstrapResult, all::AllObs) = Makie.convert_arguments(P,x,getfield(m,:observable))

"Band Conversion"
Makie.convert_arguments(P::Type{<:Band},x::AbstractVector, m::BootstrapResult, all::AllObs) = Makie.convert_arguments(P,x,(getfield(m,:observable) .- vec(std(matrix(m),dims=1))), (getfield(m,:observable) .+vec(std(matrix(m),dims=1))) )
function Makie.convert_arguments(P::Type{<:Band}, m::BootstrapResult, all::AllObs) 
    obs = getfield(m,:observable)
    std_val = vec(std(matrix(m),dims=1))
    return Makie.convert_arguments(P,obs, obs .- std_val, obs .+ std_val)
end


"Lines Conversion"
Makie.convert_arguments(P::Type{<:Lines}, x::AbstractVector, m::BootstrapResult, all::AllObs) = Makie.convert_arguments(P,x,getfield(m,:observable))
Makie.convert_arguments(P::Type{<:Lines}, m::BootstrapResult, all::AllObs) = Makie.convert_arguments(P,getfield(m,:observable))
Makie.convert_arguments(P::Type{<:Lines}, m::AbstractVector{<:BootstrapResult}, x::AbstractVector; col=1  ) = Makie.convert_arguments(P, [getfield(k,:observable)[col] for k in m],x )
Makie.convert_arguments(P::Type{<:Lines}, x::AbstractVector, m::AbstractVector{<:BootstrapResult}; col=1  ) = Makie.convert_arguments(P,x, [getfield(k,:observable)[col] for k in m] )

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
    * a AbstractVector of `bins` colors
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
"""
    analyse(bootstrap::BootstrapResult)

Plots a histogram and Q-Q Plot of `bootstrap`.
Using qqnorm and boothist.

Attributes
==========

 `col = 1`:    If `bootstrap` has multiple observables, `col` selects which observables to plot. Type `Int`

 `QQmarkercolor=Makie.wong_colors()[2]`:  Color of Marker used in qqnorm

 `QQcolor = Makie.wong_colors()[1]`:  Color of the line used in qqnorm

 `kw_args...`:  Other Arguments which are forwarded to boothist, and qqnorm, for a list see those

"""
function Hadron.inspect(m::BootstrapResult, ;col=1, QQmarkercolor=Makie.wong_colors()[2], QQcolor = Makie.wong_colors()[1],kw_args...)

    # Label(fig[1, 1:2, Top()],  uppercasefirst(repr(getfield(m,:f)[col])[2:end]), valign = :top,
    # font = :bold,
    # padding = (5, 5, 30, 5))
    boothist(m,col=col,axis=(title = "Histogram",xlabel=repr(getfield(m,:f)[col]) ,ylabel="Hits"),kw_args... )
        current_figure()|> display
    qqnorm(m,col=col,markercolor=QQmarkercolor, color=QQcolor, qqline=:fit, axis=(title = "Q-Q Plot for $(repr(getfield(m,:f)[col]))",),kw_args...)
        current_figure()|> display
    return 
end

function Hadron.inspect(m::TSBootstrapResult)

    scatter( m[:Blocksize], m[:σ],axis=(xlabel="Blocklength", ylabel = L"\sigma"))
    errorbars!(m[:Blocksize], m[:σ], m[:δσ], whiskerwidth = 10)
    current_figure()|> display
    return 
end

function Hadron.inspect(m::GammaErrorReturnTye;pdf =false)

      hist(m.data,axis=( title = "Histogram",xlabel="observable" ,ylabel="Hits"))
      display(current_figure())
           if pdf
            save("Histogram.pdf",current_figure())
        end
      GammaFbb = m.Gamma/m.Gamma[1]
      err = Hadron.gammaerror(GammaFbb,m.N,m.Wₘₐₓ,100)
      Wopt = m.Wₒₚₜ
      Wmax = m.Wₘₐₓ
      vlines([Wopt],color=Makie.wong_colors()[2],axis=( ylabel=L"\Gamma(t)", xlabel=L"t"))
      hlines!([0],color=:black)
      scatter!( 0:Wmax,GammaFbb[1:Wmax+1])
      errorbars!( 0:Wmax,GammaFbb[1:Wmax+1], err[1:Wmax+1],whiskerwidth = 3)
        display(current_figure())
        if pdf
            save("Autocorrelationfunction.pdf",current_figure())
        end

        vlines([Wopt+1],color=Makie.wong_colors()[4],axis=(xlabel=L"W",ylabel=L"\tau_{int}(W)"))
        hlines!([m.τᵢₙₜW[Wopt]],color=Makie.wong_colors()[2])
      scatter!(m.τᵢₙₜW[1:Wmax])
    errorbars!(1:Wmax,m.τᵢₙₜW[1:Wmax],m.ΔτᵢₙₜW[1:Wmax],whiskerwidth = 3 )
display(current_figure())
     if pdf
            save("Window.pdf",current_figure())
        end
      return 
end

end