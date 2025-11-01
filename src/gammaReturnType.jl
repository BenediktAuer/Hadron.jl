using LaTeXStrings
abstract type AbstractGammaErrorReturnType end

struct GammaErrorReturnTye <: AbstractGammaErrorReturnType
observable
Δvalue
ΔΔvalue
τᵢₙₜ
Δτᵢₙₜ
Wₒₚₜ
Wₘₐₓ
τᵢₙₜW
ΔτᵢₙₜW
Q
S
N
R
nrep
data
Gamma
ΔGamma
primary
end

function Base.summary(uwerr::AbstractGammaErrorReturnType)
    println("Analysis bassed on $(uwerr.N) measurements" )
if uwerr.R>1
    println("split in $(uwerr.R) replica with $(uwerr.Nrep) measurements, respectively")
    
end
    if uwerr.primary == 1
    println(L"W_{opt}           = %$(uwerr.Wₒₚₜ)")
    println("value               = $(uwerr.observable)")
    println(L"\sigma(x)         = %$(uwerr.Δvalue)")
    println(L"\sigam(\sigma)(x) = %$(uwerr.ΔΔvalue)")
    println(L"\tau_{\text{int}} = %$(uwerr.τᵢₙₜ)")
    println(L"\sigma(\tau_{\text{int}}) = %$(uwerr.Δτᵢₙₜ)")
    if uwerr.R >1
        println("Q = $(uwerr.Q)")
    end
else
    println(L"W_{opt}           = %$(uwerr.Wₒₚₜ) for the %$(length(uwerr.Wₒₚₜ), observables)")

end

end
