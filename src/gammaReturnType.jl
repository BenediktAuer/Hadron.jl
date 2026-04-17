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
    println("split in $(uwerr.R) replica with $(uwerr.nrep) measurements, respectively")
    
end
    if uwerr.primary == 1
    display(L"W_{opt}           = %$(uwerr.Wₒₚₜ)")
    display(L"\text{x}               = %$(uwerr.observable)")
    display(L"\sigma(x)         = %$(uwerr.Δvalue)")
    display(L"\sigma(\sigma(x)) = %$(uwerr.ΔΔvalue)")
    display(L"\tau_{\text{int}} = %$(uwerr.τᵢₙₜ)")
    display(L"\sigma(\tau_{\text{int}}) = %$(uwerr.Δτᵢₙₜ)")
    if uwerr.R >1
        display("Q = $(uwerr.Q)")
    end
else
    display(L"W_{opt}           = %$(uwerr.Wₒₚₜ) for the %$(length(uwerr.Wₒₚₜ), observables)")

end

end
