using Distributions: cdf, Gamma

abstract type AbstractUWError end

struct UWError{T<:Number, V<:AbstractVecOrMat{T}}<: AbstractUWError
data::V
nrep::AbstractVector{Int64}
end

function UWError(data::V) where {V<:AbstractVecOrMat{<:Number}}
    nrep = [length(data)]
    UWError(data,nrep)
end

function errorestimate(uwer::UWError; S::Float64=1.5, skip::Int =0 )
    data = uwer.data[skip+1:end]
    nrep = uwer.nrep .-skip
    N = length(data)
    if any(x->x<1,nrep) || sum(nrep) != N
        throw(ArgumentError("Inconsistent N and Nrep"))
    end
    R = length(nrep)
    mx = mean(data)
    mxr =  Iterators.map(mean,varpartition(data,nrep))
    fb = sum(mxr .* nrep)/N #weighted mean  of replica mmeans
delpro = data .- mx
Wmax = 0
Gint = 0.0
flag = false
if S!=0
    Wmax = fld(minimum(nrep),2)
    flag = true
end
GammaFbb = zeros(Float64, Wmax)
GammaFbb[1] = mean(delpro .^2)
     if GammaFbb[1] == 0 
    error(ErrorException("Error, no fluctuations!"))
     end
# findWopt!!!!(Wmax,Gint,flag,GammaFbb,delpro,nrep,N,R,S)
W::Int64 = 1
    while W<=Wmax
        GammaFbb[W+1] = 0
        i0=1
        for r in 1:R
            i1 = i0-1+nrep[r]
           GammaFbb[W+1] = GammaFbb[W+1] + sum( (delpro[i0:(i1-W)]) .*  delpro[(i0+W):i1]) 
            i0 = i0 +nrep[r]
        end
        GammaFbb[W+1] /= (N-R*W)
        if flag
            Gint +=GammaFbb[W+1]/GammaFbb[1]
            if Gint<0
                tauW = 5.e-16
            else
                tauW = S/log((Gint+1)/Gint)
            end
            gW = exp(-W/tauW)-tauW/sqrt(W*N)
            if gW<0
                Wopt = W
                Wmax = min(Wmax, 2*Wopt)
                flag =false
                
            end
        end
        W+=1
    end
if flag
    @info "Windowing condition failed"
    Wopt = Wmax
    
end
CFbbopt = GammaFbb[1] +2*sum(@view GammaFbb[2:(Wopt+1)])
if CFbbopt<=0
    @error "Gamma pathological: error^2 <0"
    return
end
@. GammaFbb +=CFbbopt/N
dGamma = gammaerror(GammaFbb, N , Wmax, 100)
CFbbopt = GammaFbb[1] +2*sum(GammaFbb[2:(Wopt+1)])
sigmaF = sqrt(CFbbopt/N)
 tauintFbb = cumsum(GammaFbb) ./ GammaFbb[1] .-0.5
if R>1
    bF = (fb-mx)/(R-1)
    mx = mx-bF
    if abs(bF)>sigmaF/4
        @info "a $(bF/sigmaF) bias of the mean has been cancelled"
        
    end
@.    mxr -= bF*N/nrep
@. fb -= bF*R
end
value = mx 
dvalue = sigmaF
ddvalue = dvalue*sqrt((Wopt+0.5)/N)
dtauintofW = tauintFbb[1:(Wmax+1)] .* sqrt.((0:Wmax) ./N) .*2
tauint =tauintFbb[Wopt+1]
dtauint = tauint*2*sqrt((Wopt-tauint+0.5)/N)
Qval=NaN
if R>1
    chisqr = sum((mxr .-Fb).^2*nrep)/CFbbopt
    Qval = 1- cdf(Gamma((R-1)/2,1),chisqr/2)
end
return GammaErrorReturnTye(value,dvalue,ddvalue,tauint,dtauint,Wopt,Wmax, tauintFbb[1:(Wmax+1)],dtauintofW[1:(Wmax+1)],Qval,S,N,R,nrep,data,GammaFbb,dGamma,true)
end

function findWopt!!!!(Wmax,Gint,flag,GammaFbb,delpro,nrep,N,R,S)
    W::Int64 = 1
    while W<=Wmax
        GammaFbb[W+1] = 0
        i0=1
        for r in 1:R
            i1 = i0-1+nrep[r]
           GammaFbb[W+1] = GammaFbb[W+1] + sum( (delpro[i0:(i1-W)]) .*  delpro[(i0+W):i1]) 
            i0 = i0 +nrep[r]
        end
        GammaFbb[W+1] /= (N-R*W)
        if flag
            Gint +=GammaFbb[W+1]/GammaFbb[1]
            if Gint<0
                tauW = 5.e-16
            else
                tauW = S/log((Gint+1)/Gint)
            end
            gW = exp(-W/tauW)-tauW/sqrt(W*N)
            if gW<0
                Wopt = W
                Wmax = min(Wmax, 2*Wopt)
                flag =false
                
            end
        end
        W+=1
    end
    
end
function gammaerror(gamma, N,W,Λ)
    err = zeros(Float64,W+1)
    gamma[(W+2):min((2*W+W+1),end)] .= 0
    for t in 0:W
        k = max(1,(t-W)):min((t+W), length(gamma)-W-1)
        err[t+1]  = sum(gamma[k.+t.+1]+(gamma[abs.(t.-k).+1]- 2 .* gamma[t+1] .*gamma[k.+1]).^2      )
        err[t+1] = sqrt(err[t+1]/N)
    end
    err[1] = 0.001
    return err

end