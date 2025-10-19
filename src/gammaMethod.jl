abstract type AbstractUWError end

struct UWError{T::Number, V::AbstractVecOrMat{T}}<: AbstractUWError

data::V
nrep::AbstractVector{Int64}
end

function UWError(data::V) where {V<:AbstractVecOrMat{<:Number}}
    nrep = [length(data)]
    UWError(data,nrep)
end

function errorestimate(uwer::UWError; S::Float64=1.5 )
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
    Wmax = fld(min(nrep),2)
    flag = true
end
GammaFbb = zeros(Float64, Wmax)
GammaFbb[1] = mean(delpro .^2)
     if GammaFbb[1] == 0 
    error(ErrorException("Error, no fluctuations!"))
     end
findWopt!!!!(Wmax,Gint,flag,GammaFbb,delpro,nrep,N)
if flag
    @info "Windowing condition failed"
    Wopt = Wmax
    
end
CFbbopt = GammaFbb[1] +2*sum(@view GammaFbb[2:(Wopt+1)])
if CFbbopt<=0
    @error "Gamma pathological: error^2 <0"
    return
end
end

function findWopt!!!!(Wmax,Gint,flag,GammaFbb,delpro,nrep,N)
    W::Int64 = 1
    while W<=Wmax
        i0=1
        for r in 1:length(nrep)
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