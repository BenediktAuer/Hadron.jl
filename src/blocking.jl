
using Statistics: mean
  
function blocking(data,l)
    if l==1
        return data
    end
    if ndims(data) ==1

        return map( mean, Iterators.partition(data, l))
    end
    if ndims(data) != 2 
        error("Currently only 1D data can be blocked!")
    end
      nrows, ncols = size(data)
    N = (nrows ÷ l) * l                # largest multiple of l ≤ nrows
    nblocks = N ÷ l
    result = Array{Float64}(undef, nblocks, ncols)

    @inbounds for j in 1:nblocks
        row_start = (j-1)*l + 1
        row_end   = j*l
        for c in 1:ncols
            s = zero(Float64)
            @simd for r in row_start:row_end
                s += data[r,c]
            end
            result[j,c] = s / l
        end
    end

    return result

end