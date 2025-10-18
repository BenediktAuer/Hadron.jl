function skip_rows(A::AbstractArray{T,1}, skip) where {T}
    A[skip+1:last(axes(A,1))]
end

function skip_rows(A::AbstractArray{T,2}, skip) where {T}
    A[ skip+1:last(axes(A,1)), :]
end