struct MatrixLogNormal{D<:MatrixNormal}
    normal::D
end

function MatrixLogNormal(M::AbstracArray, U::AbstractArray, V::AbstractArray)
    matnorm = MatrixNormal(M, U, V)
    MatrixLogNormal(matnorm)
end

function location(::Type{MatrixLogNormal}, s::Symbol, m::AbstractMatrix, U::AbstractMatrix, V::AbstractMatrix)
    @assert size(m) == (size(V, 2), size(U, 1))
    _location!(D, Val{s}, m, U, V, similar(m))
end

function _location!(::Type{MatrixLogNormal}, ::Type{Val{:mean}}, mn::AbstractMatrix, U::AbstractMatrix, V::AbstractMatrix, μ::AbstractMatrix)
    II, JJ = size(mn)
    for jdx in 1:JJ
        for idx in 1:II
            μ[idx, jdx] = mn[idx, jdx] - V[idx, idx] * U[jdx, jdx] / 2
        end
    end
    μ
end

