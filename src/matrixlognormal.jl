struct MatrixLogNormal{T<:Real,M<:AbstractMatrix,C<:AbstractPDMat} <:
    Distribution{Matrixvariate,Continuous}
    normal::MatrixNormal{T, M, C}

    function MatrixLogNormal(normal::MatrixNormal{T, M, C}) where {T, M, C}
        new{T, M, C}(normal)
    end
end

function MatrixLogNormal(M, U, V)
    MatrixLogNormal(MatrixNormal(M, U, V))
end

rand(MLN::MatrixLogNormal) = exp.(rand(MLN.normal))

function location(::Type{D}, s::Symbol,
                  M::AbstractMatrix,
                  U::AbstractMatrix, V::AbstractMatrix) where {D<:MatrixLogNormal}
    loc = similar(M)
    @inbounds for jdx in size(loc, 2)
        M[:, jdx] .= location(MvLogNormal, s,
                              M[:, jdx],
                              diagm(V[jdx, jdx] .* diag(U)))
    end
    loc
end

function location(::Type{D}, s::Symbol,
                  M::AbstractMatrix,
                  U::AbstractPDMat, V::AbstractPDMat) where {D<:MatrixLogNormal}
    location(D, s, M, U.mat, V.mat)
end

