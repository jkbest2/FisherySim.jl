"""
Catchability describes a catchability surface for a vessel. Accepts either an
Array with the same dimensions as the DiscreteFisheryDomain or a function that
accepts a single argument (e.g. 2-element tuple or vector) of coordinates and
returns the catchability coefficient for that location. Note that at present
the function form is only used to form the array version under a particular
DiscreteFisheryDomain.
"""
struct Catchability{Tq}
    catchability::Tq

    Catchability(q::Tq) where Tq = new{Tq}(q)
end

## If q is constant over the domain, just pass the value. Don't need to pass the
## DiscreteFisheryDomain, but it's allowed to keep the interface consistent.
# Catchability(q::Tq) where Tq<:Real = Catchability{Tq}(q)
Catchability(q::Real, Ω::DiscreteFisheryDomain) = Catchability(q)

function Catchability(q::A, Ω::DiscreteFisheryDomain) where A <: AbstractArray
    size(q) == size(Ω) ||
        throw(DimensionMismatch("Catchability array and domain dimensions must match."))
    Catchability(q)
end

function Catchability(f::Function, Ω::DiscreteFisheryDomain)
    q = [f(loc) for loc in Ω.locs]
    Catchability(q, Ω)
end

"""
    Catchability(DD::DomainDistribution)

Construct a `Catchability` object as a random draw from a `DomainDistribution`.
"""
function Catchability(DD::DomainDistribution)
    Catchability(rand(DD))
end

"""
    Catchability(DD::DomainDistribution, C::Catchability)

Construct a `Catchability` using the spatially constant catchability from `C`
with multiplicative error based on a random draw from a `DomainDistribution`.
"""
function Catchability(DD::DomainDistribution, C::Catchability{<:Real})
    q = C.catchability .* rand(DD)
    Catchability(q)
end

function getindex(Q::Catchability{Tq} where Tq<:Real,
                  loc::NamedTuple{(:l, :t), Tuple{Ti, Ti}} where Ti<:Integer)
    Q.catchability
end

function getindex(Q::Catchability{Tq} where Tq<:AbstractMatrix,
                  loc::NamedTuple{(:l, :t), Tuple{Ti, Ti}} where Ti<:Integer)
    Q.catchability[loc.l]
end

function getindex(Q::Catchability{Tq} where Tq<:AbstractVector{<:AbstractMatrix},
                  loc::NamedTuple{(:l, :t), Tuple{Ti, Ti}} where Ti<:Integer)
    Q.catchability[loc.t][loc.l]
end

# function getindex(Q::Catchability{Tq}, i, j) where Tq<:Real
#     Q.catchability
# end
# function getindex(Q::Catchability{Tq}, i) where Tq<:Real
#     Q.catchability
# end
# function getindex(Q::Catchability{Tq}, i, j) where Tq<:AbstractArray{<:Any,2}
#     Q.catchability[i, j]
# end
# function getindex(Q::Catchability{Tq}, i) where Tq<:AbstractArray{<:Any,2}
#     Q.catchability[i]
# end
# function getindex(Q::Catchability{Tq}, i, j, t) where Tq<:AbstractVector{<:AbstractMatrix}
#     Q.catchability[t][i, j]
# end
# function getindex(Q::Catchability{Tq}, i, t) where Tq<:AbstractVector{<:AbstractMatrix}
#     Q.catchability[t][i]
# end
