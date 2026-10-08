
@doc """raw
    rational_point_set(X::Scheme) -> RationalPointSet

For a scheme over ``k`` return its set of ``k``-rational points.
"""
@attr RationalPointSet{typeof(base_scheme(X)), S} function rational_point_set(X::S) where {S<:Scheme}
  return RationalPointSet(base_scheme(X), X)
end

domain(X::RationalPointSet) = X.domain
codomain(X::RationalPointSet) = X.codomain
coefficient_ring(X::RationalPointSet) = OO(domain(X))


function Base.show(io::IO,::MIME"text/plain", X::RationalPointSet)
  io = pretty(io)
  println(io, "Set of rational points")
  println(io, Indent(), "over " , Lowercase(), domain(X), Dedent())
  print(io, Indent(), "of " , Lowercase(), codomain(X))
end

function Base.show(io::IO, X::RationalPointSet)
  io = pretty(io)
  print(io, "Set of rational points ")
  print(io, "over " , Lowercase(), domain(X))
  print(io, " of " , Lowercase(), codomain(X))
end


function (S::RationalPointSet{<:Any,<:AbsAffineScheme})(c::Vector; check::Bool=true)
  k = coefficient_ring(S)
  c = k.(c)
  return AffineRationalPoint(S, c; check=check)
end

function (S::RationalPointSet{<:Any,<:AbsProjectiveScheme})(c::Vector; check::Bool=true)
  k = coefficient_ring(S)
  c = k.(c)
  return ProjectiveRationalPoint(S, c; check=check)
end

function (S::AbsRationalPointSet)(p::AbsRationalPoint)
  return S(coordinates(p))
end
