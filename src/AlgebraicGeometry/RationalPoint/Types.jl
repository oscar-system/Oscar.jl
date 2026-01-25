struct RationalPointSet{P<:AbsAffineScheme, T<:Scheme} <: AbsRationalPointSet{P,T}
  domain::P
  codomain::T
end

@doc raw"""
    FiniteRationalPointSet{P<:AbsAffineScheme, X<:Scheme, C<:Vector)

Iterator over a finite set of RationalPoints, each in terms of a vector of coordinates.
"""
struct FiniteRationalPointSet{P<:AbsAffineScheme, T<:Scheme, C<:Vector} <:AbsRationalPointSet{P,T}
  domain::P
  codomain::T
  coordinates_list::C

  function FiniteRationalPointSet(over::P,in::T,coord_list::C; check::Bool=true) where {P<:AbsAffineScheme,T<:Scheme,C<:Vector}
    r = new{P,T,C}(over, in, coord_list)
    @check begin
      F = coefficient_ring(ambient_coordinate_ring(in))
      all(all(Oscar.parent(x) == F for x in coord_list[i]) for i in 1:length(coord_list)) || error("coordinates do not lie in the base ring")
    end
    return r
  end
end

@doc raw"""
    AffineRationalPoint{CoeffType<:RingElem, ParentType<:RationalPointSet}

A rational point represented in terms of a vector of coordinates.

# Examples
```jldoctest
julia> A2 = affine_space(GF(2), [:x, :y]);

julia> (x, y) = coordinates(A2);

julia> X = algebraic_set(x*y);

julia> X([1, 0])
Rational point
  of V(x*y)
with coordinates (1, 0)

```
"""
struct AffineRationalPoint{S<:RingElem, T<:RationalPointSet} <: AbsAffineRationalPoint{S, T}
  # not mutable since this should be efficient
  coordinates::Vector{S}
  parent::T

  function AffineRationalPoint(parent::T, coordinates::Vector{S}; check::Bool=true) where {S <: RingElem, T<:AbsRationalPointSet}
    r = new{S, T}(coordinates, parent)
    @check begin
      X = codomain(parent)
      n = ngens(OO(ambient_space(X)))
      k = coefficient_ring(parent)

      all(Oscar.parent(i)===k for i in coordinates) || error("coordinates do not lie in the base ring")
      n == length(coordinates) || error("$(coordinates) must have length $(n)")
      _in(r, X)|| error("$(coordinates) does not lie in $(parent)")
    end
    return r
  end
end


@doc raw"""
    ProjectiveRationalPoint{CoeffType<:RingElem, ParentType<:AbsRationalPointSet}

Type for rational points in projective varieties.

# Examples
```jldoctest
julia> P2 = projective_space(QQ, 2);

julia> P2([4, 0 , 2//3])
Projective rational point
  of Projective 2-space over QQ with coordinates [s0, s1, s2]
with coordinates (4 : 0 : 2//3)

```
"""
struct ProjectiveRationalPoint{S<:RingElem, T<:AbsRationalPointSet} <: AbsProjectiveRationalPoint{S, T}
  # not mutable since this should be efficient
  coordinates::Vector{S}
  parent::T

  function ProjectiveRationalPoint(parent::T, coordinates::Vector{S}; check::Bool=true) where {S <: RingElem, T<:AbsRationalPointSet}
    r = new{S, T}(coordinates, parent)
    @check begin
      X = codomain(parent)
      n = relative_ambient_dimension(X) + 1
      k = coefficient_ring(parent)
      all(Oscar.parent(i)===k for i in coordinates) || error("coordinates do not lie in the base ring")
      n == length(coordinates) || error("$(coordinates) must have length $(n)")
      _is_projective(coordinates) || error("not a point in projective space")
      _in(r, X)|| error("$(coordinates) does not lie in $(parent)")
    end
    return r
  end
end
