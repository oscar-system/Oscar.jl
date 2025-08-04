
@doc raw"""
    groebner_walk(
      I::MPolyIdeal, 
      target::MonomialOrdering = lex(base_ring(I)),
      start::MonomialOrdering = default_ordering(base_ring(I));
      algorithm::Symbol = :standard
    )

Compute a reduced Groebner basis w.r.t. to the monomial ordering  `target` by converting it 
from a Groebner basis with respect to the ordering `start` using the Groebner Walk.

# Arguments
- `I::MPolyIdeal`: ideal one wants to compute a Groebner basis for.
- `target::MonomialOrdering=:lex`: monomial ordering one wants to compute a Groebner basis for.
- `start::MonomialOrdering=:degrevlex`: monomial ordering to begin the conversion.
- `algorithm::Symbol=:standard`: strategy of the Groebner Walk. One can choose between:
    - `standard`: Standard Walk [CLO05](@cite),
    - `generic`: Generic Walk [FJLT07](@cite),

# Examples

```jldoctest; setup=:(oldverb=get_verbosity_level(:groebner_walk);), teardown=:(set_verbosity_level(:groebner_walk, oldverb))
julia> R,(x,y) = polynomial_ring(QQ, [:x,:y]);

julia> I = ideal([y^4+ x^3-x^2+x,x^4]);

julia> groebner_walk(I, lex(R))
Gröbner basis with elements
  1: y^16
  2: x + y^12 - y^8 + y^4
with respect to the ordering
  lex([x, y])

julia> groebner_walk(I, lex(R); algorithm=:generic)
Gröbner basis with elements
  1: y^16
  2: x + y^12 - y^8 + y^4
with respect to the ordering
  lex([x, y])

julia> set_verbosity_level(:groebner_walk, 1);

julia> groebner_walk(I, lex(R))
Results for standard_walk
Crossed Cones in:
ZZRingElem[1, 1]
ZZRingElem[4, 3]
ZZRingElem[4, 1]
ZZRingElem[12, 1]
Cones crossed: 4
Gröbner basis with elements
  1: y^16
  2: x + y^12 - y^8 + y^4
with respect to the ordering
  lex([x, y])

```
"""
function groebner_walk(
  I::MPolyIdeal, 
  target::MonomialOrdering = lex(base_ring(I)),
  start::MonomialOrdering = default_ordering(base_ring(I));
  algorithm::Symbol = :standard
)
  if algorithm == :standard
    walk = (x) -> standard_walk(x, target)
  elseif algorithm == :generic
    walk = (x) -> generic_walk(x, start, target)
  else
    throw(NotImplementedError(:groebner_walk, algorithm))
  end

  Gb = groebner_basis(I; ordering=start, complete_reduction=true)
  Gb = walk(Gb)

  return Oscar.IdealGens(Gb, target; isGB=true)
end

@doc raw"""
    is_same_groebner_cone(G::Oscar.IdealGens, T::MonomialOrdering)

Check whether the leading terms of G with respect to the matrix ordering given by T agree 
with the leading terms of G with respect to the current ordering.
This means they are in the same cone of the Groebner fan. (cf. Lemma 2.2, Collart, Kalkbrener and Mall 1997)
"""
is_same_groebner_cone(G::Oscar.IdealGens, T::MonomialOrdering) = all(leading_term.(G; ordering=T) .== leading_term.(G; ordering=ordering(G)))

# converts a vector wtemp by dividing the entries with gcd(wtemp).
@doc raw"""
    convert_bounding_vector(w::Vector{T})

Scale a rational weight vector to have co-prime integer weights.
"""
convert_bounding_vector(w::Vector{T}) where {T<:Union{ZZRingElem, QQFieldElem}} = ZZ.(floor.(w//gcd(w)))

# Creates a weight ordering for R with weight cw and representing matrix T
create_ordering(R::MPolyRing, cw::Vector{L}, T::Matrix{Int}) where {L<:Number} = weight_ordering(cw, matrix_ordering(R, T))

# interreduces the Groebner basis G. 
# each element of G is replaced by its normal form w.r.t the other elements of G and the current monomial order 
# TODO reference, docstring
@doc raw"""
    autoreduce(G::Oscar.IdealGens)

Replace every element of G by the normal form with respect to the remaining elements of G and the current monomial ordering.
"""
function autoreduce(G::Oscar.IdealGens)
  generators = collect(gens(G))

  for i in 1:length(G)
    generators[i] = reduce(
      generators[i], generators[1:end .!= i]; ordering=G.ord, complete_reduction=true
    )
  end
  return Oscar.IdealGens(generators, G.ord; isGB=true)
end

#############################################
# unspecific help functions
#############################################
#TODO docstring

change_weight_vector(w::Vector{Int}, M::Matrix{Int}) = vcat(w', M[2:end, :])
change_weight_vector(w::Vector{ZZRingElem}, M::ZZMatrix) = vcat(w', M[2:end, :])

insert_weight_vector(w::Vector{Int}, M::Matrix{Int}) = vcat(w', M[1:end-1, :])
insert_weight_vector(w::Vector{ZZRingElem}, M::ZZMatrix) = vcat(w', M[1:end-1, :])

@doc raw"""
    add_weight_vector(w::Vector{ZZRingElem}, M::ZZMatrix)

Prepend the weight vector `w` as row to the matrix `M`.

"""
add_weight_vector(w::Vector{Int}, M::Matrix{Int}) = vcat(w', M)
add_weight_vector(w::Vector{ZZRingElem}, M::ZZMatrix) = ZZMatrix(vcat(w), M)

