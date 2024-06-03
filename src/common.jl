
@doc raw"""
    groebner_walk(
      I::MPolyIdeal, 
      target::MonomialOrdering = lex(base_ring(I)),
      start::MonomialOrdering = default_ordering(base_ring(I));
      perturbation_degree = 2,
      algorithm::Symbol = :standard
    )

Compute a reduced Groebner basis w.r.t. to a monomial order by converting it using the Groebner Walk.
The Groebner Walk is proposed by Collart, Kalkbrener & Mall (1997).
One can choose a strategy of:
- Standard Walk (:standard) computes the Walk like as presented in Cox, Little & O'Shea (2005).
- Generic Walk (:generic) computes the Walk as presented in Fukuda, Jensen, Lauritzen & Thomas (2005).
- Perturbed Walk (:perturbed, with p = degree of the perturbation) computes the Walk as presented in Amrhein, Gloor & Küchlin (1997).

# Arguments
- `I::MPolyIdeal`: ideal one wants to compute a Groebner basis for.
- `target::MonomialOrdering=:lex`: monomial order one wants to compute a Groebner basis for.
- `start::MonomialOrdering=:degrevlex`: monomial order to begin the conversion.
- `perturbationDegree::Int=2`: the perturbation degree for the perturbed Walk.
- `algorithm::Symbol=standard`: strategy of the Groebner Walk. One can choose between:
    - `standard`: Standard Walk,
    - `generic`: Generic Walk,
    - `perturbed`: Perturbed Walk,

# Examples

```jldoctest
julia> R,(x,y) = polynomial_ring(QQ, ["x","y"]);

julia> I = ideal([y^4+ x^3-x^2+x,x^4]);

julia> groebner_walk(I, lex(R))
Gröbner basis with elements
1 -> x + y^12 - y^8 + y^4
2 -> y^16
with respect to the ordering
lex([x, y])

julia> groebner_walk(I, lex(R); algorithm=:perturbed)
Gröbner basis with elements
1 -> x + y^12 - y^8 + y^4
2 -> y^16
with respect to the ordering
lex([x, y])

julia> julia> set_verbosity_level(:groebner_walk, 1);
julia> groebner_walk(I, lex(R))
Results for standard_walk
Crossed Cones in: 
ZZRingElem[4, 3]
ZZRingElem[4, 1]
ZZRingElem[12, 1]
ZZRingElem[1, 0]
Cones crossed: 4
Gröbner basis with elements
1 -> x + y^12 - y^8 + y^4
2 -> y^16
with respect to the ordering
lex([x, y])

julia> groebner_walk(I, lex(R); algorithm=:perturbed)
perturbed_walk results
Crossed Cones in: 
[4, 3]
[4, 1]
[5, 1]
[12, 1]
[1, 0]
Cones crossed: 5
Gröbner basis with elements
1 -> y^16
2 -> x + y^12 - y^8 + y^4
with respect to the ordering
matrix_ordering([x, y], [1 0; 0 1])
```
"""
function groebner_walk(
  I::MPolyIdeal, 
  target::MonomialOrdering = lex(base_ring(I)),
  start::MonomialOrdering = default_ordering(base_ring(I));
  perturbation_degree = length(gens(base_ring(I))), # meaning, n=#gens(R)
  algorithm::Symbol = :standard
)
  if algorithm == :standard
    walk = (x) -> standard_walk(x, target)
  elseif algorithm == :generic
    walk = (x) -> generic_walk(x, start, target)
  elseif algorithm == :perturbed
    walk = (x) -> perturbed_walk(x, start, target, perturbation_degree)
  else
    throw(NotImplementedError(:groebner_walk, algorithm))
  end

  Gb = groebner_basis(I; ordering=start, complete_reduction=true)
  Gb = walk(Gb)

  return Oscar.IdealGens(Gb, target; isGB=true)
end

exponent_vectors = f->leading_exponent_vector.(monomials(f))

function weight_ordering(w::Vector{ZZRingElem}, o::MonomialOrdering)
    i = _support_indices(o.o)
    m = ZZMatrix(1, length(w), w)
    return MonomialOrdering(base_ring(o), MatrixOrdering(i, m, false))*o
end

# returns 'true' if the leading terms of G w.r.t the matrixordering T are the same as the leading terms of G with the current ordering.
same_cone(G::Oscar.IdealGens, T::MonomialOrdering) = all(leading_term.(G; ordering=T) .== leading_term.(G; ordering=ordering(G)))

# converts a vector wtemp by dividing the entries with gcd(wtemp).
convert_bounding_vector(w::Vector{T}) where {T<:Union{ZZRingElem, QQFieldElem}} = ZZ.(floor.(w//gcd(w)))

# TODO: This comment is factually not correct, it's just the weight ordering of a matrix ordering
# returns a copy of the PolynomialRing I, equipped with the ordering weight_ordering(cw)*matrix_ordering(T).
create_ordering(R::MPolyRing, cw::Vector{L}, T::Matrix{Int}) where {L<:Number} = weight_ordering(cw, matrix_ordering(R, T))

# interreduces the Groebner basis G. 
# each element of G is replaced by its normal form w.r.t the other elements of G and the current monomial order 
# TODO reference, docstring
# interreduces the Groebner basis G.
function autoreduce(G::Oscar.IdealGens)
  generators = collect(gens(G))

  for i in 1:length(gens(G))
    generators[i] = reduce(
      generators[i], generators[1:end .!= i]; ordering=G.ord, complete_reduction=true
    )
  end
  return Oscar.IdealGens(generators, G.ord; isGB=true)
end

#############################################
# unspecific help functions
#############################################

change_weight_vector(w::Vector{Int}, M::Matrix{Int}) = vcat(w', M[2:end, :])
change_weight_vector(w::Vector{ZZRingElem}, M::ZZMatrix) = vcat(w', M[2:end, :])

insert_weight_vector(w::Vector{Int}, M::Matrix{Int}) = vcat(w', M[1:end-1, :])
insert_weight_vector(w::Vector{ZZRingElem}, M::ZZMatrix) = vcat(w', M[1:end-1, :])

add_weight_vector(w::Vector{Int}, M::Matrix{Int}) = vcat(w', M)
add_weight_vector(w::Vector{ZZRingElem}, M::ZZMatrix) = ZZMatrix(vcat(w), M)

