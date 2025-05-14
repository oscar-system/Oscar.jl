@doc raw"""
    standard_walk(G::Oscar.IdealGens, start::MonomialOrdering, target::MonomialOrdering)

Compute a reduced Groebner basis w.r.t. to a monomial order by converting it using the Groebner Walk
using the algorithm proposed by [CKM97](@cite).

# Arguments
- `G::Oscar.IdealGens`: Groebner basis of an ideal with respect to a starting monomial order.
- `target::MonomialOrdering`: monomial order one wants to compute a Groebner basis for.
- `start::MonomialOrdering`: monomial order to begin the conversion.
"""
function standard_walk(
  G::Oscar.IdealGens,
  target::MonomialOrdering
)
  start = ordering(G)
  start_weight = canonical_matrix(start)[1, :]
  target_weight = canonical_matrix(target)[1, :]

  return standard_walk(G, target, start_weight, target_weight)
end

@doc raw"""
   standard_walk(G::Oscar.IdealGens, target::MonomialOrdering, current_weight::Vector{ZZRingElem}, target_weight::Vector{ZZRingElem})

Compute a reduced Groebner basis w.r.t. to the monomial order `target` by converting it using the Groebner Walk
using the algorithm proposed by [CKM97](@cite).

# Arguments
- `G::Oscar.IdealGens`: Groebner basis of an ideal with respect to a starting monomial order.
- `target::MonomialOrdering`: monomial order one wants to compute a Groebner basis for.
- `current_weight::Vector{ZZRingElem}`: weight vector representing the starting monomial order.
- `target_weight::Vector{ZZRingElem}`: weight vector representing the target weight.
"""
function standard_walk(
  ::Type{Oscar.IdealGens},
  G::Oscar.IdealGens,
  target::MonomialOrdering,
  current_weight::Vector{ZZRingElem},
  target_weight::Vector{ZZRingElem};
)
  @vprintln :groebner_walk "Results for standard_walk"
  @vprintln :groebner_walk "Crossed Cones in: "

  @v_do :groebner_walk steps = 0

  while current_weight != target_weight
    @vprintln :groebner_walk current_weight
    G = standard_step(G, current_weight, target)

    current_weight = next_weight(G, current_weight, target_weight)

    @v_do :groebner_walk steps += 1
    @vprintln :groebner_walk 2 G
    @vprintln :groebner_walk 2 "======="
  end

  @vprint :groebner_walk "Cones crossed: "
  @vprintln :groebner_walk steps

  return G
end

standard_walk(
  G::Oscar.IdealGens,
  target::MonomialOrdering,
  current_weight::Vector{ZZRingElem},
  target_weight::Vector{ZZRingElem};
) = gens(standard_walk(Oscar.IdealGens, G, target, current_weight, target_weight))

###############################################################
# The standard step is used for the strategy standard.
###############################################################

function standard_step(G::Oscar.IdealGens, w::Vector{ZZRingElem}, target::MonomialOrdering)
  current_ordering = ordering(G)
  next = weight_ordering(w, target)

  Gw = ideal(initial_forms(G, w))
  @vprint :groebner_walk 3 "GB of initial forms: "
  @vprintln :groebner_walk 3 Gw

  H = groebner_basis(Gw; ordering=next, complete_reduction=true)
  H = Oscar.groebner_lift(gens(H), next, gens(G), current_ordering)

  @vprint :groebner_walk 10 "Lifted GB of initial forms: "
  @vprintln :groebner_walk 10 H

  return Oscar.IdealGens(H, next)
end
standard_step(G::Oscar.IdealGens, w::Vector{Int}, T::Matrix{Int}) = standard_step(G, ZZ.(w), create_ordering(base_ring(G), w, T))

@doc raw"""
    initial_form(f::MPolyRingElem, w::Vector{ZZRingElem})

Return the initial form of `f` with respect to the weight vector `w`.
"""
function initial_form(f::MPolyRingElem, w::Vector{ZZRingElem})
  R = parent(f)

  ctx = MPolyBuildCtx(R)

  E = exponents(f)
  WE = dot.(Ref(w), Vector{ZZRingElem}.(E))
  maxw = -inf

  for (e, c, we) in zip(E, coefficients(f), WE)
    if we > maxw
      finish(ctx)
      maxw = we
    end
    if we == maxw
      push_term!(ctx, c, e)
    end
  end

  return finish(ctx)
end

@doc raw"""
    initial_forms(G::Oscar.IdealGens, w::Vector{ZZRingElem})

Return the initial form of each element in `G` with respect to the weight vector `w`.
"""
initial_forms(G::Oscar.IdealGens, w::Vector{ZZRingElem}) = initial_form.(G, Ref(w))
initial_forms(G::Oscar.IdealGens, w::Vector{Int}) = initial_form.(G, Ref(ZZ.(w)))

@doc raw"""
    next_weight(G::Oscar.IdealGens, current::Vector{ZZRingElem}, target::Vector{ZZRingElem})

Return the point furthest along the line segment conv(current,target) still in the starting cone 
as described in Algorithm 5.2 on pg. 437 of "Using algebraic geometry" (Cox, Little, O'Shea, 2005).

# Arguments
- G, a reduced GB w.r.t the current monomial order
- current, a weight vector in the current Gröbner cone (corresponding to G)
- target a target vector in the Gröbner cone of the target monomial order
"""
function next_weight(G::Oscar.IdealGens, current::Vector{ZZRingElem}, target::Vector{ZZRingElem})
  V = bounding_vectors(G)
  @vprint :groebner_walk 5 "Bounding vectors: "
  @vprintln :groebner_walk 5 V

  C = dot.(Ref(current), V)
  T = dot.(Ref(target), V)

  tmin = minimum(c//(c-t) for (c,t) in zip(C,T) if t<0; init=1)

  @vprint :groebner_walk 5 "Next rational weights: "
  @vprintln :groebner_walk 5 (QQ.(current) + tmin * QQ.(target-current))

  result = QQ.(current) + tmin * QQ.(target-current)
  if !iszero(result)
    return convert_bounding_vector(QQ.(current) + tmin * QQ.(target-current))
  else
    return zeros(ZZRingElem, length(result))
  end
end

@doc raw"""
    bounding_vectors(I::Oscar.IdealGens)

Return a list of "bounding vectors" of a Gröbner basis of `I`, as pairs of 
"exponent vector of leading monomial" and "exponent vector of tail monomial".
The bounding vectors form an H-description of the Gröbner cone. 
(cf. p. 437 [CLO05](@cite)) 
"""
function bounding_vectors(I::Oscar.IdealGens)
  # TODO: Marked Gröbner basis
  gens_by_terms = terms.(I; ordering=ordering(I))
  
  v = map(Iterators.peel.(gens_by_terms)) do (lead,tail)
      Ref(leading_exponent_vector(lead)) .- leading_exponent_vector.(tail)
  end

  return unique!(reduce(vcat, v))
end

