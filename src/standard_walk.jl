###############################################################
# Implementation of the standard walk.
###############################################################
@doc raw"""
    standard_walk(G::Oscar.IdealGens, start::MonomialOrdering, target::MonomialOrdering)

Compute a reduced Groebner basis w.r.t. to a monomial order by converting it using the Groebner Walk
using the algorithm proposed by Collart, Kalkbrener & Mall (1997).

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
using the algorithm proposed by Collart, Kalkbrener & Mall (1997).

# Arguments
- `G::Oscar.IdealGens`: Groebner basis of an ideal with respect to a starting monomial order.
- `target::MonomialOrdering`: monomial order one wants to compute a Groebner basis for.
- `current_weight::Vector{ZZRingElem}`: weight vector representing the starting monomial order.
- `target_weight::Vector{ZZRingElem}`: weight vector representing the target weight.
"""
function standard_walk(
  G::Oscar.IdealGens,
  target::MonomialOrdering,
  current_weight::Vector{ZZRingElem},
  target_weight::Vector{ZZRingElem};
)
  @vprintln :groebner_walk "Results for standard_walk"
  @vprintln :groebner_walk "Crossed Cones in: "

  @v_do :groebner_walk steps = 0

  while current_weight != target_weight
    G = standard_step(G, current_weight, target)

    current_weight = next_weight(G, current_weight, target_weight)

    @v_do :groebner_walk steps += 1
    @vprintln :groebner_walk current_weight
    @vprintln :groebner_walk 2 G
  end

  @vprint :groebner_walk "Cones crossed: "
  @vprintln :groebner_walk steps

  return G
end

###############################################################
# The standard step is used for the strategies standard and perturbed.
###############################################################

@doc raw"""
    standard_step(G::Oscar.IdealGens, w::Vector{ZZRingElem}, target::MonomialOrdering)

TODO

# Arguments
- `G::Oscar.IdealGens`: Groebner basis of an ideal.
- `w::Vector{ZZRingElem}`: TODO
- `target::MonomialOrdering`: TODO
"""
function standard_step(G::Oscar.IdealGens, w::Vector{ZZRingElem}, target::MonomialOrdering)
  current_ordering = ordering(G)
  next = weight_ordering(w, target)

  Gw = ideal(initial_forms(G, w))

  H = groebner_basis(Gw; ordering=next, complete_reduction=true)
  H = lift(G, current_ordering, H, next)

  @vprintln :groebner_walk 3 Gw
  @vprintln :groebner_walk 3 H

  return interreduce_walk(H)
end
standard_step(G::Oscar.IdealGens, w::Vector{Int}, T::Matrix{Int}) = standard_step(G, ZZ.(w), create_ordering(base_ring(G), w, T))

@doc raw"""
    initial_form(f::MPolyRingElem, w::Vector{ZZRingElem})

Returns the initial form of `f` with respect to the weight vector `w`.
"""
function initial_form(f::MPolyRingElem, w::Vector{ZZRingElem})
  R = parent(f)

  ctx = MPolyBuildCtx(R)

  E = exponent_vectors(f)
  WE = dot.(Ref(w), Vector{ZZRingElem}.(E))
  maxw = maximum(WE)

  for (e, c, we) in zip(E, coefficients(f), WE)
    if we == maxw
      push_term!(ctx, c, e)
    end
  end

  return finish(ctx)
end

@doc raw"""
    initial_form(G::Oscar.IdealGens, w::Vector{ZZRingElem})

Returns the initial form of each element in `G` with respect to the weight vector `w`.
"""
initial_forms(G::Oscar.IdealGens, w::Vector{ZZRingElem}) = initial_form.(G, Ref(w))
initial_forms(G::Oscar.IdealGens, w::Vector{Int}) = initial_form.(G, Ref(ZZ.(w)))

@doc raw"""
    next_weight(G::Oscar.IdealGens, current::Vector{ZZRingElem}, target::Vector{ZZRingElem})

Returns the point furthest along the line segment conv(current,target) still in the starting cone 
as described in Algorithm 5.2 on pg. 437 of "Using algebraic geometry" (Cox, Little, O'Shea, 2005).

# Arguments
- G, a reduced GB w.r.t the current monomial order
- current, a weight vector in the current Gröbner cone (corresponding to G)
- target a target vector in the Gröbner cone of the target monomial order
"""
function next_weight(G::Oscar.IdealGens, current::Vector{ZZRingElem}, target::Vector{ZZRingElem})
  V = bounding_vectors(G)
  tmin = minimum(c//(c-t) for (c,t) in zip(dot.(Ref(current), V), dot.(Ref(target), V)) if t<0; init=1)

  @vprintln :groebner_walk 3 (QQ.(current) + tmin * QQ.(target-current))

  return QQ.(current) + tmin * QQ.(target-current) |> convert_bounding_vector
end

@doc raw"""
    bounding_vectors(I::Oscar.IdealGens)

Returns a list of "bounding vectors" of a Gröbner basis of `I`, as pairs of 
"exponent vector of leading monomial" and "exponent vector of tail monomial".
The bounding vectors form an H-description of the Gröbner cone. (cf. "Using algebraic geometry", pg. 437 (CLO, 2005)) TODO: consistent citations, compare with OSCAR
"""
function bounding_vectors(I::Oscar.IdealGens)
  # TODO: rename this to "BoundingVectors" or something similar (as in M2 implementation/master's thesis)
  # TODO: Marked Gröbner basis
  lead_exp = leading_term.(I; ordering=ordering(I)) .|> exponent_vectors .|> first
  # TODO: are leading terms being computed twice? (Once in leadexpv, once in tailexpvs) One instead could simply subtract leading terms, no? 
  tail_exps = zip(gens(I) .|> exponent_vectors, lead_exp) .|> splat((e, lead) -> filter(!=(lead), e))
  # tail_exps = tail.(I; ordering=ordering(I)) .|> exponent_vectors
  
  v = zip(lead_exp, tail_exps) .|> splat((l, t) -> Ref(l).-t)

  return unique!(reduce(vcat, v))
end