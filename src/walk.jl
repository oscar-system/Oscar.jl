
global infoLevel = 0

###############################################################
# Implementation of different variants of the Groebner Walk.
# The Groebner Walk is proposed by Collart, Kalkbrener & Mall (1997).
###############################################################
@doc raw"""
    groebner_walk(
        I::MPolyIdeal; 
        targetOrder::Symbol = lex(base_ring(I)), 
        startOrder::MonomialOrdering = default_ordering(base_ring(I));
        walktype::Symbol = :standard, 
        perturbationDegree::Int = 2
    )
    groebner_walk(
        G::Oscar.IdealGens,
        targetOrder::Union{Matrix{N}, MatElem{N}},
        startOrder::Union{Matrix{N}, MatElem{N}},
        walktype::Symbol = :standard,
        perturbationDegree::Int = 2
    ) where N

Compute a reduced Groebner basis w.r.t. to a monomial order by converting it using the Groebner Walk.
The Groebner Walk is proposed by Collart, Kalkbrener & Mall (1997).
One can choose a strategy of:
- Standard Walk (:standard) computes the Walk like as presented in Cox, Little & O´Shea (2005).
- Generic Walk (:generic) computes the Walk as presented in Fukuda, Jensen, Lauritzen & Thomas (2005).
- Perturbed Walk (:perturbed, with p = degree of the perturbation) computes the Walk ass presented in Amrhein, Gloor & Küchlin (1997).
- Tran's Walk (:tran) computes the Walk like as presented in Tran (2000).
- Fractal Walk (:fractalcombined) computes the Walk like as presented in Amrhein & Gloor (1998) with multiple extensions. The target monomial order has to be lex. This version uses the Buchberger algorithm to skip weight vectors with entries bigger than Int32.
- Fractal Walk (:fractal) computes the Walk as presented in Amrhein & Gloor (1998). Perturbs only the target vector.

# Arguments
- `I::MPolyIdeal`: ideal one wants to compute a Groebner basis for.
- `G::Oscar.IdealGens`: generators of an ideal one wants to compute a Groebner basis for.
- `target::MonomialOrdering=:lex`: monomial order one wants to compute a Groebner basis for.
- `start::MonomialOrdering=:degrevlex`: monomial order to begin the conversion.
- `walktype::Symbol=standard`: strategy of the Groebner Walk. One can choose a strategy of:
    - `standard`: Standard Walk,
    - `perturbed`: Perturbed Walk,
    - `tran`: Tran´s Walk,
    - `generic`: Generic Walk,
    - `fractal`: standard-version of the Fractal Walk,
    - `fractalcombined`: combined Version of the Fractal Walk. Target monomial order needs to be lex,
- `perturbationDegree::Int=2`: perturbationdegree for the perturbed Walk.

# Examples

```jldoctest
julia> R,(x,y) = polynomial_ring(QQ, ["x","y"]);

julia> I = ideal([y^4+ x^3-x^2+x,x^4]);

julia> groebner_walk(I, degrevlex(R), lex(R), :standard)
Gröbner basis with elements
1 -> y^16
2 -> x + y^12 - y^8 + y^4
with respect to the ordering
matrix_ordering([x, y], [1 0; 0 1])

julia> groebner_walk(I, degrevlex(R), lex(R), :perturbed, 2)
Gröbner basis with elements
1 -> y^16
2 -> x + y^12 - y^8 + y^4
with respect to the ordering
matrix_ordering([x, y], [1 0; 0 1])

julia> groebnerwalk(I, [1 1; 0 -1], [1 0; 0 1], :standard)
standard_walk results
Crossed Cones in: 
[4, 3]
[4, 1]
[12, 1]
[1, 0]
Cones crossed: 4
Gröbner basis with elements
1 -> y^16
2 -> x + y^12 - y^8 + y^4
with respect to the ordering
matrix_ordering([x, y], [1 0; 0 1])

julia> groebnerwalk(I, [1 1; 0 -1], [1 0; 0 1], :perturbed, 2)
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
  perturbation_degree = 2,
  algorithm::Symbol = :standard
)
  if algorithm == :standard
    walk = (x) -> standard_walk(x, start, target)
  elseif algorithm == :generic
    walk = (x) -> generic_walk(x, start, target)
  elseif algorithm == :perturbed
    walk = (x) -> perturbed_walk(x, start, target, perturbation_degree)
  else
    throw(NotImplementedError(:groebner_walk, algorithm))
  end

  Gb = groebner_basis(I; ordering=start, complete_reduction=true)
  Gb = walk(Gb)

  @vprintln :groebner_walk "Cones crossed: " 
  # delete_step_counter()

  return Oscar.IdealGens(gens(Gb), target; isGB=true)
end

function groebner_walk(
  G::Union{Oscar.IdealGens,MPolyIdeal},
  T::Union{Matrix{N},MatElem{N}},
  S::Union{Matrix{N},MatElem{N}},
  walktype::Symbol=:standard,
  p::Int=2,
) where {N}
  if walktype == :standard
    walk = (x) -> standard_walk(x, S, T)
  elseif walktype == :generic
    walk = (x) -> generic_walk(x, S, T)
  elseif walktype == :perturbed
    walk = (x) -> perturbed_walk(x, S, T, p)
  elseif walktype == :fractal
    walk = (x) -> fractal_walk(x, S, T)
  elseif walktype == :fractal_start_order
    walk = (x) -> fractal_walk_start_order(x, S, T)
    #    elseif walktype == :fractal_lex
    #       walk = (x) -> fractal_walk_lex(x, S, T)
    #   elseif walktype == :fractal_look_ahead
    #       walk = (x) -> fractal_walk_look_ahead(x, S, T)
  elseif walktype == :tran
    walk = (x) -> tran_walk(x, S, T)
  elseif walktype == :fractal_combined
    walk = (x) -> fractal_walk_combined(x, S, T)
  end
  delete_step_counter()
  S = Matrix{Int}(S)
  T = Matrix{Int}(T)
  # Make sure G is a fully reduced Groebner Basis
  Gb = groebner_basis(
    ideal(gens(G)); ordering=matrix_ordering(base_ring(G), S), complete_reduction=true
  )
  Gb = walk(Gb)

  @vprintln :groebner_walk "Cones crossed: " 
  delete_step_counter()
  return Oscar.IdealGens(gens(Gb), matrix_ordering(base_ring(G), T); isGB=true)
end

###########################################
# Counter for the steps
###########################################
counter = 0
function delete_step_counter()
  global counter
  temp = counter
  counter = 0
  return temp
end
function getstep_counter()
  global counter
  return counter
end
function raise_step_counter()
  global counter = getstep_counter() + 1
end

###############################################################
# Generic-version of the Groebner Walk.
###############################################################

function generic_walk(G::Oscar.IdealGens, start::MonomialOrdering, target::MonomialOrdering)
  Lm = leading_term.(G; ordering=start)
  V = next_gamma(G, Lm, zeros(ZZRingElem, 1), start, target)

  @vprintln :groebner_walk "Results for generic_walk"
  @vprintln :groebner_walk "Facets crossed for: "

  while !isempty(V)
    G, Lm = generic_step(G, Lm, V, target)
    
    # TODO: increase step_counter here
    @vprintln :groebner_walk V
    @vprintln :groebner_walk 2 G

    G = IdealGens(G, target)
    V = next_gamma(G, Lm, V, start, target)
  end

  return G
end
generic_walk(G::Oscar.IdealGens, S::Matrix{Int}, T::Matrix{Int}) = generic_walk(G, monomial_ordering(base_ring(G), S), monomial_ordering(base_ring(G), T))

function generic_step(
  G::Oscar.IdealGens, Lm::Vector{T}, v::Vector{ZZRingElem}, ord::MonomialOrdering
) where {T<:MPolyRingElem}
  facet_Generators = facet_initials(G, Lm, v)
  H = groebner_basis(
    ideal(facet_Generators); ordering=ord, complete_reduction=true, algorithm=:buchberger
  )
  H, Lm = lift_generic(gens(G), Lm, gens(H), ord)
  G = interreduce(H, Lm, ord)
  return G, Lm
end

###############################################################
# The Fractal Walk
###############################################################

##########################################
# global weightvectors
##########################################
pTargetWeights = []
pStartWeights = []
firstStepMode = false

###############################################################
# Combined version of the extensions of the Fractal Walk.
# This version
# - checks if the starting weight vector represents the monomial order and pertubes it if necessary.
# - analyses the Groebner basis Gw of the initialforms and uses the Buchberger-algorithm if the generators of Gw are binomial or less.
# - skips a step in top level in the last step.
# - checks if an entry of an intermediate weight vector is bigger than int32. In case of that the Buchberger-Algorithm is used to compute the Groebner basis of the ideal of the initialforms.
###############################################################
function fractal_walk_combined(G::Oscar.IdealGens, S::Matrix{Int}, T::Matrix{Int})
  @vprintln :groebner_walk "fractal_walk_combined results"
  @vprintln :groebner_walk "Crossed Cones in: "

  global pTargetWeights = [perturbed_vector(G, T, i) for i in 1:nvars(base_ring(G))]
  return fractal_walk_combined(G, S, T, S[1, :], pTargetWeights, 1)
end

function fractal_walk_combined(
  G::Oscar.IdealGens,
  S::Matrix{Int},
  T::Matrix{Int},
  currweight::Vector{Int},
  pTargetWeights::Vector{Vector{Int}},
  p::Int,
)
  R = base_ring(G)
  G.isGB = true
  w = currweight
  ordAlt = G.ord

  # Handling the weight of the start order.
  if (p == 1)
    if !ismonomial(initial_forms(G, w))
      global pStartWeights = [perturbed_vector(G, S, i) for i in 1:nvars(R)]
      global firstStepMode = true
    end
  end
  if firstStepMode
    w = pStartWeights[p]
  else
    w = currweight
  end

  # main loop
  while true
    t = next_weight_fractal_walk(G, w, pTargetWeights[p])

    # Handling the final step in the current depth.
    # next_weight_fractal_walk may return 0 if the target vector does not lie in the cone of T while G already defines the Groebner basis w.r.t. T.
    # -> Checking if G is already a Groebner basis w.r.t. T solves this problem and reduces computational effort since next_weight_fractal_walk returns 1 in the last step on every local path.
    if t == 1 && p != 1
      if same_cone(G, T)
        @vprintln :groebner_walk ("depth $p: in cone ", currweight, ".")...

        # Check if a target weight of pTargetWeights[p] and pTargetWeights[p-1] lies in the wrong cone.
        if !inCone(G, T, pTargetWeights, p)
          global pTargetWeights = [perturbed_vector(G, T, i) for i in 1:nvars(R)]
          @vprintln :groebner_walk ("depth $p: not in cone ", pTargetWeights[p], ".")...
        end
        return G
      end
    elseif t == [0] # The Groebner basis w.r.t. the target weight and T is already computed.
      if inCone(G, T, pTargetWeights, p)
        @vprintln :groebner_walk ("depth $p: in cone ", pTargetWeights[p], ".")...
        return G
      end
      global pTargetWeights = [perturbed_vector(G, T, i) for i in 1:nvars(R)]
      @vprintln :groebner_walk ("depth $p: not in cone ", pTargetWeights[p], ".")...
      continue
    end

    # skip a step for target monomial order lex.
    if t == 1 && p == 1
      @vprintln :groebner_walk ("depth $p: recursive call in ", pTargetWeights[p])...
      return fractal_walk_combined(G, S, T, w, pTargetWeights, p + 1)
    else
      w = w + t * (pTargetWeights[p] - w)
      w = convert_bounding_vector(w)
      Gw = ideal(initial_forms(G, w))

      # handling the current weight with regards to Int32-entries. If an entry of w is bigger than Int32 use the Buchberger-algorithm.
      if !checkInt32(w)
        println(
          "depth $p: Some entries of $w are bigger than Int32. Trying to find another weight,",
        )
        w, b = truncw(G, w, gens(Gw))
        if !b
          println("depth $p: Doing a direct conversion to the target monomial ordering.")
          ordNew = matrix_ordering(R, T)
          w = T[1, :]
          G = groebner_basis(
            ideal(G); ordering=ordNew, complete_reduction=true, algorithm=:buchberger
          )

          if !inCone(G, T, pTargetWeights, p)
            global pTargetWeights = [perturbed_vector(G, T, i) for i in 1:nvars(R)]
            @vprintln :groebner_walk ("depth $p: not in cone ", pTargetWeights[p], ".")...
          end
          return G
        end
      end
      ordNew = create_ordering(R, w, T)
      # converting the Groebner basis
      if (p == nvars(R) || isbinomial(gens(Gw)))
        H = groebner_basis(
          Gw; ordering=ordNew, complete_reduction=true, algorithm=:buchberger
        )

        @vprintln :groebner_walk ("depth $p: conversion in ", w, ".")...
        raise_step_counter()
      else
        @vprintln :groebner_walk "depth $p: recursive call in $w."
        H = fractal_walk_combined(
          Oscar.IdealGens(R, gens(Gw), ordAlt),
          S,
          T,
          deepcopy(currweight),
          pTargetWeights,
          p + 1,
        )
        global firstStepMode = false
      end
    end
    #H = liftGW2(G, ordAlt, Gw,H, ordNew)
    H = lift_fractal_walk(G, H, ordNew)
    G = interreduce_walk(H)
    ordAlt = ordNew
    currweight = w
  end
end

###############################################################
# Plain version of the Fractal Walk.
# This version checks if an entry of an intermediate weight vector is bigger than int32. In case of that the Buchberger-Algorithm is used to compute the Groebner basis of the ideal of the initialforms.
###############################################################

function fractal_walk(G::Oscar.IdealGens, S::Matrix{Int}, T::Matrix{Int})
  @vprintln :groebner_walk "FractalWalk_standard results"
  @vprintln :groebner_walk "Crossed Cones in: "

  global pTargetWeights = [perturbed_vector(G, T, i) for i in 1:nvars(base_ring(G))]
  return fractal_recursion_step(G, S, T, S[1, :], pTargetWeights, 1)
end

function fractal_recursion_step(
  G::Oscar.IdealGens,
  S::Matrix{Int},
  T::Matrix{Int},
  currweight::Vector{Int},
  pTargetWeights::Vector{Vector{Int}},
  p::Int,
)
  R = base_ring(G)
  G.isGB = true
  w = currweight
  ordAlt = G.ord
  while true
    t = next_weight_fractal_walk(G, w, pTargetWeights[p])

    # Handling the final step in the current depth.
    # next_weight_fractal_walk may return 0 if the target vector does not lie in the cone of T while G already defines the Groebner basis w.r.t. T.
    # -> Checking if G is already a Groebner basis w.r.t. T solves this problem and reduces computational effort since next_weight_fractal_walk returns 1 in the last step on every local path.        if t == 1 && p != 1
    if t == 1 && p != 1
      if same_cone(G, T)
        @vprintln :groebner_walk ("depth $p: in cone ", currweight, ".")...

        # Check if a target weight of pTargetWeights[p] and pTargetWeights[p-1] lies in the wrong cone.
        if !inCone(G, T, pTargetWeights, p)
          global pTargetWeights = [perturbed_vector(G, T, i) for i in 1:nvars(R)]
          @vprintln :groebner_walk ("depth $p: not in cone ", pTargetWeights[p], ".")...
        end
        return G
      end
    elseif t == [0] # The Groebner basis w.r.t. the target weight and T is already computed.
      if inCone(G, T, pTargetWeights, p)
        @vprintln :groebner_walk ("depth $p: in cone ", pTargetWeights[p], ".")...
        return G
      end
      global pTargetWeights = [perturbed_vector(G, T, i) for i in 1:nvars(R)]
      @vprintln :groebner_walk ("depth $p: not in cone ", pTargetWeights[p], ".")...
      continue
    end

    w = w + t * (pTargetWeights[p] - w)
    w = convert_bounding_vector(w)
    Gw = ideal(initial_forms(G, w))

    # Handling the current weight with regards to Int32-entries. If an entry of w is bigger than Int32 use the Buchberger-algorithm.
    if !checkInt32(w)
      println(
        "depth $p: Some entries of $w are bigger than Int32. Trying to find another weight,"
      )
      w, b = truncw(G, w, gens(Gw))
      if !b
        println("depth $p: Doing a direct conversion to the target monomial ordering.")
        ordNew = matrix_ordering(R, T)
        w = T[1, :]
        G = groebner_basis(
          Gw; ordering=ordNew, complete_reduction=true, algorithm=:buchberger
        )

        if !inCone(G, T, pTargetWeights, p)
          global pTargetWeights = [perturbed_vector(G, T, i) for i in 1:nvars(R)]
          @vprintln :groebner_walk ("depth $p: not in cone ", pTargetWeights[p], ".")...
        end
        return G
      end
    end
    ordNew = create_ordering(R, w, T)
    # Converting the Groebner basis
    if p == nvars(R)
      H = groebner_basis(
        Gw; ordering=ordNew, complete_reduction=true, algorithm=:buchberger
      )

      @vprintln :groebner_walk ("depth $p: conversion in ", w, ".")...
      raise_step_counter()
    else
      @vprintln :groebner_walk "depth $p: recursive call in $w."
      H = fractal_recursion_step(
        Oscar.IdealGens(R, gens(Gw), ordAlt),
        S,
        T,
        deepcopy(currweight),
        pTargetWeights,
        p + 1,
      )
    end
    #H = liftGW2(G, R, Gw, H, Rn)
    H = lift_fractal_walk(G, H, ordNew)
    G = interreduce_walk(H)
    ordAlt = ordNew
    currweight = w
  end
end

###############################################################
# Extends the plain Fractal Walk by checking the start order.
###############################################################

function fractal_walk_start_order(G::Oscar.IdealGens, S::Matrix{Int}, T::Matrix{Int})
  @vprintln :groebner_walk "fractal_walk_withStartorder results"
  @vprintln :groebner_walk "Crossed Cones in: "


  global pTargetWeights = [perturbed_vector(G, T, i) for i in 1:nvars(
    Singular.base_ring(G)
  )]
  return fractal_walk_recursive_startorder(G, S, T, S[1, :], pTargetWeights, 1)
end

function fractal_walk_recursive_startorder(
  G::Oscar.IdealGens,
  S::Matrix{Int},
  T::Matrix{Int},
  currweight::Vector{Int},
  pTargetWeights::Vector{Vector{Int}},
  p::Int,
)
  R = base_ring(G)
  G.isGB = true
  ordAlt = G.ord

  # Handling the starting weight.
  if (p == 1)
    if !ismonomial(initial_forms(G, currweight))
      global pStartWeights = [perturbed_vector(G, S, i) for i in 1:nvars(R)]
      global firstStepMode = true
    end
  end
  if firstStepMode
    w = pStartWeights[p]
  else
    w = currweight
  end

  while true
    t = next_weight_fractal_walk(G, w, pTargetWeights[p])

    # Handling the final step in the current depth.
    # next_weight_fractal_walk may return 0 if the target vector does not lie in the cone of T while G already defines the Groebner basis w.r.t. T.
    # -> Checking if G is already a Groebner basis w.r.t. T solves this problem and reduces computational effort since next_weight_fractal_walk returns 1 in the last step on every local path.        if t == 1 && p != 1
    if t == 1 && p != 1
      if same_cone(G, T)
        @vprintln :groebner_walk ("depth $p: in cone ", currweight, ".")...

        # Check if a target weight of pTargetWeights[p] and pTargetWeights[p-1] lies in the wrong cone.
        if !inCone(G, T, pTargetWeights, p)
          global pTargetWeights = [perturbed_vector(G, T, i) for i in 1:nvars(R)]
          @vprintln :groebner_walk ("depth $p: not in cone ", pTargetWeights[p], ".")...
        end
        return G
      end
    elseif t == [0] # The Groebner basis w.r.t. the target weight and T is already computed.
      if inCone(G, T, pTargetWeights, p)
        @vprintln :groebner_walk ("depth $p: in cone ", pTargetWeights[p], ".")...
        return G
      end
      global pTargetWeights = [perturbed_vector(G, T, i) for i in 1:nvars(R)]
      @vprintln :groebner_walk ("depth $p: not in cone ", pTargetWeights[p], ".")...
      continue
    end

    w = w + t * (pTargetWeights[p] - w)
    w = convert_bounding_vector(w)
    Gw = ideal(initial_forms(G, w))

    # Handling the current weight with regards to Int32-entries. If an entry of w is bigger than Int32 use the Buchberger-algorithm.
    if !checkInt32(w)
      w, b = truncw(G, w, gens(Gw))
      if !b
        ordNew = matrix_ordering(R, T)
        w = T[1, :]
        G = groebner_basis(
          Gw; ordering=ordNew, complete_reduction=true, algorithm=:buchberger
        )

        if !inCone(G, T, pTargetWeights, p)
          global pTargetWeights = [perturbed_vector(G, T, i) for i in 1:nvars(R)]
          @vprintln :groebner_walk ("depth $p: not in cone ", pTargetWeights[p], ".")...
        end
        return G
      end
    end
    ordNew = create_ordering(R, w, T)

    # Converting the Groebner basis
    if p == nvars(R)
      H = groebner_basis(
        Gw; ordering=ordNew, complete_reduction=true, algorithm=:buchberger
      )

      @vprintln :groebner_walk  ("depth $p: conversion in ", w, ".")...
      raise_step_counter()
    else
      @vprintln :groebner_walk "depth $p: recursive call in $w."
      H = fractal_walk_recursive_startorder(
        Oscar.IdealGens(R, gens(Gw), ordAlt),
        S,
        T,
        deepcopy(currweight),
        pTargetWeights,
        p + 1,
      )
      global firstStepMode = false
    end
    #H = liftGW2(G, R, Gw, H, Rn)
    H = lift_fractal_walk(G, H, ordNew)
    G = interreduce_walk(H)
    ordAlt = ordNew
    currweight = w
  end
end

###############################################################
# Tran´s version of the Groebner Walk.
# Returns the intermediate Groebner basis if an entry of an intermediate weight vector is bigger than int32.
###############################################################

function tran_walk(G::Oscar.IdealGens, S::Matrix{Int}, T::Matrix{Int})
  @vprintln :groebner_walk "tran_walk results"
  @vprintln :groebner_walk "Crossed Cones in: "

  currweight = S[1, :]
  tarweight = T[1, :]
  R = base_ring(G)
  if !ismonomial(initial_forms(G, currweight))
    currweight = perturbed_vector(G, S, nvars(R))
  end

  while true
    w = next_weight(G, currweight, tarweight)

    # return the Groebner basis if an entry of w is bigger than int32.
    if !checkInt32(w)
      w, b = truncw(G, w, initial_forms(G, w))
      !b && throw(
        error(
          "Some entries of the intermediate weight-vector $w are bigger than int32. Choose an other algorithm instead.",
        ),
      )
    end
    if w == tarweight
      if same_cone(G, T)
        @vprintln :groebner_walk ("Cones crossed: ", counter)...
        return G
      elseif inSeveralCones(initial_forms(G, tarweight))
        tarweight = representation_vector(G, T)
        continue
      end
    end
    G = standard_step_without_int32_check(G, w, T)

    @vprintln :groebner_walk w
    @vprintln :groebner_walk 2 G

    currweight = w
    raise_step_counter()
  end
end

###############################################################
# Standard step without checking of the entries of a given weight vector.
###############################################################

function standard_step_without_int32_check(
  G::Oscar.IdealGens, w::Vector{Int}, T::Matrix{Int}
)
  R = base_ring(G)
  ordAlt = G.ord
  ordNew = create_ordering(R, w, T)
  Gw = initial_forms(G, w)
  H = groebner_basis(
    ideal(Gw); ordering=ordNew, complete_reduction=true, algorithm=:buchberger
  )
  #H = liftGW2(G, R, Gw, H, Rn)
  H = lift(G, ordAlt, H, ordNew)
  return interreduce_walk(H)
end

#################################################################
# Procedures of the fractal walk.
# The fractal walk is proposed by Amrhein & Gloor (1998).
#################################################################

# lifts the Groebner basis G to the Groebner basis w.r.t. in the Fractal Walk like it´s done in Fukuda et. al (2005).
function lift_fractal_walk(
  G::Oscar.IdealGens, H::Oscar.IdealGens, orderingNew::MonomialOrdering
)
  R = base_ring(G)
  G.isGB = true
  G = Oscar.IdealGens(
    R,
    [
      gen - Oscar.IdealGens(
        [reduce(gen, gens(G); ordering=G.ord, complete_reduction=true)], H.ord
      )[1] for gen in gens(H)
    ],
    orderingNew,
  )
  G.isGB = true
  return G
end

# returns ´true´ if all polynomials of the given array are monomials.
function ismonomial(Gw::Vector{T}) where {T<:MPolyRingElem}
  for g in Gw
    if length(coefficients(g)) > 1
      return false
    end
  end
  return true
end

# returns ´true´ if all polynomials of the given array are binomials or less.
function isbinomial(Gw::Vector{T}) where {T<:MPolyRingElem}
  for g in Gw
    if length(coefficients(g)) > 2
      return false
    end
  end
  return true
end

# returns the next t to compute the next weight vector w(t) = w + t * (tw - w) like it´s done in Amrhein & Gloor (1998). This Method is NOT tested sufficiently.
function nextT(
  G::Oscar.IdealGens, w::Array{T,1}, tw::Array{K,1}
) where {T<:Number,K<:Number}
  if (w == tw)
    return [0]
  end
  tmin = 2
  t = 0
  for g in gens(G)
    a = Singular.leading_exponent_vector(g)
    d = Singular.exponent_vectors(tail(g))
    for v in d
      frac = (dot(w, a) - dot(w, v) + dot(tw, v) - dot(tw, a))
      if frac != 0
        t = (dot(w, a) - dot(w, v))//frac
      end
      if t > 0 && t < tmin
        tmin = t
      end
    end
  end
  if tmin <= 1
    return tmin
  else
    return [0]
  end
end

# returns the next t to compute the next weight vector w(t) = w + t * (tw - w) like it´s done in Cox, Little & O'Sheao (2005).
function next_weight_fractal_walk(
  G::Oscar.IdealGens, cweight::Array{T,1}, tweight::Array{K,1}
) where {T<:Number,K<:Number}
  if (cweight == tweight)
    return [0]
  end
  tmin = 1
  for v in difference_lead_tail(G)
    tw = dot(tweight, v)
    if tw < 0
      cw = dot(cweight, v)
      t = cw//(cw - tw)
      if t < tmin
        tmin = t
      end
    end
  end

  # BigInt is needed to prevent overflows in the conversion of the weight vectors.
  return BigInt(numerator(tmin))//BigInt(denominator(tmin))
end

# returns 'true' if the leading terms of G w.r.t the matrixordering T are the same as the leading terms of G w.r.t the weighted monomial ordering with weight vector of pvecs[p] (pvecs[p-1]) and the matrixordering T.
function inCone(G::Oscar.IdealGens, T::Matrix{Int}, pvecs::Vector{Vector{Int}}, p::Int)
  # pvecs(1) equals the first row of T
  if p == 1
    return true
  end
  ord = matrix_ordering(base_ring(G), T)
  cvzip = zip(gens(G), initial_forms(G, pvecs[p - 1]), initial_forms(G, pvecs[p]))
  for (g, in, in2) in cvzip
    if !isequal(leading_term(g; ordering=ord), leading_term(in; ordering=ord)) ||
      !isequal(leading_term(g; ordering=ord), leading_term(in2; ordering=ord))
      return false
    end
  end
  return true
end

#################################################################
# Procedures of the generic walk.
# The generic walk is proposed by Fukuda, Lauritzen & Thomas (2005).
#################################################################

# returns the initials of the polynomials w.r.t. the vector v.
function facet_initials(
  G::Oscar.IdealGens, lm::Vector{<:MPolyRingElem}, v::Vector{ZZRingElem}
)
  initials = Vector{MPolyRingElem}()

  # TODO: this really wants marked Gröbner bases
  ctx = MPolyBuildCtx(base_ring(G))
  for (g,el) in zip(G, exponent_vectors.(lm) .|> first)
    for (c,e) in zip(coefficients(g), exponent_vectors(g))
      if el == e || is_parallel(ZZ.(el - e), v)
        push_term!(ctx, c, e)
      end
    end

    push!(initials, finish(ctx))
  end

  return initials
end

# returns the differences of the exponent vectors of the leading terms and the polynomials of the generators of I.
#Comment: I shouldn't need the ordering T here. (Lm already consists of monomials)

function difference_lead_tail(
  G::Oscar.IdealGens, Lm::Vector{L}, T::MonomialOrdering
) where {L<:MPolyRingElem}
  lead_exp = leading_term.(Lm; ordering=T) .|> exponent_vectors .|> first
  
  v = zip(lead_exp, exponent_vectors.(G)) .|> splat((l, t) -> Ref(l).-t)

  return unique!(reduce(vcat, v))
end

difference_lead_tail(::Type{ZZRingElem}, G::Oscar.IdealGens, Lm::Vector{<:MPolyRingElem}, T::MonomialOrdering) = Vector{ZZRingElem}.(difference_lead_tail(G, Lm, T))

# function difference_lead_tail(
#   G::Oscar.IdealGens, Lm::Vector{L}, T::Union{Matrix{N},MatElem{N}}
# ) where {L<:MPolyRingElem,N}
#   v = Vector{Int}[]
#   for i in 1:length(G)
#     ltu = collect(
#       AbstractAlgebra.exponent_vectors(
#         leading_term(Lm[i]; ordering=matrix_ordering(base_ring(G), T))
#       ),
#     )[1]
#     for e in AbstractAlgebra.exponent_vectors(G[i])
#       if ltu != e
#         push!(v, ltu .- e)
#       end
#     end
#   end
#   return unique!(v)
# end

# returns true if the vector u is parallel to the vector v.
function is_parallel(u::Vector{T}, v::Vector{T}) where T
  count = 1
  x = 0
  for i in 1:length(u)
    if u[i] == 0
      if v[count] == 0
        count += +1
      else
        return false
      end
    else
      x = v[count]//u[i]
      count += 1
      break
    end
  end
  if count > length(v)
    return true
  end
  for i in count:length(v)
    @inbounds if v[i] != x * u[i]
      return false
    end
  end
  return true
end

# performs the lifting in the generic Walk like it´s proposed by Fukuda et al. (2005).
function lift_generic(
  G::Vector{T}, Lm::Vector{T}, H::Vector{T}, ord::MonomialOrdering
) where {T<:MPolyRingElem}
  R = parent(first(G))
  Newlm = Array{elem_type(R),1}(undef, 0)
  liftPolys = Array{elem_type(R),1}(undef, 0)
  for g in H
    push!(Newlm, leading_term(g; ordering=ord))
    push!(liftPolys, g - reduce_walk(g, G, Lm, ord))
  end
  return liftPolys, Newlm
end

# returns all v \in V if v<0 w.r.t. the ordering represented by T and v>0 w.r.t the ordering represented by S.
function filter_by_ordering(S::ZZMatrix, T::ZZMatrix, V::Vector{Vector{ZZRingElem}})
  pred = v->(
    less_than_zero(T, v) && 
    greater_than_zero(S, v)
  )
  return unique!(filter(pred, V))
end

# function filter_by_ordering(S::Matrix{Int}, T::Matrix{Int}, V::Vector{Vector{Int}})
#   btz = Set{Vector{Int}}()
#   for v in V
#     if less_than_zero(T, v) && greater_than_zero(S, v)
#       push!(btz, v)
#     end
#   end
#   return btz
# end

# returns all v \in V if w<v w.r.t. the facet-preorder.
# function filter_lf(w::Vector{Int}, S::Matrix{Int}, T::Matrix{Int}, V::Set{Vector{Int}})
#   btz = Set{Vector{Int}}()
#   for v in V
#     if less_facet(w, v, S, T)
#       push!(btz, v)
#     end
#   end
#   return btz
# end

function filter_lf(w::Vector{ZZRingElem}, start::ZZMatrix, target::ZZMatrix, V::Vector{Vector{ZZRingElem}})
  pred = v->less_facet(w, v, start, target)

  return unique!(filter(pred, V))
end

# computes the next vector in the generic walk.
function next_gamma(
  G::Oscar.IdealGens, 
  Lm::Vector{<:MPolyRingElem}, 
  w::Vector{ZZRingElem}, 
  start::MonomialOrdering, 
  target::MonomialOrdering
)
  S = canonical_matrix(start)
  T = canonical_matrix(target)
  
  V = filter_by_ordering(S, T, difference_lead_tail(ZZRingElem, G, Lm, target))

  if (w != zeros(ZZRingElem, 1))
    V = filter_lf(w, S, T, V) #TODO
  end

  if isempty(V)
    return V
  end

  op = (minV,v) ->
    if less_facet(v, minV, S, T)
      return v
    else return minV
  end

  return foldl(op, V; init=first(V))
end

# tests if v>0 w.r.t. the ordering M.
# TODO What is the definition?
# This should be the ordering on Q^n induced by M: ( u <_M v iff Mu <_lex Mv)?  
function greater_than_zero(M::ZZMatrix, v::Vector{ZZRingElem})
  nrows, _ = size(M)

  for i in 1:nrows
    d = dot(M[i,:], v)

    if d != 0
      return d > 0
    end
  end

  return false
end

# tests if v<0 w.r.t. the ordering M.
# less_than_zero(M::ZZMatrix, v::Vector{ZZRingElem}) = new_less_than_zero(M, v)
function less_than_zero(M::ZZMatrix, v::Vector{ZZRingElem})
  nrows, _ = size(M)

  for i in 1:nrows
    d = dot(M[i,:], v)

    if d != 0
      return d < 0
    end
  end

  return false
end

# tests if u<v w.r.t. the facet-preorder represented by the matrices S and T.
function less_facet(u::Vector{ZZRingElem}, v::Vector{ZZRingElem}, S::ZZMatrix, T::ZZMatrix)
  for i in 1:size(T, 1)
    for j in 1:size(S, 1)
      @inbounds Tuv = dot(T[i, :], u) * dot(S[j, :], v)
      @inbounds Tvu = dot(T[i, :], v) * dot(S[j, :], u)
      if Tuv != Tvu
        return Tuv < Tvu
      end
    end
  end
  return false
end

# returns the multiple m for all terms q in p with lm * m = q.
function divides_walk(p::MPolyRingElem, lm::MPolyRingElem, S::MPolyRing)
  div = false
  newpoly = MPolyBuildCtx(S)
  for term in terms(p)
    (b, c) = divides(term, lm)
    if b
      push_term!(
        newpoly, first(coefficients(c)), first(AbstractAlgebra.exponent_vectors(c))
      )
      div = true
    end
  end
  return finish(newpoly), div
end

# returns p reduced by the Groebner basis G w.r.t. the leading terms Lm.
function reduce_walk(
  p::MPolyRingElem, G::Vector{T}, Lm::Vector{T}, ord::MonomialOrdering
) where {T<:MPolyRingElem}
  for i in 1:length(G)
    (q, b) = divides_walk(p, Lm[i], parent(p))
    if b
      return reduce_walk(p - (q * G[i]), G, Lm, ord)
    end
  end
  return p
end

# this function interreduces the Groebner basis G w.r.t. the leading terms Lm with tail-reduction.
function interreduce(
  G::Vector{T}, Lm::Vector{T}, ord::MonomialOrdering
) where {T<:MPolyRingElem}
  for i in 1:Singular.length(G)
    G[i] = reduce_walk(G[i], G[1:end .!= i], Lm[1:end .!= i], ord)
  end
  return G
end

###############################################################
# Several Procedures for the Groebner Walk
###############################################################

# multiplies every entry of the given weight w with 0.1 as long as it stays on the same halfspace as w.
function truncw(G::Oscar.IdealGens, w::Vector{Int}, inw::Vector{T}) where {T<:MPolyRingElem}
  while !checkInt32(w)
    for i in 1:length(w)
      w[i] = round(w[i] * 0.10)
    end
    w = convert_bounding_vector(w)
    if inw != initial_forms(G, w)
      # initials are different - return unrounded weight
      return w, false
    end
  end
  # converted to Vector w of the same face
  return w, true
end

# computes a p-perturbed vector from the matrix M.
function perturbed_vector(G::Oscar.IdealGens, M::ZZMatrix, p::Int)
  n = size(M, 1)
  rows = [M[i,:] for i in 1:p]

  m = maximum.(Ref(abs), rows)
  m_sum = sum(m)
  max_deg = maximum(total_degree.(G)) # TODO: I think this is total degree

  e = max_deg * m_sum + 1
  w = M[1, :] * e^(p - 1)
  for i in 2:p
    w += e^(p - i) * M[i, :]
  end

  w = sum(rows .* (p .- Vector(1:p)))

  return convert_bounding_vector(w)
end

# returns 'true' if the leading terms of G w.r.t the matrixorder T are the same as the leading terms of G w.r.t the weighted monomial order with weight vector t and matrix T.
#function inCone(G::Oscar.IdealGens, T::Union{Matrix{N}, MatElem{N}}, t::Vector{Int})
#    R = change_order(G.base_ring, T)
#    I = Singular.Ideal(R, [change_ring(x, R) for x in gens(G)])
#    cvzip = zip(Singular.gens(I), initials(R, Singular.gens(I), t))
#    for (g, ing) in cvzip
#        if !isequal(Singular.leading_exponent_vector(g), Singular.leading_exponent_vector(ing))
#            return false
#        end
#    end
#    return true
#end

# returns 'true' if the leading terms of G w.r.t the matrixordering T are the same as the leading terms of G w.r.t the weighted monomial order with weight vector t and the matrix order T.
#function inCone(G::Oscar.IdealGens, t::Vector{Int})
#    cvzip = zip(Singular.gens(G), initials(base_ring(G), gens(G), t))
#    for (g, ing) in cvzip
#        if !isequal(leading_exponent_vector(g), leading_exponent_vector(ing))
#            return false
#        end
#    end
#    return true
#end

# returns 'true' if the leading terms of G w.r.t the matrixordering T are the same as the leading terms of G with the current ordering.
same_cone(G::Oscar.IdealGens, T::MonomialOrdering) = all(leading_term.(G; ordering=T) .== leading_term.(G; ordering=ordering(G)))

# TODO: Actual docstring
#= Lifting step from Proposition 3.2 of "The generic Gröbner walk" (Fukuda et al., 2005)
  
  Input:  - G , the reduced G.B of I w.r.t current
          - current, the current monomial order
          - H, the reduced G.B of inw(I) w.r.t the next weight vector w 
          - target, the next monomial order 
  
  Output: - an inclusion minimal G.B of target obtained by subtracting normal forms

  QUESTION: why do we need "target" in Oscar.IdealGens(...)? 

  COMMENT: I think "target" is inappropriately named. It is rather "next_ordering" (i.e the target order, refined by w)
=#
function lift(
  G::Oscar.IdealGens, # momentane GB
  current::MonomialOrdering,
  H::Oscar.IdealGens, # soll GB von initial forms sein
  target::MonomialOrdering,
)
  G = Oscar.IdealGens(
    [
      gen - Oscar.IdealGens(
        [reduce(gen, gens(G); ordering=current, complete_reduction=true)], target
      )[1] for gen in gens(H)
    ],
    target;
    isGB=true,
  )

  return G
end

# lifts the Groebner basis G to the Groebner basis w.r.t. the Ring Rn like it´s done in Collart et al. (1997).
# FIXME: Needs to be fixed.
function liftGW2(
  G::Oscar.IdealGens,
  orderingAlt::MonomialOrdering,
  inG::MPolyIdeal{T},
  H::Oscar.IdealGens,
  ordering::MonomialOrdering,
) where {T<:MPolyRingElem}
  H = gens(H)
  R = base_ring(G)
  G = gens(G)
  for i in 1:length(H)
    q = divrem(gens(Oscar.IdealGens([H[i]], orderingAlt)), gens(inG))
    H[i] = R(0)
    for j in 1:length(gens(inG))
      println(gens(inG)[j])
      println(G[j])
      H[i] = H[i] + q[1][1][j] * G[j]
    end
  end
  return Oscar.IdealGens(filter(!iszero, H), ordering)
end

# divisionalgorithm that returns q with q_1*f_1 + ... + q_s *f_s=p.
function division_algorithm(p::T, f::Vector{T}, R::MPolyRing) where {T<:MPolyRingElem}
  q = Array{Singular.elem_type(R),1}(undef, length(f))
  for i in 1:length(f)
    q[i] = R(0)
  end
  while !isequal(p, R(0))
    i = 1
    div = false
    while (div == false && i <= length(f))
      b, m = divides(leading_term(p), leading_term(f[i]))
      if b
        q[i] = q[i] + m
        p = p - (m * f[i])
        div = true
      else
        i = i + 1
      end
    end
    if div == false
      p = p - leading_term(p)
    end
  end
  return q
end

# converts a vector wtemp by dividing the entries with gcd(wtemp).
function convert_bounding_vector(w::Vector{T}) where {T<:Union{ZZRingElem, QQFieldElem}}
  c = gcd(w)
  if c != 0
    return ZZ.(floor.(w//c))
  else
    w
  end
end

# returns a copy of the PolynomialRing I, equipped with the ordering weight_ordering(cw)*matrix_ordering(T).
create_ordering(R::MPolyRing, cw::Vector{L}, T::Matrix{Int}) where {L<:Number} = weight_ordering(cw, matrix_ordering(R, T))

# interreduces the Groebner basis G. 
# each element of G is replaced by its normal form w.r.t the other elements of G and the current monomial order 
# TODO reference, docstring
# interreduces the Groebner basis G.
function interreduce_walk(G::Oscar.IdealGens)
  generators = collect(gens(G))

  I = 0
  for i in 1:length(gens(G))
    generators[i] = reduce(
      generators[i], generators[1:end .!= i]; ordering=G.ord, complete_reduction=true
    )
  end
  return Oscar.IdealGens(generators, G.ord)
end

#############################################
# unspecific help functions
#############################################

function ident_matrix(n::Int)
  M = zeros(Int, n, n)
  for i in 1:n
    M[i, i] = 1
  end
  return M
end

function anti_diagonal_matrix(n::Int)
  M = zeros(Int, n, n)
  for i in 1:n
    M[i, n + 1 - i] = -1
  end
  return M
end

# Singular.isequal depends on order of generators
function equalitytest(G::Oscar.IdealGens, K::Oscar.IdealGens)
  if length(gens(G)) != length(gens(K))
    return false
  end
  generators = Singular.gens(G)
  count = 0
  for gen in generators
    for r in Singular.gens(K)
      if gen * first(coefficients(leading_term(r; ordering=G.ord))) -
         r * first(coefficients(leading_term(gen; ordering=G.ord))) == 0
        count += 1
        break
      end
    end
  end
  if count == length(gens(G))
    return true
  end
  return false
end

function ordering_as_matrix(w::Vector{Int}, ord::Symbol)
  if length(w) > 2
    if ord == :lex
      return [
        w'
        ident_matrix(length(w))[1:(length(w) - 1), :]
      ]
    end
    if ord == :deglex
      return [
        w'
        ones(Int, length(w))'
        ident_matrix(length(w))[1:(length(w) - 2), :]
      ]
    end
    if ord == :degrevlex
      return [
        w'
        ones(Int, length(w))'
        anti_diagonal_matrix(length(w))[1:(length(w) - 2), :]
      ]
    end
    if ord == :revlex
      return [
        w'
        anti_diagonal_matrix(length(w))[1:(length(w) - 1), :]
      ]
    end
  else
    error("not implemented")
  end
end

function change_weight_vector(w::Vector{Int}, M::Matrix{Int})
  return [
    w'
    M[2:length(w), :]
  ]
end

function insert_weight_vector(w::Vector{Int}, M::Matrix{Int})
  return [
    w'
    M[1:(length(w) - 1), :]
  ]
end

function add_weight_vector(w::Vector{Int}, M::Matrix{Int})
  return [
    w'
    M
  ]
end

add_weight_vector(w::Vector{ZZRingElem}, M::ZZMatrix) = ZZMatrix(add_weight_vector(Int.(w), Matrix{Int}(M)))

function ordering_as_matrix(ord::Symbol, nvars::Int)
  if ord == :lex
    return ident_matrix(nvars)
  end
  if ord == :deglex
    return [
      ones(Int, nvars)'
      ident_matrix(nvars)[1:(nvars - 1), :]
    ]
  end
  if ord == :degrevlex
    return [
      ones(Int, nvars)'
      anti_diagonal_matrix(nvars)[1:(nvars - 1), :]
    ]
  end
  if ord == :revlex
    return [
      w'
      anti_diagonal_matrix(length(w))[1:(length(w) - 1), :]
    ]
  else
    error("not implemented")
  end
end

function checkInt32(w::Vector{Int})
  for i in 1:length(w)
    if tryparse(Int32, string(w[i])) == nothing
      return false
    end
  end
  return true
end

# computes the representation of the matrixorder defined by T.
function representation_vector(G::Oscar.IdealGens, T::Matrix{Int})
  n = size(T)[1]
  M = maximum(T)
  
  d0 = maximum(g -> deg(g, n), gens(g))

  d = M * (2 * d0^2 + (n + 1) * d0)
  w = d^(n - 1) * T[1, :]
  for i in 2:n
    w = w + d^(n - i) * T[i, :]
  end
  return w
end

# checks if Gw is an initialform of an ideal corresponding to a face of the Groebner fan with dimension less than n-1.
function inSeveralCones(Gw::Vector{T}) where {T<:MPolyRingElem}
  counter = 0
  for g in Gw
    if size(collect(Singular.coefficients(g)))[1] > 2
      return true
    end
    if size(collect(Singular.coefficients(g)))[1] == 2
      counter = counter + 1
    end
  end
  if counter > 1
    return true
  end
  return false
end


#Create a copy of the "lift" function (to export it without conflicts)

function lift2(
  G::Oscar.IdealGens, # momentane GB
  current::MonomialOrdering,
  H::Oscar.IdealGens, # soll GB von initial forms sein
  target::MonomialOrdering,
)
  G = Oscar.IdealGens(
    [
      gen - Oscar.IdealGens(
        [reduce(gen, gens(G); ordering=current, complete_reduction=true)], target
      )[1] for gen in gens(H)
    ],
    target;
    isGB=true,
  )

  return G
end