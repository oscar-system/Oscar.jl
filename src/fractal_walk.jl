###############################################################
# Plain version of the Fractal Walk.
# This version checks if an entry of an intermediate weight vector is bigger than int32. In case of that the Buchberger-Algorithm is used to compute the Groebner basis of the ideal of the initialforms.
###############################################################

function fractal_walk(G::Oscar.IdealGens, start::MonomialOrdering, target::MonomialOrdering)
  @vprintln :groebner_walk "FractalWalk_standard results"
  @vprintln :groebner_walk "Crossed Cones in: "

  S = canonical_matrix(start)
  T = canonical_matrix(target)

  pTargetWeights = [perturbed_vector(G, T, i) for i in 1:nvars(base_ring(G))]
  return fractal_recursion_step(G, S, T, S[1, :], pTargetWeights, 1)
end

function fractal_recursion_step(
  G::Oscar.IdealGens,
  S::ZZMatrix,
  T::ZZMatrix,
  current_weight::Vector{ZZRingElem},
  pTargetWeights::Vector{Vector{ZZRingElem}},
  p::Int
)
  R = base_ring(G)
  G.isGB = true
  w = current_weight
  old_ordering = ordering(G)
  while true
    t = next_weight_fractal_walk(G, w, pTargetWeights[p])

    # Handling the final step in the current depth.
    # next_weight_fractal_walk may return 0 if the target vector does not lie in the cone of T while G already defines the Groebner basis w.r.t. T.
    # -> Checking if G is already a Groebner basis w.r.t. T solves this problem and reduces computational effort since next_weight_fractal_walk returns 1 in the last step on every local path.        if t == 1 && p != 1
    if t == 1 && p != 1
      if same_cone(G, T)
        @vprintln :groebner_walk ("depth $p: in cone ", current_weight, ".")...

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
        deepcopy(current_weight),
        pTargetWeights,
        p + 1,
      )
    end
    #H = liftGW2(G, R, Gw, H, Rn)
    H = lift_fractal_walk(G, H, ordNew)
    G = autoreduce(H)
    ordAlt = ordNew
    current_weight = w
  end
end

