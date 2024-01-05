export annihilator
export is_projective


function annihilator(M::SubquoModule)
  R = base_ring(M)
  F = FreeMod(R, 1)
  I = ideal(R, [one(R)])
  for v in gens(M)
    h = hom(F, M, [v])
    K, _ = kernel(h)
    # this is a hack, because getindex is broken for K
    g = Oscar.as_matrix(F, ambient_representatives_generators(K))
    I = intersect(I, ideal(R, minors(g, 1)))
  end
  return I
end

@doc raw"""
    is_projective(M::SubquoModule)

Given a subquotient ``M = (A + B)/B`` over a ring ``R`` return a triple 
`(success, π, σ)` where `success` is `true` or `false` depending on 
whether or not ``M`` is projective and maps ``π : Rʳ ↔ M : σ`` 
where ``π`` is a projection onto ``M`` and ``σ`` a section of ``π`` 
so that ``Rʳ ≅ M ⊕ N`` splits as a direct sum via the projector 
``σ ∘ π``.
"""
function is_projective(M::SubquoModule; check::Bool=true)
  ### Assemble a presentation of M over its base ring
  # TODO: Eventually replace by the methods for free 
  # presentation from the module package. 
  R = base_ring(M)
  r = ngens(M)
  Rr = FreeMod(R, r)
  proj = hom(Rr, M, gens(M))

  if check
    iszero(annihilator(M)) || return false, proj, nothing
  end

  K, inc = kernel(proj)
  s = ngens(K)
  Rs = FreeMod(R, s)
  a = compose(hom(Rs, K, gens(K)), inc)
  # This is the presentation matrix
  A = matrix(a)

  success, P, k = _is_projective_without_denominators(A)
  !success && return false, hom(Rr, M, [zero(M) for x in gens(Rr)]), hom(M, Rr, [zero(Rr) for x in gens(M)])
  k == 0 || error("invalid numerator required")
  return true, hom(Rr, M, gens(M)), hom(M, Rr, [sum([P[i, j]*Rr[i] for i in 1:ngens(Rr)]) for j in 1:ncols(P)])
end

function is_projective(M::FreeMod; check::Bool=true)
  return true, identity_map(M), identity_map(M)
end
###
# Given a presentation matrix A for a module M over a ring R, 
# this procedure checks whether M is projective and in the affirmative 
# case returns a triple (true, Q, k) where 1//(unit^k) * Q is the 
# projector. 
#
# The optional argument `unit` indicates that the ring R was obtained 
# as the localization R = A[unit⁻¹] from another ring A. 
# The output is then returned in a way that Q lifts to A. 
#
# Note: It is highly advised to feed this function only matrices 
# whose entries do not have denominators. Otherwise, the program will 
# most likely run very slow.
function _is_projective_without_denominators(A::MatElem; 
    unit::RingElem=one(base_ring(A)),
    task::Symbol=:with_projector
  )
  R = base_ring(A)
  m = nrows(A)
  n = ncols(A)

  # The condition for end of recursion:
  if iszero(A) 
    return true, identity_matrix(R, n), 0
  end

  entry_list = [(A[i,j], i, j) for i in 1:m for j in 1:n if !iszero(A[i,j])]
  I = ideal(R, [a[1] for a in entry_list])
  if !(one(R) in I)
    # If I is not the unit ideal, it might still be the case that this holds true for the 
    # restriction to the connected components of Spec(R).
    R isa Union{MPolyRing, MPolyQuoRing, MPolyLocRing, MPolyQuoLocRing} || error("method not implemented")
    # TODO: Work this out without schemes on the purely algebraic side.
    # This is a temporary hotfix to address a particular boundary case, see #1882. 
    # The code below is also not generic and should eventually be adjusted.

    X = Spec(R)
    U = connected_components(X)
    l = length(U)
    l == 1 && return false, zero_matrix(R, n, n), 0
    projectors = MatElem[] # The local projectors on each component
    expon = Int[]
    for k in 1:l
      L, inc = localization(R, complement_equation(U[k])) # We can not use OO(U[k]) directly, because
                                                          # we'd be missing the map in that case. 
      success, p, k = _is_projective_without_denominators(change_base_ring(L, A), unit=L(unit), task=task)
      !success && return false, zero_matrix(R, n, n), 0
      q = zero_matrix(R, nrows(p), ncols(p))
      for i in 1:nrows(p)
        for j in 1:ncols(p)
          q[i, j] = preimage(inc, p[i, j])
        end
      end
      push!(projectors, q)
      push!(expon, k)
    end

    c = coordinates(one(R), ideal(R, [complement_equation(U[k])^(expon[k]+1) for k in 1:length(expon)]))
    result = sum([c[k]*projectors[k] for k in 1:length(projectors)])

    # Copied from below
    d = lcm(_lifted_denominator.(result))
    if isone(d)
      return true, result, 0
    end

    inner, outer = ppio(d, _lifted_numerator(unit))
    (mpow, upow) = Oscar._minimal_power_such_that(_lifted_numerator(unit), 
                                                  x->(divides(x, inner)[1]))

    # pull u^mpow from each entry of the matrix:
    L = zero_matrix(R, n, n)
    for i in 1:n
      for j in 1:n
        L[i, j] = R(_lifted_numerator(result[i, j])*
                    divides(upow, inner)[2]*
                    divides(d, _lifted_denominator(result[i, j]))[2]
                   )*inv(R(outer*(_lifted_denominator(unit)^mpow)))
      end
    end
    return true, L, mpow
  end


  lambda1 = coordinates(one(R), I) # coefficients for the first powers
  involved_entries = [k for k in 1:length(lambda1) if !iszero(lambda1[k])]

  # For every 'involved entry' u = aᵢⱼ, localize at u and form 
  # the submatrix B = Bᵢⱼ. For B we obtain the projectors Pᵢⱼ 
  # by recursion. It will, in general have a denominator uᵏ for 
  # some k and is returned as a pair (Qᵢⱼ, k) for the given 
  # unit u and a matrix Qᵢⱼ ∈ Rⁿˣⁿ.
  
  # Get rid of the zero entries.
  entry_list = entry_list[involved_entries]
  lambda1 = lambda1[1, involved_entries]
  sub_results = Tuple{Bool, <:MatrixElem, Int}[]
  sub_localizations = Tuple{<:Ring, <:Map}[]
  for (u, i, j) in entry_list
    # Assemble the submatrix Bᵢⱼ
    B = copy(A)
    for k in 1:i-1
      multiply_row!(B, u, k)
      add_row!(B, -A[k, j], i, k)
    end
    for k in i+1:nrows(B)
      multiply_row!(B, u, k)
      add_row!(B, -A[k, j], i, k)
    end
    Asub = B[[k for k in 1:m if k != i], [k for k in 1:n if k !=j]]
    Rloc, inc = localization(R, u)
    push!(sub_localizations, (Rloc, inc))
    # We expect a pair (Q, k) consisting of a matrix Q defined over Rloc, 
    # but liftable to a matrix over R without effort. The local projector 
    # over Rloc is then 1//uᵏ ⋅ Q. 
    push!(sub_results, _is_projective_without_denominators(map_entries(Rloc, Asub), unit=u, task=task))
    B = last(sub_results)[2]
    k = last(sub_results)[3]
    if last(sub_results)[1] == false
      return false, zero_matrix(R, n, n), 0
    end
  end
  powers = [k for (_, _, k) in sub_results]
  # Find the coefficients cᵣ so that ∑ᵣ cᵣ ⋅ uᵣᵏ⁽ʳ⁾ = 1 ∈ R.
  # TODO: This can be optimized
  if task == :with_projector
    c = coordinates(one(R), ideal(R, [u^(k+1) for ((u, _, _), k) in zip(entry_list, powers)]))

    result = zero_matrix(R, n, n)
    for ((u, i, j), (_, Q, k), (Rloc, inc), lambda) in zip(entry_list, sub_results, sub_localizations, c)
      # Lift the matrix Q to an n×n-matrix over R
      Qinc = zero_matrix(R, n, n)
      for r in 1:j-1
        for s in 1:j-1
          Qinc[r, s] = preimage(inc, Q[r, s])
        end
        for s in j+1:n
          Qinc[r, s] = preimage(inc, Q[r, s-1])
        end
      end
      for r in j+1:n
        for s in 1:j-1
          Qinc[r, s] = preimage(inc, Q[r-1, s])
        end
        for s in j+1:n
          Qinc[r, s] = preimage(inc, Q[r-1, s-1])
        end
      end

      # The matrix Qinc needs to be fed with a vector v whose j-th entry 
      # has been deleted via the localization at u. This can only be done 
      # over R for u⋅v, so we multiply by u, keeping in mind that this will 
      # become a unit in the localization. 
      P = u*identity_matrix(R, n)
      P[j, j] = 0
      for l in 1:n
        l == j && continue
        P[j, l] = -A[i, l]
      end
      result = result + lambda*P*Qinc
    end

    # pull the denominator u from the result and make the matrix liftable. 
    d = reduce(lcm, _lifted_denominator.(result))
    if isone(d)
      return true, result, 0
    end

    inner, outer = ppio(d, _lifted_numerator(unit))
    (mpow, upow) = Oscar._minimal_power_such_that(_lifted_numerator(unit), 
                                                  x->(divides(x, inner)[1]))

    # pull u^mpow from each entry of the matrix:
    L = zero_matrix(R, n, n)
    for i in 1:n
      for j in 1:n
        L[i, j] = R(_lifted_numerator(result[i, j])*
                    divides(upow, inner)[2]*
                    divides(d, _lifted_denominator(result[i, j]))[2]
                   )*inv(R(outer*(_lifted_denominator(unit)^mpow)))
      end
    end
    return true, L, mpow
  elseif task==:without_projector 
    return true, zero_matrix(R, n, n), 0
  end
  error("task could not be identified")
end

# helper functions to unify the interface
_lifted_denominator(a::RingElem) = one(parent(a))
_lifted_denominator(a::AbsLocalizedRingElem) = denominator(a)
_lifted_denominator(a::MPolyQuoRingElem) = one(base_ring(parent(a)))
_lifted_denominator(a::MPolyQuoLocRingElem) = lifted_denominator(a)
_lifted_numerator(a::RingElem) = a
_lifted_numerator(a::AbsLocalizedRingElem) = numerator(a)
_lifted_numerator(a::MPolyQuoRingElem) = lift(a)
_lifted_numerator(a::MPolyQuoLocRingElem) = lifted_numerator(a)

########################################################################
# Various localization routines for localizing at powers of elements   #
#                                                                      #
# This deserves special constructors, because we can deliver maps for  # 
# lifting which is not possible in general.                            #
########################################################################
function localization(A::MPolyQuoRing, f::MPolyQuoRingElem)
  R = base_ring(A)
  U = MPolyPowersOfElement(R, [lift(f)])
  W = MPolyLocRing(R, U)
  L = MPolyQuoLocRing(R, modulus(A), U, A, W)
  function func(a::MPolyQuoRingElem)
    parent(a) == A || error("element does not belong to the correct ring")
    return L(a, check=false)
  end
  function func_inv(a::MPolyQuoLocRingElem{<:Any, <:Any, <:Any, <:Any, 
                                                 <:MPolyPowersOfElement}
    )
    L == parent(a) || error("element does not belong to the correct ring")
    iszero(numerator(a)) && return zero(A)
    isone(lifted_denominator(a)) && return A(lifted_numerator(a))
    success, c = divides(numerator(a), denominator(a))
    if !success 
      error("lifting not possible")
    end
    return c
  end
  return L, MapFromFunc(A, L, func, func_inv)
end

function localization(A::MPolyLocRing, f::MPolyLocRingElem)
  R = base_ring(A)
  d = numerator(f)
  U = MPolyPowersOfElement(R, [d])
  L = MPolyLocRing(R, U*inverted_set(A))
  function func(a::MPolyLocRingElem)
    parent(a) == A || error("element does not belong to the correct ring")
    return L(a, check=false)
  end
  function func_inv(a::MPolyLocRingElem{<:Any, <:Any, <:Any, <:Any, 
                                              <:MPolyPowersOfElement}
    )
    L == parent(a) || error("element does not belong to the correct ring")
    isone(denominator(a)) && return A(numerator(a))
    iszero(numerator(a)) && return zero(A)
    i, o = ppio(denominator(a), d)
    return A(divexact(numerator(f), i), o, check=false)
  end
  return L, MapFromFunc(A, L, func, func_inv)
end

function localization(A::MPolyQuoLocRing, f::MPolyQuoLocRingElem)
  R = base_ring(A)
  d = lifted_numerator(f)
  U = MPolyPowersOfElement(R, [d])
  L = MPolyQuoLocRing(R, modulus(underlying_quotient(A)), U*inverted_set(A))
  function func(a::MPolyQuoLocRingElem)
    parent(a) == A || error("element does not belong to the correct ring")
    return L(a)
  end
  function func_inv(a::MPolyQuoLocRingElem{<:Any, <:Any, <:Any, <:Any, 
                                              <:MPolyPowersOfElement}
    )
    L == parent(a) || error("element does not belong to the correct ring")
    isone(lifted_denominator(a)) && return A(lifted_numerator(a))
    iszero(numerator(a)) && return zero(A)
    i, o = ppio(lifted_denominator(a), d)
    success, c = divides(numerator(a), parent(numerator(a))(i))
    if !success 
      # last resort:
      return A(lifted_numerator(a), lifted_denominator(a))
    end
    return A(lift(c), o, check=false)
  end
  return L, MapFromFunc(A, L, func, func_inv)
end

