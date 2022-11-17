export is_projective, annihilator

function annihilator(M::SubQuo)
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

iszero(I::Ideal) = all(x->(iszero(x)), gens(I))

@Markdown.doc """
    is_projective(M::SubQuo)

Given a subquotient ``M = (A + B)/B`` over a ring ``R`` return a triple 
`(success, π, σ)` where `success` is `true` or `false` depending on 
whether or not ``M`` is projective and maps ``π : Rʳ ↔ M : σ`` 
where ``π`` is a projection onto ``M`` and ``σ`` a section of ``π`` 
so that ``Rʳ ≅ M ⊕ N`` splits as a direct sum via the projector 
``σ ∘ π``.
"""
function is_projective(M::SubQuo; check::Bool=true)
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
  !success && return false, zero_map(Rr, M), zero_map(M, Rr)
  k == 0 || error("invalid numerator required")
  return true, hom(Rr, M, gens(M)), hom(M, Rr, [sum([P[i, j]*Rr[i] for i in 1:ngens(Rr)]) for j in 1:ncols(P)])
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
function _is_projective_without_denominators(A::MatElem; unit::RingElem=one(base_ring(A)))
  R = base_ring(A)
  m = nrows(A)
  n = ncols(A)

  # The condition for end of recursion:
  if iszero(A) 
    return true, one(MatrixSpace(R, n, n)), 0
  end

  entry_list = [(A[i,j], i, j) for i in 1:m for j in 1:n if !iszero(A[i,j])]
  I = ideal(R, [a[1] for a in entry_list])
  if !(one(R) in I)
    return false, zero(MatrixSpace(R, n, n)), 0
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
    Rloc, inc = Localization(R, u)
    push!(sub_localizations, (Rloc, inc))
    # We expect a pair (Q, k) consisting of a matrix Q defined over Rloc, 
    # but liftable to a matrix over R without effort. The local projector 
    # over Rloc is then 1//uᵏ ⋅ Q. 
    push!(sub_results, _is_projective_without_denominators(map_entries(Rloc, Asub), unit=u))
    B = last(sub_results)[2]
    k = last(sub_results)[3]
    if last(sub_results)[1] == false
      return false, zero(MatrixSpace(R, n, n)), 0
    end
  end
  powers = [k for (_, _, k) in sub_results]
  # Find the coefficients cᵣ so that ∑ᵣ cᵣ ⋅ uᵣᵏ⁽ʳ⁾ = 1 ∈ R.
  # TODO: This can be optimized
  c = coordinates(one(R), ideal(R, [u^(k+1) for ((u, _, _), k) in zip(entry_list, powers)]))

  result = zero(MatrixSpace(R, n, n))
  for ((u, i, j), (_, Q, k), (Rloc, inc), lambda) in zip(entry_list, sub_results, sub_localizations, c)
    # Lift the matrix Q to an n×n-matrix over R
    Qinc = zero(MatrixSpace(R, n, n))
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
    P = u*one(MatrixSpace(R, n, n))
    P[j, j] = 0
    for l in 1:n
      l == j && continue
      P[j, l] = -A[i, l]
    end
    result = result + lambda*P*Qinc
  end

  # pull the denominator u from the result and make the matrix liftable. 
  d = lcm(_lifted_denominator.(result))
  if isone(d)
    return true, result, 0
  end

  inner, outer = ppio(d, _lifted_numerator(unit))
  (mpow, upow) = Oscar._minimal_power_such_that(_lifted_numerator(unit), 
                                                x->(divides(x, inner)[1]))

  # pull u^mpow from each entry of the matrix:
  L = zero(MatrixSpace(R, n, n))
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

# helper functions to unify the interface
_lifted_denominator(a::RingElem) = one(parent(a))
_lifted_denominator(a::AbsLocalizedRingElem) = denominator(a)
_lifted_denominator(a::MPolyQuoElem) = one(base_ring(parent(a)))
_lifted_denominator(a::MPolyQuoLocalizedRingElem) = lifted_denominator(a)
_lifted_numerator(a::RingElem) = a
_lifted_numerator(a::AbsLocalizedRingElem) = numerator(a)
_lifted_numerator(a::MPolyQuoElem) = lift(a)
_lifted_numerator(a::MPolyQuoLocalizedRingElem) = lifted_numerator(a)

########################################################################
# Various localization routines for localizing at powers of elements   #
#                                                                      #
# This deserves special constructors, because we can deliver maps for  # 
# lifting which is not possible in general.                            #
########################################################################
function Localization(A::MPolyQuo, f::MPolyQuoElem)
  R = base_ring(A)
  U = MPolyPowersOfElement(R, [lift(f)])
  W = MPolyLocalizedRing(R, U)
  L = MPolyQuoLocalizedRing(R, modulus(A), U, A, W)
  function func(a::MPolyQuoElem)
    parent(a) == A || error("element does not belong to the correct ring")
    return L(a, check=false)
  end
  function func_inv(a::MPolyQuoLocalizedRingElem{<:Any, <:Any, <:Any, <:Any, 
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
  return L, MapFromFunc(func, func_inv, A, L)
end

function Localization(A::MPolyLocalizedRing, f::MPolyLocalizedRingElem)
  R = base_ring(A)
  d = numerator(f)
  U = MPolyPowersOfElement(R, [d])
  L = MPolyLocalizedRing(R, U*inverted_set(A))
  function func(a::MPolyLocalizedRingElem)
    parent(a) == A || error("element does not belong to the correct ring")
    return L(a, check=false)
  end
  function func_inv(a::MPolyLocalizedRingElem{<:Any, <:Any, <:Any, <:Any, 
                                              <:MPolyPowersOfElement}
    )
    L == parent(a) || error("element does not belong to the correct ring")
    isone(denominator(a)) && return A(numerator(a))
    iszero(numerator(a)) && return zero(A)
    i, o = ppio(denominator(a), d)
    return A(divexact(numerator(f), i), o, check=false)
  end
  return L, MapFromFunc(func, func_inv, A, L)
end

function Localization(A::MPolyQuoLocalizedRing, f::MPolyQuoLocalizedRingElem)
  R = base_ring(A)
  d = lifted_numerator(f)
  U = MPolyPowersOfElement(R, [d])
  L = MPolyQuoLocalizedRing(R, modulus(underlying_quotient(A)), U*inverted_set(A))
  function func(a::MPolyQuoLocalizedRingElem)
    parent(a) == A || error("element does not belong to the correct ring")
    return L(a)
  end
  function func_inv(a::MPolyQuoLocalizedRingElem{<:Any, <:Any, <:Any, <:Any, 
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
  return L, MapFromFunc(func, func_inv, A, L)
end

