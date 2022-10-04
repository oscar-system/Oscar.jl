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

  success, P = _is_projective(A)
  !success && return false, zero_map(Rr, M), zero_map(M, Rr)
  return true, hom(Rr, M, gens(M)), hom(M, Rr, [sum([P[i, j]*Rr[i] for i in 1:ngens(Rr)]) for j in 1:ncols(P)])
end

### This function assumes that the entry aᵢⱼ is a unit in 𝒪(X).
# Row reduction is performed, then the projector is computed 
# recursively for the lower block matrix.
# The return value is aᵢⱼ ⋅ P for the actual projector P.
function _is_projective(A::MatrixElem, i::Int, j::Int)
  R = base_ring(A)

  m = nrows(A)
  n = ncols(A)

  ### perform the row-reduction for the given unit entry
  B = copy(A)
  u = A[i, j]
#  uinv = inv(u)
#  multiply_row!(B, uinv, i)
#  for k in 1:i-1
#    add_row!(B, -A[k, j], i, k)
#  end
#  for k in i+1:nrows(B)
#    add_row!(B, -A[k, j], i, k)
#  end
  for k in 1:i-1
    multiply_row!(B, u, k)
    add_row!(B, -A[k, j], i, k)
  end
  for k in i+1:nrows(B)
    multiply_row!(B, u, k)
    add_row!(B, -A[k, j], i, k)
  end
  Asub = B[[k for k in 1:m if k != i], [k for k in 1:n if k !=j]]

  ### Assemble aᵢⱼ times the projector
  P = u*one(MatrixSpace(R, n, n))
  P[j, j] = 0
  for l in 1:n
    l == j && continue
    P[j, l] = -B[i, l]
  end

  ### Harvest the projector on the subspace
  success, Psub = _is_projective(Asub)
  # We must assume Psub to have denominators coming 
  # from the coefficients of the linear combination 
  # for 1 ∈ R[aᵢⱼ⁻¹].
  !success && return false, P

  Qinc = zero(P)
  for r in 1:j-1
    for s in 1:j-1
      Qinc[r, s] = Psub[r, s]
    end
    for s in j+1:n
      Qinc[r, s] = Psub[r, s-1]
    end
  end
  for r in j+1:n
    for s in 1:j-1
      Qinc[r, s] = Psub[r-1, s]
    end
    for s in j+1:n
      Qinc[r, s] = Psub[r-1, s-1]
    end
  end

  # Find the minimal power of u = aᵢⱼ such that u^k * Qinc
  # is in the image of R → R[u⁻¹].

  d = lcm(lifted_denominator.(Qinc))
  @show d
  (k, v) = Oscar._minimal_power_such_that(lifted_numerator(u), x->!(divides(d, x)[1]))
  @show k
  @show v

  return true, P*Qinc, v
end

function _is_projective(A::MatElem)
  R = base_ring(A)

  m = nrows(A)
  n = ncols(A)

  ### Checking for one possible end of recursion: The matrix presents 
  # a free module
  if iszero(A)
    return true, one(MatrixSpace(R, n, n))
  end

  entry_list = [(A[i,j], i, j) for i in 1:m for j in 1:n if !iszero(A[i,j])]
  I = ideal(R, [a[1] for a in entry_list])
  one(R) in I || return false, zero(MatrixSpace(R, n, n))

  c = coordinates(one(R), I)
  involved_entries = [k for k in 1:length(c) if !iszero(c[k])]
  if typeof(base_ring(A))<:MPolyQuoLocalizedRing
    @show c
    @show [length(monomials(lifted_denominator(a))) for a in c]
  end

  ### Localize at the hypersurfaces given by the involved entries 
  # and collect the local projectors. 
  rec_results = []
  restriction_maps = []
  for k in involved_entries
    L, res = Localization(R, entry_list[k][1])
    AU = map_entries(res, A)
    push!(rec_results, _is_projective(AU, entry_list[k][2], entry_list[k][3]))
    push!(restriction_maps, res)
    if last(rec_results)[1] == false
      return false, zero(MatrixSpace(R, n, n))
    end
  end

  ### Lift all the local projectors to the top ring and return their 
  # linear combination.
  Q = zero(MatrixSpace(R, n, n))
  v = [b for (_, _, b) in rec_results]
  @show length(involved_entries)
  @show length(v)
  # This must exist but can probably be tuned:
  @show "finding coefficients"
  cnew = coordinates(one(R), ideal(R, v))
  @show "done"
  for (k, l, (_, Psub, b), res) in zip(involved_entries, cnew, rec_results, restriction_maps)
    u = entry_list[k][1]
    @show k
    @show v
    w = divides(R(b), R(u))[2]
    Q += l * map_entries(x->preimage(res, x), base_ring(Psub)(w)*Psub)
  end

  return true, Q
end

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
      @show denominator(a)
      @show total_degree(denominator(a))
      @show length(monomials(denominator(a)))
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
  L = MPolyQuoLocalizedRing(R, modulus(quotient_ring(A)), U*inverted_set(A))
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
    @show total_degree(lifted_denominator(a))
    @show length(monomials(lifted_denominator(a)))
    @show d
    i, o = ppio(lifted_denominator(a), d)
    success, c = divides(numerator(a), parent(numerator(a))(i))
    if !success 
      @show denominator(a)
      @show total_degree(lifted_denominator(a))
      @show length(monomials(lifted_denominator(a)))
      # last resort:
      return A(lifted_numerator(a), lifted_denominator(a))
    end
    return A(lift(c), o, check=false)
  end
  return L, MapFromFunc(func, func_inv, A, L)
end

