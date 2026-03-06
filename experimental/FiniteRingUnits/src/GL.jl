function _presentation_sln(K, n)
  @assert n >= 3
  @assert K isa Oscar.Hecke.fqPolyRepField || Oscar.is_absolute(K)
  p = Int(characteristic(K))
  f = Oscar.defining_polynomial(K)
  m = degree(K)
  a(b, i) = K isa Oscar.Hecke.fqPolyRepField ? Int(Oscar.coeff(b, i)) : Int(lift(ZZ, Oscar.coeff(b, i)))

  symbs = String[]
  _d = Dict()
  l = 1
  for i in 1:n
    for j in 1:n
      if i == j
        continue
      end
      for r in 0:m-1
        push!(symbs, "x_$i,$j$(r == 0 ? "(1)" : "(a^$r)")")
        _d[(i, j, r)] = l
        l += 1
      end
    end
  end
  F = free_group(symbs)
  g = gens(F)
  d = Dict(k => g[_d[k]] for k in keys(_d))
  rels = elem_type(F)[]
  ll = 1
  for i in 1:n
    for j in 1:n
      if i == j
        continue
      end
      for r in 0:m-1
        #@info ll, 1, i, j, r
        push!(rels, d[i, j, r]^p)
        ll += 1
      end
      for r in 0:m-1
        for s in 0:m-1
          #@info ll, 2, i, j, r, s
          push!(rels, Oscar.comm(d[i, j, r], d[i, j, s]))
          ll += 1
        end
      end
      for k in 1:n
        for l in 1:n
          if k == l || j == k || i == l
            continue
          end
          for r in 0:m-1
            for s in 0:m-1
              #@info ll, 3, i, j, k, l, r, s
              push!(rels, Oscar.comm(d[i, j, r], d[k, l, s]))
              ll += 1
            end
          end
        end
      end

      for k in 1:n
        #i, j, k
        if i == j || j == k || i == k
          continue
        end
        for r in 0:m-1
          for s in 0:m-1
            #@info ll, 4, i, j, k, r, s
            b = gen(K)^r * gen(K)^s
            push!(rels, Oscar.comm(d[i, j, r], d[j, k, s]^-1)^-1 * prod(d[i, k, o]^(-a(b, o)) for o in 0:m-1))
            ll += 1
          end
        end
      end
    end
  end
  #@info rels
  filter!(!is_one, rels)
  return F, d, _d, rels
end

function _effective_presentation_of_slnq(GM)
  n = degree(GM)
  if n >= 3
    _effective_presentation_of_slnq_largen(GM)
  elseif n == 2
    _effective_presentation_of_slnq_2(GM)
  else
    error("asds")
  end
end

function _effective_presentation_of_slnq_largen(GM)
  K = base_ring(GM)
  n = degree(GM)
  @assert n >= 3
  m = degree(K)
  #@info "construction presentation"
  F, d, _d, rels = _presentation_sln(K, n)
  #@info "done"
  ge = Vector{Oscar.dense_matrix_type(K)}(undef, length(d))
  a = gen(K)
  I = Oscar.identity_matrix(K, n)
  for (i, j, r) in keys(d)
    for r in 0:m-1
      II = deepcopy(I)
      II[i, j] = a^r
      ge[_d[i, j, r]] = II
    end
  end
  geinv = inv.(ge)
  G, FtoG = quo(F, rels)

  dislog = m -> begin
    es = _write_as_product_of_elementary_matrices(m)
    o = one(G)
    for e in es
      o = o * FtoG(_elementary_matrix_to_gen(e, F, d))
    end
    return o
  end
  return EffectivePresentation(GM, G, x -> dislog(matrix(x)), w -> GM(map_word(w, ge)))
end

function effective_presentation(F::FinField)
  A, AtoF = Oscar.unit_group(F)
  # a is primitive
  AtoG = isomorphism(FPGroup, A)
  G = codomain(AtoG)
  return EffectivePresentation(F, G, x -> AtoG(preimage(AtoF, x)), y -> AtoF(preimage(AtoG, y)))
end

function _effective_presentation_of_glnq(G)
  n = degree(G)
  if n == 1
    return effective_presentation(G)
  end
  K = base_ring(G)
  H = Oscar.SL(n, K)
  A = _effective_presentation_of_slnq(H)
  #@info A
  C = effective_presentation(K)
  #@info C
  f = x -> G(x) # SL -> GL
  g = x -> Oscar.det(Oscar.matrix(x)) # GL -> K^*
  fpreim = x -> begin
    @assert Oscar.det(Oscar.matrix(x)) == 1
    H(Oscar.matrix(x))
  end
  gpreim = x -> begin
    z = Oscar.identity_matrix(K, n)
    z[1, 1] = x
    @assert Oscar.det(z) == x
    G(z)
  end
  return extension(A, C, G, f, fpreim, g, gpreim)
end

###

function _write_as_product_of_elementary_matrices(Nred)
  Nred = deepcopy(Nred)
  N = deepcopy(Nred)
  k = nrows(Nred)
  Nred2 = deepcopy(Nred)

  if !isone(Oscar.det(Nred))
    throw(ArgumentError("Matrix must have determinant one"))
  end

  trafos = typeof(Nred)[]

  for i in 1:k
    Nred, tra = _normalize_column(Nred, i)
    append!(trafos, tra)
  end
  for T in trafos
    Nred2 = T * Nred2
  end
  @assert Nred2 == Nred

  Nredtr = transpose(Nred)

  trafos_tr = typeof(Nred)[]

  for i in 1:k
    Nredtr, tra = _normalize_column(Nredtr, i)
    append!(trafos_tr, tra)
  end
  Nredtr, trafos3 = _normalize_diagonal(Nredtr)
  append!(trafos_tr, trafos3)
  @assert isone(Nredtr)

  res = typeof(Nred)[]

  for T in trafos
    push!(res, _inv_elementary_matrix(T))
  end

  for T in reverse(trafos_tr)
    push!(res, transpose(_inv_elementary_matrix(T)))
  end

  @assert reduce(*, res; init = one(N)) == N
  return res
end

function _elementary_matrix_to_gen(M, F, dd)
  res = one(F)
  n = nrows(M)
  F = base_ring(M)
  d = degree(F)
  for i in 1:n
    for j in 1:n
      if i != j && !iszero(M[i, j])
        # this is e_{ij}(M[i, j])
        for k in 0:d-1
          a = M[i, j]
          if is_zero(Oscar.coeff(a, k))
            continue
          end
          res = res * dd[i, j, k]^(F isa Oscar.fqPolyRepField ? Int(Oscar.coeff(a, k)) : Int(lift(ZZ, Oscar.coeff(a, k))))
        end
        return res
      end
    end
  end
end

function _normalize_diagonal(N)
  n = nrows(N)
  trafos = typeof(N)[]
  R = base_ring(N)
  for i in n:-1:2
    a = N[i, i]
    if is_one(a)
      continue
    end
    inva = inv(a)
    E1 = elementary_matrix(R, n, i - 1, i, -one(R))
    E2 = elementary_matrix(R, n, i, i - 1,  one(R))
    E3 = elementary_matrix(R, n, i - 1, i, -one(R))
    E4 = elementary_matrix(R, n, i - 1, i, a)
    E5 = elementary_matrix(R, n, i, i - 1,  -inva)
    E6 = elementary_matrix(R, n, i - 1, i, a)
    N = E6 * E5 * E4 * E3 * E2 * E1 * N
    push!(trafos, E1)
    push!(trafos, E2)
    push!(trafos, E3)
    push!(trafos, E4)
    push!(trafos, E5)
    push!(trafos, E6)
  end
  @assert isone(N)
  return N, trafos
end

function _inv_elementary_matrix(M)
  n = nrows(M)
  N = Oscar.identity_matrix(base_ring(M), n)
  for i in 1:n
    for j in 1:n
      if i != j && !iszero(M[i, j])
        N[i, j] = -M[i, j]
      end
    end
  end
  @assert isone(N * M)
  return N
end

function _normalize_column(N, i)
  n = nrows(N)
  R = base_ring(N)
  trafos = typeof(N)[]
  if is_unit(N[i, i])
    ainv = inv(N[i, i])
    for j in n:-1:(i + 1)
      if is_zero(N[j, i])
        continue
      end
      E = elementary_matrix(R, n, j, i, -ainv * N[j, i])
      #@show N
      N = Oscar.mul!(N, E, N)
      #@show N
      push!(trafos, E)
    end
    return N, trafos
  else
    for j in (i + 1):n
      if is_unit(N[j, i])
        E1 = elementary_matrix(R, n, i, j, one(R))
        N = Oscar.mul!(N, E1, N)
        push!(trafos, E1)
        E2 = elementary_matrix(R, n, j, i, -one(R))
        N = Oscar.mul!(N, E2, N)
        push!(trafos, E2)
        E3 = elementary_matrix(R, n, i, j, one(R))
        N = Oscar.mul!(N, E3, N)
        push!(trafos, E3)
        @assert is_unit(N[i, i])
        N, trafos2 = _normalize_column(N, i)
        append!(trafos, trafos2)
        return N, trafos
      end
    end
  end
  error("Something went wrong")
end

function elementary_matrix(R, n, i, j, a)
  @assert i != j
  M = Oscar.identity_matrix(R, n)
  M[i, j] = a
  return M
end

####### SL(2, q)
#
# I so pale
#
# SL(2, 2^e)

function _effective_presentation_of_slnq_2(GM)
  K = base_ring(GM)
  if order(K) == 2 || characteristic(K) != 2
    return effective_presentation(GM)
  end
  @assert degree(GM) == 2
  q = order(K)
  FU, mFU = unit_group(K)
  @assert ngens(FU) == 1
  w = mFU(FU[1])
  w2m = Oscar.minpoly(w^2)
  u = i -> Int(lift(ZZ, Oscar.coeff(w2m, i)))
  _x = preimage(mFU, 1 + w^2).coeff[1]
  _y = preimage(mFU, w^2).coeff[1]
  g, o = Oscar.gcdinv(_y, q - 1)
  @assert is_one(g)
  m = Int(mod(_x * o, q - 1))
  @assert w^(2*m) == 1 + w^2
  @assert characteristic(K) == 2
  @assert K isa Oscar.Nemo.fqPolyRepField || Oscar.is_absolute(K)
  e = degree(K)
  F = free_group(["t", "d", "U"])
  t, d, U = gens(F)
  rels = elem_type(F)[]
  push!(rels, (U*t)^3)
  push!(rels, U^2)
  push!(rels, (U*d)^2)
  push!(rels, (t*d)^Int(q-1))
  push!(rels, t^2)
  push!(rels, (t^(d^m))^(-1) * Oscar.comm(t, d))
  push!(rels, prod((t^(u(i)))^(d^i) for i in 0:degree(w2m)))
  G, FtoG = quo(F, rels)

  _t = GM(matrix(K, 2, 2, [1, 1, 0, 1]))
  _d = GM(matrix(K, 2, 2, [w^-1, 0, 0, w]))
  _U = GM(matrix(K, 2, 2, [0, 1, -1, 0]))

  D = _find_standard_elementary_matrix(K, w, _t, _d, _U, t, d, U)

  dislog = m -> begin
    es = _write_as_product_of_elementary_matrices(m)
    o = one(G)
    for e in es
      o = o * FtoG(_elementary_matrix_to_gen(e, F, D))
    end
    return o
  end

  expo = w -> begin
    return map_word(w, [_t, _d, _U])
  end

  return EffectivePresentation(GM, G, x -> dislog(matrix(x)), expo)
end

function _find_standard_elementary_matrix(K, w, _t, _d, _U, t, d, U)
  e = degree(K)
  a = gen(K)
  w2 = w^2
  D = Dict()
  for i in 0:(e - 1)
    o = Oscar.disc_log(w2, a^i)
    #@info _d^-o * _t * _d^o 
    D[1, 2, i] = d^-o * t * d^o
    @assert map_word(D[1, 2, i], [_t, _d, _U]) == parent(_t)(matrix(K, 2, 2, [1, a^i, 0, 1]))

    #@info _d^o * (_U^-1 * _t * _U) * _d^-o
    D[2, 1, i] = d^o * (U^-1 * t * U) * d^-o
    @assert map_word(D[2, 1, i], [_t, _d, _U]) == parent(_t)(matrix(K, 2, 2, [1, 0, a^i, 1]))
  end
  return D
end
