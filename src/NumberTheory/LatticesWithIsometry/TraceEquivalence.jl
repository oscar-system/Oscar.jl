export trace_lattice, inverse_trace_lattice

function _cyclotomic_tower(n::Integer)
  K,a = CyclotomicRealSubfield(n, cached=false)
  Kt, t = PolynomialRing(K, "t", cached=false)
  E,b = number_field(t^2-a*t+1, "b", cached=false)
  set_attribute!(E, :cyclo, ZZ(n))
  return E, b
end

function _is_cyclotomic_field(E::Oscar.Field)
  f = get_attribute(E, :cyclo)
  if f === nothing
    return (false, fmpz(1))
  else
    return (true, f)
  end
end

function trace_lattice(L::Hecke.AbsLat, order::Integer = 2)
  #@req order > 0 "The order must be positive"
  n = degree(L)
  E = base_field(L)
  G = gram_matrix(rational_span(L))

  if E == QQ || E == ZZ
    if order == 1
      f = identity_matrix(E, n)
    elseif order == 2
      f = -identity_matrix(E, n)
    else
      error("For ZLat the order must be 1 or 2")
    end
    return L,f
  end

  bool, m = _is_cyclotomic_field(E)

  #@req bool "Base field must be cyclotomic"

  Eabs, EabstoE = Hecke.absolute_simple_field(E)
  EtoEabs = pseudo_inv(EabstoE)
  d = degree(Eabs)
  Gabs = map_entries(EtoEabs, G)
  z = gen(Eabs)
  chi = minpoly(z)
  Mchi = companion_matrix(chi)
  f = block_diagonal_matrix([Mchi for i=1:n])

  coeffs = fmpq[]
  for i=1:n
    for iz=1:d
      for j=1:n
        for jz=1:d
	  g = z^(iz-jz)*Gabs[i,j]
	  push!(coeffs, trace(g))
	end 
      end
    end
  end
  gram = matrix(QQ, d*n, d*n, coeffs)
  @assert f*gram*transpose(f) == gram
  return Zlattice(gram = 1//d*gram), f
end

function _Sq(q::Integer)
  @assert q > 1
  Sq = [k for k=1:q if gcd(k,q) == 1]
  @assert length(Sq) == euler_phi(q)
  return Sq
end

function _order_max(rank::Integer)
  k = iseven(rank) ? rank : rank+1
  order_max = maximum(euler_phi_inv(k))
  order_max = maximum(filter(n -> divides(ZZ(k), euler_phi(n))[1], 1:order_max))
  return order_max
end

function inverse_trace_lattice(L::ZLat, f::fmpq_mat)
  n = findfirst(n -> is_one(f^n), 1:_order_max(rank(L)))
  #@req ((n !== nothing) || (divides(rank(L), euler_phi(n))[1])) "(L, f) is not the trace lattice of a hermitian lattice over a cyclotomic field"

  E,b = _cyclotomic_tower(n)

  m = divexact(rank(L), euler_phi(n))
  # We choose as basis for the hermitian lattice Lh the identity matrix
  gram = matrix(zeros(E, m,m))
  G = gram_matrix(L)
  Sq = _Sq(n)

  for i=1:m
    for j=1:m
      alpha = sum([G[1+(i-1)*euler_phi(n),:]*transpose(f^k)[:,1+(j-1)*euler_phi(n)] for k in Sq])[1]
      gram[i,j] = alpha
    end
  end

  s = involution(E)
  @assert transpose(map_entries(s, gram)) == gram

  Lh = hermitian_lattice(E, gram = gram)
  return Lh
end

#function _induced_image_in_Oq(qL::Hecke.TorQuadMod, f::fmpq_mat)
#  imgs = elem_type(qL)[]
#  d = denominator(f)
#  m = order(qL)
#  dinv = invmod(d,m)
#  for g in gens(qL)
#    x = lift(g)*f
#    x = QQ(dinv*d)*x
#    push!(imgs, qL(x))
#  end
#  return hom(qL, qL, imgs)
#end

function _eichler_isometry(G, u, v, y, mu)
  n = ncols(G)
  d = absolute_degree(base_ring(G))
  s = involution(base_ring(G))
  @assert u*G*transpose(map_entries(s,u)) == 0
  @assert u*G*transpose(map_entries(s,y)) == 0
  @assert v*G*transpose(map_entries(s,y)) == 0
  V = VectorSpace(base_ring(G) ,n)
  uv = u*G*transpose(map_entries(s,v))[1]
  @assert mu*uv+s(mu*uv) == (-y*G*transpose(map_entries(s,y)))[1]
  imgs = [e + ((e*G*transpose(map_entries(s,u)))[1]//s(uv))*y + ((mu*(e*G*transpose(map_entries(s,u)))[1])//(-s(uv)) - (e*G*transpose(map_entries(s,y)))[1]//uv)*u for e in gens(V)]
  f = matrix(reduce(vcat,imgs))
  @assert G == f*G*transpose(map_entries(s,f))
  return _restrict_scalars(f), one(base_ring(G)), f
end

function _restrict_scalars(M)
  n = ncols(M)
  d = absolute_degree(base_ring(M))
  N = matrix(zeros(QQ, d*n, d*n))
  bas = absolute_basis(base_ring(M))
  for i=0:n-1
    for j=0:n-1
      F[d*i+1:(i+1)*d+1, d*j+1:(j+1)*d+1] = matrix(absolute_coordinates.(M[i,j].*(bas)))
    end
  end
  return F
end

function _skew_element_cyclo(E)
  K = base_field(E)
  b = gen(E)
  s = involution(E)
  M = matrix(K, 2,1,[2, b+s(b)])
  r,k = left_kernel(M)
  @assert r == 1
  c = E(k[1]+b*k[2])
  @assert iszero(c+s(c))
  return c
end

function _symmetry(G, v, sigma)
  s = involution(parent(sigma))
  @assert v*G*transpose(map_entries(s,v)) == sigma+s(sigma)
  n = ncols(G)
  V = VectorSpace(base_ring(G), n)
  imgs = [e - (e*G*transpose(map_entries(s,v)))*inv(sigma)*v for e in gens(V)]
  g = matrix(reduce(vcat, imgs))
  return _restrict_scalars(g), det(g), g
end

function _my_rand_QQ(k::Oscar.IntegerUnion = 64)
  p = ZZ(rand(-k:k))
  q = ZZ(rand(-k:k))
  while q == 0
    q = ZZ(rand(-k:k))
  end
  return p//q
end

function _my_rand(K::AnticNumberField, k::Oscar.IntegerUnion = 64)
  d = degree(K)
  a = gen(K)
  pows = powers(a, d-1)
  z = K(0)
  for zz in pows
    z += zz*my_rand_QQ(k)
  end
  return z
end

function _my_rand(VOE::AbstractAlgebra.Generic.FreeModule{Hecke.NfRelOrdElem{nf_elem, Hecke.NfRelElem{nf_elem}}}, k::Oscar.IntegerUnion = 5)
  n = rank(VOE)
  OE = base_ring(VOE)
  v = vec(zeros(OE, 1, n))
  for i=1:n
    v[i] = rand(OE, k)
  end
  return VOE(v)
end

function _vecQQ_to_VE(VE::VE::AbstractAlgebra.Generic.FreeModule{Hecke.NfRelElem{nf_elem}}, e::Vector{fmpq)
  n = size(e)[2]
  E = base_ring(VE)
  d = absolute_degree(E)
  @assert divides(n,d)[1]
  @assert divexact(n,d) == rank(VE)
  Eabs, EbastoE = absolute_simple_field(E)
  v = vec(zeros(E, 1, divexact(n,d)))
  for i=0:rank(VE)-1
    v[i] = EabstoE(Eabs(e[i*d:(i+1)*d]))
  end
  return VE(v)
end

function image_centralizer_in_Oq(L::ZLat, f::fmpq_mat)
  Lh = inverse_trace_lattice(L, f)
  G = gram_matrix(rational_span(Lh))
  n = degree(Lh)
  @assert n > 1
  E = base_field(Lh)
  invo = involution(E)
  b = gen(E)
  ord = get_attribute(E, :cyclo)
  @req ord >= 3 "f must of order at least 3"
  prime = false

  if length(prime_divisors(ord)) == 1
    prime = true
    OE = maximal_order(E)
    D = different(OE)
    P = prime_decomposition(OE, minimum(D))[1][1]
  end

  qL = discriminant_group(L)
  OqL = orthogonal_group(qL)
  MGf = matrix_group([f])
  fqL = OqL(gens(MGf)[1])
  Oqf, OqftoOqL = centralizer(OqL, fqL)
  K = base_field(E)

  VO = FreeModule(O(E), n)
  VE = VectorSpace(E, n)
  gens = [fqL, OqL(gens(matrix_group([-f^0]))[1])]
  dets = [b^n, E(-1)^n]
  S = subgroup(Oqf, gens)
  w = _skew_element_cyclo(E)
  count = 0
  ind = 1
  flag = true

  M = gram_matrix(lll_reduction(L))
  extra = [M[:,j] for j=1:ncols(M)]
  while order(S) != order(Oqf)
    count +=1
    if flag
      iterate(qL, ind) === nothing && flag = false && @vprint :LWI 1 "Done enumerating the discriminant group\n" && continue
      s, ind = iterate(qL, ind)
    elseif length(extra) != 0
      e = pop!(extra)
      s = _vecQQ_to_VE(VE, e).v
    else
      s = _my_rand(VO).v
    end
    if iszero(s)
      continue
    end
    sigma = (s*G*transpose(map_entries(invo, s)))[1]//2
    @vprint :LWI 1 "Computing spinor norms of Oqd. Total: $(order(Oqf)). Remaining: $(order(Oqf)//order(S)). Number of tries: $(count). \n"
    GC.gc()
    if divides(order(qL), 2)[1] && sigma == 0
      sg = s*G
      (m,i) = minimum([(norm(val), index) for (index, val) in enumerate(sg) if val != 0])
      v = gen(VE, i).v
      dimker, ker = left_kernel((G*transpose(map_entries(invo, vcat(s, v)))))
      for i=1:dimker
        y = ker[i,:]
        mu = (-y*G*transpose(map_entries(invo, y)))[1]//(2*(s*G*transpose(map_entries(invo,v)))[1])
        T, det, TE = _eichler_isometry(G, s, v, y, mu)
        if gcd(denominator(T), order(qL)) == 1
          Tbar = OqL(T)
          Tbar = Oqf(Tbar)
          if Tbar notin S
            push!(gens, Tbar)
            push!(dets, E(1))
            S = sub(Oqf, gens)
          end
        end
      end
    end
    if prime
      I = D*fractional_ideal(OE, s*G)
      if valuation(I, P) <= 0
        continue
      end
    end
    for i in 1:100
      sigma1 = sigma + E(_my_rand(K))*omega
      if sigma1 == 0
        continue
      end 
      tau, determ, tauE = _symmetry(G, s, sigma1)
      if gcd(denominator(tau), order(qL)) != 1
        continue
      end
      taubar = OqL(tau)
      taubar = Oqf(taubar)
      if taubar notin S
        push!(gens, taubar)
        push!(dets, determ)
        S = sub(Oqf, gens)
      end
    end
  end
  Lhn = hermitian_lattice(E, (1//ord)*G)
end
