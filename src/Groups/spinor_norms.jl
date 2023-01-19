export image_in_Oq

@doc Markdown.doc"""
    sigma_sharp(L::ZLat, p) -> Vector{Tuple{fmpz, fmpq}}

Return generators for $\Sigma^\#(L\otimes \ZZ_p)$ of a lattice `L`.

- a list of tuples `(det_p, spin_p)` with `det = +- 1`.
"""
function sigma_sharp(L::ZLat, p)
  T = primary_part(discriminant_group(L), p)[1]
  q = Hecke.gram_matrix_quadratic(normal_form(T)[1])
  res = _sigma_sharp(rank(L), det(L), q, p)
  return [(ZZ(a[1]),QQ(a[2])) for a in res]
end



function _sigma_sharp(rkL, detL, q, p)
  lq = ncols(q)
  delta = detL * det(q)
  up = Hecke._min_nonsquare(p)
  if p != 2
    if rkL == lq
      return [(1, 1)]
    elseif rkL == lq + 1
      return [(-1, 2 * delta)]
    else
      return [(-1, 1), (1, up)]
    end
  end

  # p=2 case
  if rkL > lq + 1
    if p == 2
      return [(-1,1), (1,3), (1,7)]
    end
  end
  blocks = Hecke.collect_small_blocks(q)
  u = [0, 0, 0, 0]
  v = [0, 0, 0, 0]
  w = [[],[],[],[]]
  for b in blocks
    k = valuation(denominator(b), 2)+1
    r = ncols(b)
    if k <= 4
      if r == 2
        if b[1,1] == b[2,2] == 0
          u[k] += 1
        else
          v[k] += 1
        end
      else
        push!(w[k], numerator(b[1,1]))
      end
    end
  end
  gamma20 = [(1,3), (1,7), (-1, 1)]
  gamma21 = [(1,3), (1,7)]
  gamma22 = [(1,5)]

  if rkL > lq
    return gamma20
  end
  if u[2] + v[2] > 0
    if length(w[2]) > 0
      return gamma20
    end
    return gamma21
  end
  if length(w[2]) == 1
    epsilon = numerator(w[2][1])
    delta = detL*ZZ(2)^valuation(-detL, 2) // prod(QQ, [numerator(det(q)) for q in blocks if valuation(denominator(q), 2) > 1])
    if u[3] + v[3] + length(w[4]) > 0
      if length(w[3]) > 0
        return gamma20
      end
      return [(1, 5), (-1, delta)]
    end
    if length(w[3]) == 2
      return gamma20
    end
    if length(w[3]) == 1
      phi = numerator(w[3][1])
      if epsilon*phi % 4 == 3
        return [(1,7), (-1,delta)]
      end
      if epsilon*phi % 4 == 1
        return [(1,3), (-1,delta)]
      end
    end
    return [(-1, delta)]
  end
  if length(w[2]) == 2
    epsilon = numerator(w[2][1])
    phi = numerator(w[2][2])
    if mod(epsilon*phi , 4) == 3
      return gamma20
    end
    if length(w[3]) > 0
      return gamma20
    end
    return [(1,5), (-1,epsilon)]
  end
  @assert length(w[2]) == 0
  if u[3] + v[3] > 0
    return gamma22
  end
  if length(w[3]) == 2
    return gamma22
  end
  if length(w[3])<= 2
    return []
  end
  @assert false
end

@doc Markdown.doc"""
    reflection(gram::fmpq_mat, v::fmpq_mat) -> fmpq_mat

Return the matrix representation of the orthogonal reflection in the row vector `v`.
"""
function reflection(gram::MatElem, v::MatElem)
  n = ncols(gram)
  E = identity_matrix(base_ring(gram), n)
  c = base_ring(gram)(2) * ((v * gram * transpose(v)))[1,1]^(-1)
  ref = zero_matrix(base_ring(gram), n, n)
  for k in 1:n
    ref[k,:] = E[k,:] - c*(E[k,:] * gram * transpose(v))*v
  end
  return ref
end

@doc Markdown.doc"""
    spin(gram_diag::MatElem, isometry::MatElem, check=true) -> fmpq

Compute the spinor norm of `f`.

`gram_diag` must be a diagonal matrix.
"""
function spin(gram_diag::MatElem, isometry::MatElem, check=true)
  G = gram_diag
  f = isometry
  @assert ncols(G) == nrows(G) == ncols(f) == nrows(f) "G and f must be square matrices"
  @assert is_diagonal(G)
  if check
    @assert G == f * G * transpose(f)  "f must be an isometry"
  end

  n = ncols(G)
  R = base_ring(G)
  spinor_norm = one(R)
  for i in 1:n
    w = zero_matrix(R, 1, n); w[1,i] = 1
    v = w * f
    r = v - w
    s = r * G * transpose(r)
    if !iszero(s)
      tau = reflection(G, r)
      f = f * tau
      @assert w * f == w
      spinor_norm *= s
    else
      r1 = v + w
      s1 = r1 * G * transpose(r1)/2
      r2 = v
      s2 = r2 * G * transpose(r2)/2
      @assert !iszero(s1) && !iszero(s2)
      tau1 = reflection(G, r1)
      tau2 = reflection(G, r2)
      f = f * tau2 * tau1
      @assert w * f == w
      spinor_norm *= s1 * s2
    end
  end
  @assert isone(f)
  return spinor_norm[1,1]
end

@doc Markdown.doc"""
    det_spin(G::fmpq_mat, T::fmpq_mat, p, nu) -> Tuple{fmpq, fmpq}

Return approximations for `(det_p, spin_p)` of the approximate isometry `T`.

The algorithm is by Shimada [Shim2018](@cite)

We follow the conventions of Miranda and Morrison that the quadratic form is defined by
`Q(x) = (x G x.T)/2`. Then the spinor norm of the reflection in x is Q(x).

# Arguments
- `G::fmpq_mat`: a diagonal matrix
- `T::fmpq_mat`: an isometry up to some padic precision
- `p`: a prime number
- `nu`: an integer giving the valuation of the approximation error of `T`
"""
function det_spin(G::fmpq_mat, T::fmpq_mat, p, nu)
  p = ZZ(p)
  if p == 2
    delta = 1
  else
    delta = 0
  end
  gammaL = [valuation(d, p) for d in diagonal(G)]
  gamma = minimum(gammaL)
  l = ncols(G)
  E = parent(G)(1)
  reflection_vectors = []

  k = 1
  while k <= l
    g = T[k,:]
    # error estimates
    lambd = valuation(g, p)
    rho = min(delta + nu + gamma, 2*nu + gamma)
    sigma = min( delta + nu + gamma, delta + nu + lambd, 2*nu + gamma)
    kappa = sigma - gammaL[k] - 2*delta
    if (rho <= gammaL[k] + delta) || (kappa < 1 + 2*delta)
      # precision too low
      return ZZ(0), QQ(0)
    end
    bm = g - E[k,:]
    qm = bm * G * transpose(bm)
    if valuation(qm, p) <= gammaL[k] + 2*delta
      tau1 = reflection(G, bm)
      push!(reflection_vectors, bm)
      tau2 = E
    else
      bp = g + E[k,:]
      qp = bp * G * transpose(bp)
      @assert valuation(qp, p) <= gammaL[k] + 2*delta
      tau1 = reflection(G, bp)
      tau2 = reflection(G, E[k,:])
      push!(reflection_vectors,bp)
      push!(reflection_vectors,E[k,:])
    end
    lambdaT = valuation(T, p)
    alpha = valuation(tau1, p)
    beta = valuation(tau2, p)
    theta = gamma + min(kappa + 2*min(0, lambd), nu + min(0, lambd), 2*nu)
    nu = min(nu + alpha, lambdaT + theta - gammaL[k] - delta, nu + theta - gammaL[k] - delta) + beta
    T = T * tau1 * tau2
    k += 1
  end
  err = valuation(T - E, p)
  @assert err >= nu
  spinor_norm = prod(QQ, [(v*G*transpose(v))[1,1]*(1//2) for v in reflection_vectors])
  determinant = QQ(-1)^(length(reflection_vectors))
  v = valuation(spinor_norm, p)
  u = divexact(spinor_norm, p^v)
  if p == 2
      u = mod(u, 8)
  else
      u = mod(u, p)
  end
  spinor_norm = u * p^mod(v, 2)
  return determinant, spinor_norm
end

#    Elements of the product over `C_2` x `\QQ_p* / (\QQ_p*)^2` at primes `p`.
function _det_spin_group(primes::Vector{fmpz}; infinity = true)
  #@assert infinity
  K, _ = Hecke.rationals_as_number_field()
  # f : QQ -> K
  f = MapFromFunc(x -> K(x), x -> coeff(x, 0), QQ, K)
  OK = maximal_order(K)
  primes_as_ideals = [prime_decomposition(OK, p)[1][1] for p in primes]
  stuff = [Hecke.local_multiplicative_group_modulo_squares(P) for P in primes_as_ideals]
  grps = [s[1] for s in stuff]
  maps = Any[s[2] for s in stuff]
  if infinity
    Ainf = abelian_group(2)
    minf = MapFromFunc(x -> iszero(x[1]) ? one(K) : -one(K), x -> coeff(x, 0) > 0 ? Ainf([0]) : Ainf([1]), Ainf, K)
    push!(grps, Ainf)
    push!(maps, minf)
  end
  A, proj, inj = direct_product(grps..., task = :both)
  backwardmap = x -> sum([inj[i](maps[i]\(f(x))) for i in 1:length(maps)])
  forwardmap = function(x)
    elems = [f\(maps[i](proj[i](x))) for i in 1:length(grps)]
    elems_integral = fmpz[]
    for i in 1:(length(elems) - 1)
      push!(elems_integral, ZZ(denominator(elems[i])^2 * elems[i]))
    end
    cprimes = copy(primes)
    for i in 1:length(primes)
      if cprimes[i] == 2
        cprimes[i] = primes[i]^4
      else
        cprimes[i] = primes[i]^3
      end
    end
    y = crt(elems_integral, cprimes)
    if sign(y) == sign(elems[end])
      z = QQ(y)
    else
      z = QQ(y + sign(elems[end]) * prod(cprimes))
    end
    @assert backwardmap(z) == x
    return z
  end
  grps_det = [abelian_group(2) for i in 1:length(primes)]
  push!(grps_det, A)
  D, projD, injD = direct_product(grps_det...,task=:both)
  maps_det = [(primes[i],MapFromFunc(x-> isone(x) ? zero(grps_det[i]) : grps_det[i][1], ZZ, grps_det[i])*injD[i]) for i in 1:length(primes)]
  maps_det = Dict(maps_det)
  projd = Any[(primes[i],projD[end]*proj[i]*maps[i]*inv(f)) for i in 1:length(primes)]
  injd = Any[(primes[i],f*inv(maps[i])*inj[i]*injD[end]) for i in 1:length(primes)]
  if infinity
    push!(projd,(PosInf(), proj[end]))
    push!(injd,(PosInf(), inj[end]))
  end
  projd = Dict(projd)
  injd = Dict(injd)
  diagonal_morphism = MapFromFunc(forwardmap, backwardmap, A, QQ)
  return D, inv(diagonal_morphism)*injD[end], projd, injd, maps_det
end

@doc Markdown.doc"""
    det_spin_homomorphism(L::ZLat) -> GAPGroupHomomorphism

Return the det spin homomorphism.
"""
function det_spin_homomorphism(L::ZLat; signed=false)
  T = discriminant_group(L)
  Oq = orthogonal_group(T)
  S = prime_divisors(2 * order(domain(Oq)))
  A, diagonal, proj,inj,det_hom = _det_spin_group(S, infinity=false)

  # \Sigma^\#(L)
  # This is the the image of K under `(det, spin)` where `K` is the kernel
  # of O(L) --> O(L^v / L)
  sigma_sharp_gens = elem_type(A)[]
  for p in S
    ss = sigma_sharp(L, p)
    for s in ss
      push!(sigma_sharp_gens, det_hom[p](s[1])+inj[p](s[2]))
    end
  end

  result = Dict([(f, zero(A)) for f in gens(Oq)])


  #=
  \Gamma_S = \{(d,s) \in \Gamma_\QQ \mid s \in S \}

    MATH:

    \Gamma_S^+ = ker(\Gamma_\QQ \to \{\pm 1\}, (d,s) \mapsto \sign(ds))
  =#
  GammaS = [diagonal(QQ(s)) for s in S]
  if signed
    push!(GammaS, diagonal(QQ(-1))+sum([det_hom[p](ZZ(-1)) for p in S]))
  else
    push!(GammaS, sum([det_hom[p](ZZ(-1)) for p in S]))
    push!(GammaS, diagonal(QQ(-1)))
  end
  gens_ker = append!(sigma_sharp_gens, GammaS)
  _, j = sub(A, gens_ker)
  D, proj = cokernel(j)

  for p in S
    if p == -1
      continue
    end
    Tp,iTp = primary_part(T, p)
    Tpn,iTpn = normal_form(Tp)
    iT = inv(iTpn)*iTp
    if order(Tp) == 1
      continue
    end
    u = reduce(vcat,[b.data.coeff for b in gens(Tp)])
    Op = orthogonal_group(Tpn)
    q = gram_matrix_quadratic(Tpn)
    # Shimada, Connected Components of the Moduli of Elliptic K3 Surfaces
    # 5.2 Step 2 page 539
    # The discriminant group has full length and a term 1/2 || 3/2
    # to determine the lattice one needs to know the determinant
    tmp = (det(q)*det(L))
    val = valuation(tmp, p)
    unit = divexact(tmp, p^val)
    @assert val == 0
    if p == 2 && length(elementary_divisors(Tp))== rank && unit % 8 != 1
      n = ncols(q)
      for i in 2:n
        if q[i, i+1]==0 && q[i - 1, i]==0 && denominator(q[i, i])==2
          q[i, i] *= 5
          @goto breaked
        end
      end
      # no break
      if denominator(q[1,1])==2 && (n==1 || q[1,2]==0)
        q[1,1] *= 5
      elseif denominator(q[end,end])==2 && q[end,end-1]==0
        q[end,end] *= 5
      else
        @assert false "bug in det_spin_homomorphism"
      end
      @label breaked
    end
    M = inv(q)
    # diagonalize
    diag, t = Hecke._gram_schmidt(M,QQ)
    # t = transpose(t) ???
    q0 = q*denominator(q)
    v1 = valuation(denominator(t), p)
    v2 = valuation(denominator(inv(t)), p)
    v = -v1 -v2 # lower bound for precision loss due to diagonalization
    # compute the spin of lifts of the generators
    for f in gens(Oq)
      # take only the action on the p-part
      fp = hom(Tpn,Tpn,[preimage(iT,f(iT(g))) for g in gens(Tpn)])
      fp = matrix(Op(fp))
      prec0 = 1
      prec = 25 # initial precision
      # change to the user basis
      g = u * fp * inv(u)
      while true
        R = ResidueRing(ZZ, p^(prec+3))
        conv = MapFromFunc(x -> R(numerator(x)) * R(denominator(x)^(-1)), QQ, R)
        _g = Hecke.hensel_qf(map_entries(conv, q0), change_base_ring(R, g), prec0, prec, p)
        g = change_base_ring(ZZ, _g)
        gg = t*M*g*inv(t*M)
        det_p, spin_p = det_spin(diag, gg, p, prec + v)
        if det_p != 0
          result[f]+= inj[p](QQ(spin_p)) + det_hom[p](ZZ(det_p))
          break
        end
        prec0 = prec
        prec = 2 * prec
      end
    end
  end
  Agap, i1,_ = _isomorphic_gap_group(D)
  return hom(Oq,Agap,gens(Oq),[i1(proj(result[f])) for f in gens(Oq)],check=false)
end


@doc Markdown.doc"""
    image_in_Oq(L::ZLat) -> AutomorphismGroup{Hecke.TorQuadMod}, GAPGroupHomomorphism

Return the image of $O(L) \to O(L^\vee / L)$.

```
julia> L = 2*root_lattice(:D,4)
Quadratic lattice of rank 4 and degree 4 over the rationals

julia> imOL,inj = image_in_Oq(L);

julia> order(imOL)
1152

```

By using the strong approximation theorem, we can also deal with
indefinite lattices.

```
julia> gram = ZZ[0 2 0; 2 0 0; 0 0 4]
[0   2   0]
[2   0   0]
[0   0   4]

julia> L = Zlattice(gram=gram)
Quadratic lattice of rank 3 and degree 3 over the rationals

julia> Oq, inj = image_in_Oq(L);

julia> order(Oq)
12

```
"""
@attr function image_in_Oq(L::ZLat)::Tuple{AutomorphismGroup{Hecke.TorQuadMod}, GAPGroupHomomorphism{AutomorphismGroup{Hecke.TorQuadMod}, AutomorphismGroup{Hecke.TorQuadMod}}}
  @req iseven(L) "Implemented only for even lattices so far. If you really need this, you can rescale the lattice to make it even and then project the orthogonal group down."
  if rank(L) > 2 && !is_definite(L)
    # use strong approximation
    f = det_spin_homomorphism(L,signed=false)
    return kernel(f)
  end
  # we can compute the orthogonal group of L
  Oq = orthogonal_group(discriminant_group(L))
  G = orthogonal_group(L)
  return sub(Oq, [Oq(g, check=false) for g in gens(G)])
end

@attr function image_in_Oq_signed(L::ZLat)::Tuple{AutomorphismGroup{Hecke.TorQuadMod}, GAPGroupHomomorphism{AutomorphismGroup{Hecke.TorQuadMod}, AutomorphismGroup{Hecke.TorQuadMod}}}
  @req iseven(L) "Implemented only for even lattices so far. If you really need this, you can rescale the lattice to make it even and then project the orthogonal group down."
  @req rank(L) > 2 && !is_definite(L) "L must be indefinite of rank at least 3"
  # use strong approximation
  f = det_spin_homomorphism(L,signed=true)
  return kernel(f)
end
