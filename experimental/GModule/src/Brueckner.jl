module RepPc
using Oscar

#=TODO
 - construct characters along the way as well?
 - compare characters rather than the hom_base
 - maybe reason from theory what reps are going to be new?
 - conjugate to smallest field?
 - allow trivial stuff
=#
"""
  For K a finite field, Q, a number field or QQAb, find all
abs. irred. representations of G.

Note: the reps are NOT necessarily over the smallest field.

Note: the field is NOT extended - but it throws an error if it was too small.

Note: `group(M)` for the returned gmodules `M` will have a pcgs of `G` as
      its `gens` value, thus these generators will in general differ from
      the generators of `G`.

Implements: Brueckner, Chap 1.2.3
"""
function reps(K, G::Oscar.PcGroup)
  @req is_finite(G) "the group is not finite"
  if order(G) == 1
    F = free_module(K, 1)
    h = hom(F, F, [F[1]])
    return [gmodule(F, G, typeof(h)[h for i = gens(G)])]
  end

  pcgs = Oscar.GAPWrap.Pcgs(GapObj(G))
  @assert length(pcgs) == ngens(G)
  pcgs == GAP.Globals.fail && error("the group is not polycyclic")

  gG = [Oscar.group_element(G, x) for x = pcgs]
  s, ms = sub(G, [gG[end]])
  o = Int(order(s))
  @assert is_prime(o)
  z = roots(K(1), o)
  @assert characteristic(K) == o || length(z) == o
  F = free_module(K, 1)
  R = [gmodule(F, s, [hom(F, F, [r*F[1]])]) for r = z]
  @hassert :BruecknerSQ 2 Oscar.GrpCoh.is_consistent(R[1])

  for i=length(gG)-1:-1:1
    h = gG[i]
    ns, mns = sub(G, gG[i:end])
    @assert mns(ns[1]) == h
    p = Int(divexact(order(ns), order(s)))
    @assert is_prime(p)
    new_R = []
    todo = trues(length(R)) # which entries in `R` have to be handled
    #TODO: use extend below
    for pos in 1:length(R)
      if todo[pos]
        r = R[pos]
        F = r.M
        @assert group(r) == s
        rh = gmodule(group(r), [action(r, preimage(ms, ms(x)^h)) for x = gens(s)])
        @hassert :BruecknerSQ 2 Oscar.GrpCoh.is_consistent(rh)
        l = Oscar.GModuleFromGap.hom_base(r, rh)
        @assert length(l) <= 1
        Y = matrix(action(r, preimage(ms, h^p)))
        if length(l) == 1
          # The representation extends from the subgroup,
          # all these extensions are pairwise inequivalent.
          X = l[1]
          Xp = X^p
          #Brueckner: C*Xp == Y for some scalar C
          ii = findfirst(!is_zero, Xp)
          @assert !iszero(Y[ii])
          C = divexact(Y[ii], Xp[ii])
          @assert C*Xp == Y
          # I think they should always be roots of one here.
          # They should - but they are not:
          # Given that X is defined up-to-scalars only, at best
          # C is a root-of-1 * a p-th power:
          # Y is in the image of the rep (action matrix), hence has
          # finite order (at least if the group is finite), hence
          # det(Y) is a root-of-1, so X is defined up to scalars,
          # xX for x in the field., hence Xp = X^p is defined up
          # to p-th powers: x^p Xp, so
          # C x^p Xp = Y
          # applying det:
          # det(C x^p Xp) = C^n x^(pn) det(Xp) = det(Y) = root-of-1
          # so I think that shows that C is (up to p-th powers)
          # also a root-of-1
          #
          # However, I don't know how to use this...
          rt = roots(C, p)
          @assert characteristic(K) == p || length(rt) == p
          Y = r.ac
          for x = rt
            nw = gmodule(F, ns,  vcat([hom(F, F, x*X)], Y))
            @hassert :BruecknerSQ 2 Oscar.GrpCoh.is_consistent(nw)
            push!(new_R, nw)
          end
        else #need to extend dim
          n = dim(r)
          F = free_module(K, dim(r)*p)

          # a block permutation matrix for the element `h`
          z = zero_matrix(K, dim(F), dim(F))
          z[1:n,(p-1)*n+1:end] = Y
          #= This is wrong in Brueckner - or he's using a different
             conjugation. Max figured out what to do: the identity block
             needs to be lower left, and upper right the inverse.

             He might have been doing other conjugations s.w.
          =#
          for ii=2:p
            z[(ii-1)*n+1:ii*n, (ii-2)*n+1:(ii-1)*n] = identity_matrix(K, n)
          end
          md = [hom(F, F, z)]

          conjreps = [eltype(md)[] for i in 1:p]
          M = free_module(K, dim(r))

          # a block diagonal matrix for each generators of `s`
          for g = gens(s)
            z = zero_matrix(K, dim(F), dim(F))
            for j=1:p
              Y = action(r, g)
              m = matrix(Y)
              z[(j-1)*n+1:j*n, (j-1)*n+1:j*n] = m
              push!(conjreps[j], hom(M, M, m))
              g = preimage(ms, ms(g)^h)
            end
            push!(md, hom(F, F, z))
          end

          # Find the positions of the equiv. classes of the `h`-conjugate
          # representations.
          # In order to get pairwise equivalent representations,
          # we do not induce the representations equivalent to the
          # other conjugates.
          for j in 2:p
            Gj = gmodule(M, s, conjreps[j])
            for k in (pos+1):length(R)
              if is_isomorphic(Gj, R[k]) > 0
                todo[k] = false
              end
            end
          end

          push!(new_R, gmodule(F, ns, md))
          @hassert :BruecknerSQ 2 Oscar.GrpCoh.is_consistent(new_R[end])
        end
      end
    end
    s, ms = ns, mns
    R = new_R
  end
  return [gmodule(x.M, G, x.ac) for x = R]
end


"""
Brueckner Chap 1.3.1

Given
  mp: G ->> Q

Find a set of primes such that if there is any irreducible F_p module M
s.th. there is an epimorphism of G onto the extension of Q by M,
the p is in the set.

If `easy` is set to false, then only the primes coming from the maximal abelian quotient are returned. This is faster - but not alwasy enough.
"""
function find_primes(mp::Map{<:Oscar.GAPGroup, PcGroup}; easy::Bool = false)
  if easy
    inv = abelian_invariants(kernel(mp)[1])
    #TODO: is this a good idea? kernel might be expensive...
    @req !(0 in inv) "the kernel has infinite abelianization, so every prime is admissible; pass `primes`"
    return sort!(unique!(reduce(vcat, (prime_divisors(x) for x in inv), init = ZZRingElem[])))
  end
  G = domain(mp)
  Q = codomain(mp)
  if order(Q) == 1
    F = free_module(ZZ, 1)
    I = [gmodule(F, Q, [hom(F, F, [F[1]]) for x in gens(Q)])]
  else #TODO: repsn, reps offer choice
    if is_abelian(Q) #much faster, new Breuer method
      I = irreducible_modules(ZZ, Q)
    else
      I = reps(abelian_closure(QQ)[1], Q)
      #Brueckner, p35: irreducible here is not necessary, so the
      #  expensive find minimal field step can be omitted.
      #also: we need only one per Galois-orbit as conjugate ones will
      #coincide after restriction
      I = [gmodule(ZZ, gmodule(QQ, gmodule(CyclotomicField, x))) for x = I]
    end
  end
  lp = Set(prime_divisors(order(Q)))
  ZG = free_res(group_algebra(ZZ, G); side = :right, one_only = true)
  for i = I
    ib = gmodule(i.M, G, [action(i, mp(g)) for g = gens(G)])
    ia = gmodule(FinGenAbGroup, ib)
    X = hom(ZG, ia)
    a = map(X, 0)
    b = map(X, 1)
#    a, b = H_one_maps(ia)    
#    da = Oscar.dual(a)
#    db = Oscar.dual(b)
    #=
    R = Q/Z, then we should have
      R^l -a-> R^n -b-> R^m
    and the H^1 we want is ker(b)/im(a)
    however, actually, a and b run between Z^l's.
    Taking duals:
      Z^l <-a'- Z^n <-b'- Z^m
    should give me
      quo(im(b'), ker(a'))
    as the dual to what I want.
    ======================================
    Wrong / don't know why correct.
    2nd attempt:
    Im(a) = R^? as R is divisible, the image is a quotient (domain modulo
    kernel), hence a power of R
    Ker(b) = R^? x Torsion:
    transform b in SNF (change of basis R^n and R^m with Gl(n, Z))
    then Ker = R^(number of 0) x T and the T is the non-zero elem. divisors.
    Thus the quotient is R^? x T
    (no duality was harmed here)
    (the cohomology also works as the action (matrices) are identical for
    R and Z (namely integral) and the cohomology does not do any computation
    right until the end when images and kernels are obtained. The Maps
    are correct...)
    TODO: this is not (yet) implemented this way
    =#
    sb = snf(matrix(b))
    sa = snf(matrix(a))
    ra = length([i for i=1:nrows(sa) if !is_zero_row(sa, i)])
    rb = length([i for i=1:nrows(sb) if is_zero_row(sb, i)])
    @req ra == rb "group must be finite"
    t = prod([sb[i,i] for i=1:min(nrows(sb), ncols(sb)) if !is_zero_entry(sb, i, i)]; init = ZZ(1))
    #careful: if the rank of the image is < than the rank of the (free bit)
    #of the kernel, then there are infinite extension possible, i.e. all primes
    #Lemma 1, p27
    if t > 1
      push!(lp, prime_divisors(t)...)
    end
  end
  return lp
end

#TODO: redo, sensible strategy:
# - do not search for primes always, only if after some extension the process
#   stopped
# - possibly extend the reps by inflation to see if we actually need new ones
# - if we do this, then don't compute/ use reps that are the result of
#   inflation
"""
Given
    mQ: G ->> Q
Find all possible extensions of Q by an irreducible F_p module
that admit an epimorphism from G.
Implements the SQ-Algorithm by Brueckner, Chap 1.3

If necessary, the prime(s) p that can be used are computed as well.
"""
function brueckner(mQ::Map{<:Oscar.GAPGroup, PcGroup}; primes::Vector=[], limit::Int = typemax(Int), easy::Bool = false)
  Q = codomain(mQ)
  G = domain(mQ)
  @vprint :BruecknerSQ 1 "lifting $mQ using SQ\n"
  if length(primes) == 0
    @vprint :BruecknerSQ 1 "primes not provided, searching...\n"
    lp = find_primes(mQ; easy) 
  else
    lp = map(ZZRingElem, primes)
  end
  @vprint :BruecknerSQ 1 "using primes $lp\n"

  allR = []
  for p = lp
    _, j = ppio(exponent(Q), p)
    f = j == 1 ? 1 : modord(p, j)
    @assert (p^f-1) % j == 0
    @vprint :BruecknerSQ 2 "computing reps over GF($p, $f)\n"
    if f == 1
      @vtime :BruecknerSQ 2 I = reps(GF(Int(p)), Q)
    else
      @vtime :BruecknerSQ 2 I = reps(GF(Int(p), f), Q)
    end
    @vprint :BruecknerSQ 1 "have $(length(I)) representations\n"

    for i = I
      @vprint :BruecknerSQ 1 "starting to process module\n"
      @vprint :BruecknerSQ 2 "... transfer over min. field\n"
      @vtime :BruecknerSQ 2 ii = Oscar.GModuleFromGap.gmodule_minimal_field(i)
      @vprint :BruecknerSQ 2 "... lift...\n"
      #TODO: why do we need the module over GF(p)??? 
      #we don't
      #iii = Oscar.GModuleFromGap.gmodule(GF(Int(p)), ii)
      @vtime :BruecknerSQ 2 l = lift(ii, mQ; limit = limit - length(allR))
      @vprint :BruecknerSQ 2 "found $(length(l)) many\n"
      #TODO: in Plesken p119 has more comments what not to do
      append!(allR, [x for x in l])# if is_surjective(x)])
      if length(allR) >= limit
        return allR
      end
    end
  end
  return allR
end

function trivial_chain(C::GModule, n::Int)
  #TODO: do for other n as well...(or change the name)
  @assert n == 2
  G = C.G
  c = Dict((one(G), one(G)) => zero(C.M))
  S = elem_type(C.G)
  T = elem_type(C.M)
  return Oscar.GrpCoh.CoChain{2, S, T}(C, c, x->zero(C.M))
end

lift_list = []
"""
  mp: G ->> Q
  C a F_p[Q]-module
  Find all extensions of Q my C s.th. mp can be lifted to an epi.
"""
function lift(C::GModule, mp::Map; limit::Int = typemax(Int))
#  global lift_list
#  push!(lift_list, (C, mp))
  #m: G->group(C)
  #compute all(?) of H^2 that will describe groups s.th. m can be lifted to

  G = domain(mp)
  N = group(C)
  @assert isa(N, PcGroup)
  @assert codomain(mp) == N
  ZN = free_res(group_algebra(ZZ, N); side = :right)

  R = relators(G)
  M = C.M

  #=
    G    -->> N
    |
    V  this is needed
    V
    H    -->> N for the new group

   thus G ni g -> (n, m) for n in N and m in the module.
   for this to work, the relations in G need to be satisfied for the images


   g_i is mapped to (m(g_i), pro[i](D))
   this needs to be "collected"
  =#

  #@assert ngens(G) > 0 && length(R) > 0
  if ngens(G) == 0
    D = sub(M, [zero(M)])[1]
    pro = [zero_map(D, M)]
    inj = [zero_map(M, D)]
  else
    D, pro, inj = direct_product([M for i=1:ngens(G)]..., task = :both)
  end
  #hopefully, if ngens(G) == 0, then the group is trivial and hopefully
  #we're not getting here...
  if length(R) == 0
    K = free_module(base_ring(M), 0)
    pK = [zero_map(K, M)]
    iK = [zero_map(M, K)]
  else
    K, pK, iK = direct_product([M for i=1:length(R)]..., task = :both)
  end
  #if ngens(G) > 0 and length(R) == 0 the group is infinite
  #and we should not get here...
  S = relators(N)
  # length(S) == 0 can happen in the very 1st step: we might start
  # by G ->> trivial pc group
  if length(S) != 0
    X, pX, iX = direct_product([M for i=1:length(S)]..., task = :both)
  end

  #=following Max/ Claude:
   the lifting assignes module elems to generators of G
   in each legal extension, those have to map to zero under the relations

   For a specific group extension given by mu in H^2(N, C)
   this is a map C^#G -> C^#R (#G no gens, #R no rels of G)
   for mu zero, this map is linear
   In general we need the kernel of the map - for each mu in H^2
   Max/Claude decompose this, the map is affine: linear + rhs
   The linear part is independent of mu, depending only on C
   (realised via mu = 0)
   The rhs is linear from H^2(N, C) -> C^#R
   So this is done by evaluating at generators of H^2.
   This drops the evaluations from |H^2| -> dim H^2, so logarithmic
  =#
  #start with the map
  GG = split_extension(C)
  mGG_PC = isomorphism(PcGroup, GG)

  im = elem_type(K)[]
  gns = [GG(mp(G[i]), zero(M)) for i=1:ngens(G)]
  for h = gens(D)
    gns = [GG(mp(G[i]), pro[i](h)) for i=1:ngens(G)]
    push!(im, K([map_word(R[i], gns).m for i=1:length(R)]))
  end
  s = hom(D, K, im)

  #the H^2-> RHS is not always used, so we defer this slightly.


  function _process_trivial(mu; is_trivial::Bool = false, limit::Int)
    @assert is_trivial
    local allG = typeof(mp)[]

    k, mk = kernel(s)
    if order(k) == 1
      @show :nothing, 1
      @assert length(allG) == 0
      return allG
    end

    h1 = cohomology_group(hom(ZN, C), 1; no_image = true) #need the no-image tag..
    #want a map from k1[1] (the kern) to D (images of the gens of the
    #large group)
    #let l be a 1-chain in k1[1]
    #gen(G) -mp-> N and we map on the 1-chain evaluated there
    #note: we need Z^1, no H^1, hence the no_image to not quotient
    #out. Needs different interface
    _imgs = elem_type(D)[]
    for i = gens(h1[1])
      c1 = h1[3](i)
      push!(_imgs, sum(inj[l](c1((mp(G[l]),))) for l=1:ngens(G)))
    end
    inf = hom(h1[1], D, _imgs)
    q, mq = quo(k, image(inf)[1])
    k = q #according to Brueckner lifting with multiplicity is
    #linked to ngens(k) at this point. S is the kernel, S' the image of inf
    #thus the number a is ngens(k) - ngens(image) == ngens(q)...
    mk = pseudo_inv(mq)*mk
    if order(k) == 1
      @show :nothing, 2
      @assert length(allG) == 0
      return allG
    end
    @show :split, k
    if ngens(k) > 1
      @show :lifting_with_mult, ngens(k)
      CC, injCC = direct_product([C for i=1:ngens(k)]...; task = :sum)
      muCC = Oscar.GrpCoh.CoChain{2, PcGroupElem, elem_type(CC.M)}(CC,
        Dict{NTuple{2, PcGroupElem}, elem_type(CC.M)}(),
        x->sum(injCC[i](mu.D(x)) for i=1:ngens(k)))
      GG = extension(muCC)    
      mGG_PC = isomorphism(PcGroup, GG)
    else
      injCC = [identity_map(C.M)]
    end
    im = [GG(mp(G[i]), sum(injCC[l](pro[i](mk(k[l]))) for l=1:ngens(k))) for i=1:ngens(G)]
    #@show [map_word(R[i], im) for i=1:length(R)]
    im = map(mGG_PC, im)
    #@show [map_word(R[i], im) for i=1:length(R)]
    hm = hom(G, codomain(mGG_PC), im)
    push!(allG, hm)
    @assert is_surjective(hm)
    return allG
  end

  mu = trivial_chain(C, 2)
  allG = _process_trivial(mu; is_trivial = true, limit)
  if length(allG) >= limit || gcd(order(C.G), order(C.M)) == 1 #trivial H^2
    return allG
  end

  H2, _, z = cohomology_group(hom(ZN, C), 2)
  @show H2

  im = elem_type(K)[]
  for g = gens(H2)
    GG = extension(z(g))
    gns = [GG(mp(g), zero(C.M)) for g = gens(G)]
    push!(im, K([map_word(R[i], gns).m for i=1:length(R)]))
  end
  mH2_K = hom(H2, K, im)

  #useful is only the part of H2 that lands in the image...
  #I'd like preimage(mH2_K, image(s)[1])
  is, mis = image(s)
  ih, mih = image(mH2_K)
  b, mb = intersect(is, ih) #mb: b -> is
  k, mk = kernel(mH2_K)
  H2, mH2 = sub(H2, vcat([preimage(is, mb(g)) for g = gens(b)], map(mk, gens(k))))
  @show H2
  z = mH2*z
  mH2_K = mH2*mH2_K
 
 #= 
  ImH2 = image(mH2_K)
  Is, mIs = image(s)
  H2, mH2 = intersect(ImH2[1], Is)
  @show H2
  z = mH2*ImH2[2]*pseudo_inv(mH2_K)*z
  mH2_K = mH2*ImH2[2]
  @show z
 =#
 

  function _process(h; is_trivial::Bool = false, limit::Int)
    @assert !is_trivial
    @show h
    rhs = mH2_K(h)
    mu = z(h)

    local allG = typeof(mp)[]

    im = elem_type(K)[]
    fl, pe = has_preimage_with_preimage(s, rhs)
    @assert fl
    if !fl
      return allG
    end
    @show rhs
    k, mk = kernel(s)
    @show :non_split, k
    GG = extension(mu)
    mGG_PC = isomorphism(PcGroup, GG)
    
    for x = k
      im = [(GG(mp(G[i]), pro[i](mk(x)-pe))) for i=1:ngens(G)]
      #@show [map_word(R[i], im) for i=1:length(R)]
      im = map(mGG_PC, im)
      #@show [map_word(R[i], im) for i=1:length(R)]
      hm = hom(G, codomain(mGG_PC), im)
      push!(allG, hm)
      if length(allG) >= limit
        return allG
      end
      @assert is_surjective(hm)
    end
    return allG
  end


  #XXX: if done correctly, C is abs. irr over some Fq
  #     the endo is Fq and we need for each "Fq line" only one rep
  #     possibly the one where the 1st non-zero entry is 1
  #     cohomology elements are in C^k as well - or the representing chains
  #     are
#  @show Ends = Oscar.GModuleFromGap.hom_base(C, C)
#  @assert length(Ends) == degree(base_ring(C))

  seen = [zero(H2)]
  for h = H2
    h in seen && continue
    for x = base_ring(H2) #the line will not give "different" quotients
      iszero(x) && continue
      push!(seen, x*h)
    end

    append!(allG, _process(h; is_trivial = false, limit = limit - length(allG)))
    if length(allG) >= limit
      return allG
    end
  end

  return allG
end

function solvable_quotient(G::Oscar.GAPGroup)
  A, _ = maximal_abelian_quotient(G)
  if is_finite(A)
    # not `maximal_abelian_quotient(PcGroup, G)`: GAP hands back a pc group on
    # a non-canonical pcgs for some inputs, and `isomorphism(PcGroup, .)`
    # refuses those. Going through `FinGenAbGroup` always gives a full pc
    # group, which `reps` needs.
    B, mB = maximal_abelian_quotient(FinGenAbGroup, G)
    iso = isomorphism(PcGroup, B)
    return hom(G, codomain(iso), [iso(mB(x)) for x in gens(G)])
  else
#    error("group must be finite - this one has an infinite maximal abelian quotient")
  end
  #the old version: start deliberately with the trivial quotient
  q = cyclic_group(1)
  mp = hom(G, q, [one(q) for g in gens(G)])
end

function sq(mp::Map, primes::Vector=[]; index::Union{Integer, ZZRingElem, Nothing} = nothing)
  if index !== nothing
    lf = factor(ZZRingElem(index))
    primes = prime_divisors(ZZ(index))
    while length(primes) > 0
      @vtime :BruecknerSQ 2 nw = brueckner(mp; limit = 1, primes)
      if length(nw) == 0 
        return mp
      end
      mp = nw[1]
      for (p, k) = lf
        if p in primes && valuation(order(codomain(mp)), p) >= k
          deleteat!(primes, findfirst(isequal(p), primes))
#          @show :removing, p
        end
      end
    end
    return mp
  end
  while true
    nw = brueckner(mp; limit = 1, primes, easy = true)
    #easy == true means use primes from the max. ab. quotient only
    #as this is not enough, we need to try harder
    #TODO: fix the logic so that ONLY the new primes are actually tested
    if length(nw) == 0 && length(primes) == 0
      primes = setdiff(find_primes(mp; easy = false), find_primes(mp; easy = true))
      nw = brueckner(mp; limit = 1, easy = false)
    end
    if length(nw) == 0
      return mp
    end
    mp = nw[1]
    @vprint :BruecknerSQ 2 "found quotient of order $(order(codomain(mp)))\n"
  end
end



#= issues/ TODO
 - does one need all SQs? are all maximal ones (in the sense of 
   no further extension possible) isomorphic?
 - part c: this will extend by the modules several times (a maximal
   number of times), useful if as above, any maximal chain will do
 - use (or not) compatible pairs to construct fewer possibilities
   (if legal...)
 - gmodules as gset, support the interface? in particular orbits?
   (for fewer extensions)
 - gmodule -> matrix group (in some cases)
 - filter for special targets (which???)
 - does this work with gmodules in char 0?

According to Max: 
 - if G is finite, then the maximal quotient is unique: the quotient
   modulo G'''' the infinite derived subgroup
 - if G is infinite and has a maximal quotient then it is also unique

So we can adjust the strategy accordingly.
In particular Satz 15 c should be implemented.
=#

#= EXAMPLE
Brueckner G1
F = @free_group(:a, :b)
G, hom = quo(F, [a^2*b*a^-1*b*a^-1*b^-1*a*b^-2,  a^2*b^-1*a*b^-1*a*b*a^-1*b^2])
f1 = Oscar.RepPc.solvable_quotient(G)
f2 = Oscar.RepPc.brueckner(f1; primes = [2])
Oscar.RepPc.brueckner(f2[1]; primes = [2])
f3 = ans;
f2 = Oscar.RepPc.brueckner(f1; primes = [2])
f3 = Oscar.RepPc.brueckner(f2[1]; primes = [2])
f4 = Oscar.RepPc.brueckner(f3[1]; primes = [3])
f5 = Oscar.RepPc.brueckner(f4[1]; primes = [2])
f6 = Oscar.RepPc.brueckner(f5[1]; primes = [2])
f7 = Oscar.RepPc.brueckner(f6[1], primes = [2])

(there is a [5] step missing)

H2, mH2, _ = Oscar.GrpCoh.H_two(C; redo = true, lazy = true);
T, mT = Oscar.GrpCoh.compatible_pairs(C)
G = gmodule(T, [Oscar.hom(H2, H2, [preimage(mH2, mT(g, mH2(a))) for a = gens(H2)]) for g = gens(T)])

then one wants the orbits of "elements" in G...

gset(matrix_group([matrix(x) for x= action(G)]))
@time orbits(ans)

this works, but constracts all elements...

G = gset(T, (a, g) -> preimage(mH2, mT(g, mH2(a))), collect(H2), closed = true)
orbits(G)

Brueckner G5
F = @free_group(:a, :b, :c)
G, _ = quo(F, [c*a*inv(c), inv(b)* a*b*a* b*a*c*b*inv(a)* inv(c)*b, c*inv(b)*a*c*b*c*inv(a)])

=#

end #module RepPc

using .RepPc
