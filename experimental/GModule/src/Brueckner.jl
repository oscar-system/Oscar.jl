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

`dim_bound` restricts the search to representations of at most that
dimension. Dimensions only grow along the pc chain, so a bound prunes whole
branches rather than filtering at the end - which is the difference between
answering and not answering for a large `G`.

Note: the reps are NOT necessarily over the smallest field.

Note: the field is NOT extended - but it throws an error if it was too small.

Note: `group(M)` for the returned gmodules `M` will have a pcgs of `G` as
      its `gens` value, thus these generators will in general differ from
      the generators of `G`.

Implements: Brueckner, Chap 1.2.3
"""
function reps(K, G::Oscar.PcGroup; dim_bound::Int = typemax(Int))
  @req is_finite(G) "the group is not finite"
  @req dim_bound >= 1 "the dimension bound has to be positive"
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
          # inducing multiplies the dimension by `p` and nothing further down
          # the chain ever shrinks it, so this branch is out of range for good.
          # The `h`-conjugates of `r` have the same dimension, so they drop out
          # here too and need not be marked as done.
          n*p > dim_bound && continue

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
              if is_isomorphic(Gj, R[k])
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
    _admissible_primes(mQ::Map)

Given `mQ: G ->> Q`, the primes `p` for which `Q` has an extension by an
irreducible `F_p[Q]`-module admitting an epimorphism from `G`.

Those are exactly the primes dividing the torsion of `N^ab`, `N = ker(mQ)`: any
such module is a quotient of `N`, and conversely a non-zero `N/[N,N]N^p` has an
irreducible quotient.
"""
function _admissible_primes(mQ::Map{<:Oscar.GAPGroup, PcGroup})
  inv = abelian_invariants(kernel(mQ)[1])
  @req !(0 in inv) "the kernel has infinite abelianization, so every prime is admissible; pass `primes`"
  return sort!(unique!(reduce(vcat, (prime_divisors(x) for x in inv), init = ZZRingElem[])))
end

"""
Brueckner Chap 1.3.1

Given
  mp: G ->> Q

Find a set of primes such that if there is any irreducible F_p module M
s.th. there is an epimorphism of G onto the extension of Q by M,
the p is in the set.
"""
function find_primes(mp::Map{<:Oscar.GAPGroup, PcGroup})
  G = domain(mp)
  Q = codomain(mp)
  if order(Q) == 1
    F = free_module(ZZ, 1)
    I = [gmodule(F, Q, [hom(F, F, [F[1]]) for x in gens(Q)])]
  else #TODO: repsn, reps offer choice
#    I = irreducible_modules(ZZ, Q)
    I = reps(abelian_closure(QQ)[1], Q)
    #Brueckner, p35: irreducible here is not necessary, so the
    #  expensive find minimal field step can be omitted.
    I = [gmodule(ZZ, gmodule(QQ, gmodule(CyclotomicField, x))) for x = I]
  end
  lp = Set(prime_divisors(order(Q)))
  for i = I
    ib = gmodule(i.M, G, [action(i, mp(g)) for g = gens(G)])
    ia = gmodule(FinGenAbGroup, ib)
    a, b = Oscar.GrpCoh.H_one_maps(ia)
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
    q = cokernel(b)[1]
#    q = quo(kernel(da)[1], image(db)[1])[1]
    t = torsion_subgroup(q)[1]
    if order(t) > 1
      push!(lp, prime_divisors(order(t))...)
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
function brueckner(mQ::Map{<:Oscar.GAPGroup, PcGroup}; primes::Vector=[], limit::Int = typemax(Int))
  Q = codomain(mQ)
  G = domain(mQ)
  @vprint :BruecknerSQ 1 "lifting $mQ using SQ\n"

  allR = []

  # collect lifts for the given primes, from modules of dimension in
  # `lo+1:hi`; `true` once `limit` of them are known
  function _extend_by(lp::Vector{ZZRingElem}, lo::Int, hi::Int)
    @vprint :BruecknerSQ 1 "using primes $lp, dimensions $(lo+1) to $hi\n"
    for p in lp
      _, j = ppio(exponent(Q), p)
      f = j == 1 ? 1 : modord(p, j)
      @assert (p^f-1) % j == 0
      @vprint :BruecknerSQ 2 "computing reps over GF($p, $f)\n"
      if f == 1
        @vtime :BruecknerSQ 2 I = reps(GF(Int(p)), Q; dim_bound = hi)
      else
        @vtime :BruecknerSQ 2 I = reps(GF(Int(p), f), Q; dim_bound = hi)
      end
      @vprint :BruecknerSQ 1 "have $(length(I)) representations\n"

      for i in I
        dim(i) > lo || continue     # already tried in an earlier round
        @vprint :BruecknerSQ 1 "starting to process module\n"
        @vprint :BruecknerSQ 2 "... transfer over min. field\n"
        @vtime :BruecknerSQ 2 ii = Oscar.GModuleFromGap.gmodule_minimal_field(i)
        @vprint :BruecknerSQ 2 "... lift...\n"
        #TODO: why do we need the module over GF(p)???
        iii = Oscar.GModuleFromGap.gmodule(GF(Int(p)), ii)
        @vtime :BruecknerSQ 2 l = lift(iii, mQ; limit = limit - length(allR))
        @vprint :BruecknerSQ 2 "found $(length(l)) many\n"
        #TODO: in Plesken p119 has more comments what not to do
        append!(allR, [x for x in l])# if is_surjective(x)])
        length(allR) >= limit && return true
      end
    end
    return false
  end

  #= The dimension bound prunes the pc chain rather than the answer, so asking
     for small modules first is far cheaper than asking for all of them and
     is usually where the next layer is anyway. Walk it upward, and pay for
     the complete set only if nothing turned up. A caller that wants every
     lift has to pay for it either way and goes straight there.
  =#
  function _search(lp::Vector{ZZRingElem})
    limit == typemax(Int) && return _extend_by(lp, 0, typemax(Int))

    lo = 0
    for hi in (1, 2, 4)
      _extend_by(lp, lo, hi) && return true
      lo = hi
    end
    return _extend_by(lp, lo, typemax(Int))
  end

  if length(primes) > 0
    _search(map(ZZRingElem, primes))
    return allR
  end

  #= `_admissible_primes` needs the kernel, hence an enumeration of the |Q|
     cosets. A caller that stops early can often avoid that: the primes
     dividing |Q| are free to name, and one of them usually does lift. Only
     the exact set is guaranteed complete, so an exhaustive caller goes
     straight there.
  =#
  cheap = limit == typemax(Int) ? ZZRingElem[] : sort(prime_divisors(order(Q)))
  _search(cheap) && return allR

  @vprint :BruecknerSQ 1 "primes not provided, searching...\n"
  _search(setdiff(_admissible_primes(mQ), cheap))
  return allR
end

"""
  mp: G ->> Q
  C a F_p[Q]-module
  Find all extensions of Q my C s.th. mp can be lifted to an epi.
"""
function lift(C::GModule, mp::Map; limit::Int = typemax(Int))
  #m: G->group(C)
  #compute all(?) of H^2 that will describe groups s.th. m can be lifted to

  G = domain(mp)
  N = group(C)
  @assert isa(N, PcGroup)
  @assert codomain(mp) == N
  # the surjectivity argument in `_process` needs `M` to have no proper
  # non-zero submodule
  @req dim(C) == 1 || is_irreducible(C) "the module has to be irreducible"

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

  D, pro, inj = direct_product([M for i in 1:ngens(G)]..., task = :both)
  # `direct_product` needs at least one factor; a presentation without relators
  # imposes no conditions on the derivations
  K = is_empty(R) ? free_module(base_ring(M), 0) :
                    direct_product([M for i in 1:length(R)]..., task = :none)

  # |Z^1(N, M)|, the number of lifts that miss `M`; independent of the cocycle
  ordZN = ngens(N) == 0 ? ZZ(1) :
                          order(kernel(Oscar.GrpCoh.H_one_maps(C)[2])[1])

  # the canonical lifts of the generators of `G` into an extension
  function _lifted_gens(ext)
    GG, _, _, GMtoGG = ext
    gns = [GMtoGG([x for x in Oscar.GAPWrap.ExtRepOfObj(GapObj(h))], zero(M)) for h in gens(N)]
    return [map_word(mp(g), gns, init = one(GG)) for g in gens(G)]
  end

  # by how much the relators of `G` miss being satisfied there
  function _defect(ext, gns)
    GG, GGinj, GGpro, _ = ext
    rel = [map_word(r, gns, init = one(GG)) for r in R]
    @assert all(x->isone(GGpro(x)), rel)
    return K([preimage(GGinj, x) for x in rel])
  end

  #= Replacing `gns[i]` by `gns[i]*m` moves the relator defects by a map that
     is linear in the `m` and built only from the action of `N` on `M` and the
     relators: the cocycle enters the defects as a constant. So `s` below, and
     with it Z^1(G, M), is the same for every class, and is worth computing
     once - it costs ngens(D) evaluations of every relator in `GG`, which for
     a large quotient dwarfs everything else here.
  =#
  ext0 = Oscar.GrpCoh.split_extension(PcGroup, C)
  gns0 = _lifted_gens(ext0)
  @hassert :BruecknerSQ 1 is_zero(_defect(ext0, gns0))
  s = hom(D, K, [K([preimage(ext0[2], map_word(r, [gns0[i] * ext0[2](pro[i](h)) for i in 1:ngens(G)])) for r in R]) for h in gens(D)])
  k, mk = kernel(s)

  # `pe` solves s(pe) = defect, so the twist -pe kills the relators
  function _process(ext, pe; is_trivial::Bool = false, limit::Int)
    GG, GGinj, GGpro, _ = ext
    res = typeof(mp)[]
    @assert isa(GG, PcGroup)

    gns = _lifted_gens(ext)
    @hassert :BruecknerSQ 1 s(pe) == _defect(ext, gns)

    #= The lifts form a torsor under Z^1(G, M) = `k`. Such a lift misses `M`
       iff its image is a complement to `M` in `GG`, since the image meets `M`
       in a submodule of the irreducible `M`. Only a split extension has
       complements, and there the offending lifts are exactly the subspace
       inf(Z^1(N, M)) of `k`. Hence:
         - for a non-trivial class no surjectivity test is needed at all;
         - for the trivial one, a surjective lift exists iff `k` is bigger than
           that subspace, and then some generator of `k` lies outside it.
       Trying the generators first matters: `k` is enumerated along an rref
       basis, so the subspace can occupy a long prefix.
    =#
    is_trivial && order(k) == ordZN && return res

    function _try(x)
      hm = hom(G, GG, [gns[i] * GGinj(pro[i](-pe + mk(x))) for i in 1:ngens(G)])
      if is_trivial
        is_surjective(hm) || return false
      else
        @hassert :BruecknerSQ 1 is_surjective(hm)
      end

      push!(res, hm)
      return length(res) >= limit
    end

    tried = elem_type(k)[]
    if is_trivial
      for x in gens(k)
        push!(tried, x)
        _try(x) && return res
      end
    end

    for x in k
      x in tried && continue
      _try(x) && return res
    end
    @hassert :BruecknerSQ 1 length(res) == order(k) - (is_trivial ? ordZN : 0)
    return res
  end

  allG = _process(ext0, zero(D); is_trivial = true, limit)
  if length(allG) >= limit || gcd(order(C.G), order(C.M)) == 1 #trivial H^2
    return allG
  end

  H2, z, _ = Oscar.GrpCoh.H_two(C; lazy = true)

  #= By Schur, E = End_{F_p[N]}(M) is a field, and each of its units is an
     automorphism of `M` commuting with the action of `N`, so it pairs with
     the identity on `N`. The isomorphism such a pair induces between the
     extensions for `h` and `u*h` is then the identity on `N`, so it carries
     lifts of `mp` to lifts of `mp`, bijectively and preserving both
     surjectivity and the kernel: one class per E-line through 0 describes
     every quotient that the whole line does. H^2 is an E-vector space, so
     that line is the F_p-span of the images of `h` under an F_p-basis of E.

     E-lines are as far as this goes. `Oscar.GrpCoh.compatible_pairs` gives
     coarser orbits on H^2, but they buy nothing: (a, b) fixes the surjection
     `mp` only for b = id, and those pairs are exactly the units of E. For
     b != id the isomorphism E(h) -> E((a,b)*h) induces `b` on `N`, so lifts
     of `mp` there are lifts of b^-1*mp here and have to be recovered by
     transporting `mp`. The units of E act freely on H^2 minus 0, hence meet
     no stabiliser, so an orbit of size s needs s/(|E|-1) transported maps -
     exactly the number of E-lines it contains. The number of lifts to
     compute is therefore the same, and all that orbits would save is
     building each extension once per orbit instead of once per line, which
     does not pay for the `automorphism_group(M)` inside `compatible_pairs`.
     (Thm 15, part b & c) (and the weird lemma)
  =#
  p = Int(characteristic(base_ring(M)))
  S = elem_type(N)
  T = elem_type(M)
  endo = [hom(M, M, b) for b in Oscar.GModuleFromGap.hom_base(C, C)]

  #= The action of E on H^2 is F_p-linear, so get it once as maps rather than
     transporting a cochain through `z` for every line: that is one pass over
     the generators of H^2 instead of one per line. Pushing the cochain of a
     generator forward along every basis element of E before moving on keeps
     the values it memoised while being evaluated.

     Over the prime field the line through `h` is spanned by `h`, and no
     cochain has to be transported at all.
  =#
  endo_H2 = if length(endo) == 1
    [id_hom(H2)]
  else
    imgs = map(gens(H2)) do g
      c = z(g)
      [preimage(z, Oscar.GrpCoh.CoChain{2, S, T}(C, Dict{NTuple{2, S}, T}(),
                                                 x -> b(c(x[1], x[2]))))
       for b in endo]
    end
    [hom(H2, H2, [imgs[i][j] for i in 1:ngens(H2)]) for j in 1:length(endo)]
  end

  function _line(h)
    line = [zero(H2)]
    for f in endo_H2
      g = f(h)
      line = [x + l*g for x in line for l in 0:p-1]
    end
    return line
  end

  #= The defects are linear in the cocycle as well - a relator evaluates to a
     sum of cocycle values moved around by the action - so read them off a map
     built from the generators of H^2. Deciding whether a class lifts is then
     linear algebra, and only the classes that do lift need their extension
     built.
  =#
  rhs_H2 = hom(H2, K, elem_type(K)[_defect(e, _lifted_gens(e)) for e in
                       (Oscar.GrpCoh.extension(PcGroup, z(g)) for g in gens(H2))])

  seen = Set{elem_type(H2)}()

  for h in H2
    is_zero(h) && continue
    h in seen && continue
    union!(seen, _line(h))

    fl, pe = has_preimage_with_preimage(s, rhs_H2(h))
    fl || continue

    append!(allG, _process(Oscar.GrpCoh.extension(PcGroup, z(h)), pe; is_trivial = false, limit = limit - length(allG)))
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
  end

  # no maximal finite abelian quotient to start from
  q = cyclic_group(1)
  return hom(G, q, [one(q) for g in gens(G)])
end

function sq(mp::Map, primes::Vector=[]; index::Union{Integer, ZZRingElem, Nothing} = nothing)
  if index === nothing
    @req is_finite(maximal_abelian_quotient(domain(mp))[1]) "infinite abelianization: there is no maximal finite solvable quotient; pass `index`"
  end

  if index !== nothing
    lf = factor(ZZRingElem(index))
    primes = prime_divisors(ZZ(index))
    while length(primes) > 0
      @vtime :BruecknerSQ 1 nw = brueckner(mp; limit = 1, primes)
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
    nw = brueckner(mp; limit = 1, primes)
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
F = @free_group(:a, :b)
G, hom = quo(F, [a^2*b*a^-1*b*a^-1*b^-1*a*b^-2,  a^2*b^-1*a*b^-1*a*b*a^
-1*b^2])
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

=#

end #module RepPc

using .RepPc
