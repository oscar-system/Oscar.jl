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
function reps(K, G::Oscar.GAPGroup)
  @req is_finite(G) "the group is not finite"
  if order(G) == 1
    F = free_module(K, 1)
    h = hom(F, F, [F[1]])
    return [gmodule(F, G, typeof(h)[h for i = gens(G)])]
  end

  pcgs = GAP.Globals.Pcgs(GapObj(G))
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
        rh = gmodule(group(r), [action(r, preimage(ms, x^h)) for x = gens(s)])
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
          ii = findfirst(x->!iszero(x), Xp)
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
            for k in (pos+1):length(R)
              if length(Oscar.GModuleFromGap.hom_base(
                          gmodule(M, s, conjreps[j]), R[k])) > 0
                todo[k] = false
                continue
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
  return R
end


"""
Brueckner Chap 1.3.1

Given
  mp: G ->> Q

Find a set of primes suth that are any irreducible F_p module M
s.th. there is an epimorphism of G onto the extension of Q by M,
the p is in the set.
"""
function find_primes(mp::Map{<:Oscar.GAPGroup, PcGroup})
  G = domain(mp)
  Q = codomain(mp)
  if order(Q) == 1
    F = free_module(ZZ, 1)
    I = [gmodule(F, Q, [hom(F, F, [F[1]]) for x in gens(Q)])]
  else
    I = irreducible_modules(ZZ, Q)
  end
  lp = Set(collect(keys(factor(order(Q)).fac)))
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
      push!(lp, collect(keys(factor(order(t)).fac))...)
    end
  end
  return lp
end

"""
Given
    mQ: G ->> Q
Find all possible extensions of Q by an irreducible F_p module
that admit an epimorphism from G.
Implements the SQ-Algorithm by Brueckner, Chap 1.3

If necessary, the prime(s) p that can be used are computed as well.
"""
function brueckner(mQ::Map{<:Oscar.GAPGroup, PcGroup}; primes::Vector=[])
  Q = codomain(mQ)
  G = domain(mQ)
  @vprint :BruecknerSQ 1 "lifting $mQ using SQ\n"
  if length(primes) == 0
    @vprint :BruecknerSQ 1 "primes not provided, searching...\n"
    lp = find_primes(mQ)
  else
    lp = map(ZZRingElem, primes)
  end
  @vprint :BruecknerSQ 1 "using primes $lp\n"

  allR = []
  for p = lp
    _, j = ppio(order(Q), p)
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
      iii = Oscar.GModuleFromGap.gmodule(GF(Int(p)), ii)
      @vtime :BruecknerSQ 2 l = lift(iii, mQ)
      @vprint :BruecknerSQ 2 "found $(length(l)) many\n"
      append!(allR, [x for x in l])# if is_surjective(x)])
    end
  end
  return allR
end

"""
  mp: G ->> Q
  C a F_p[Q]-module
  Find all extensions of Q my C s.th. mp can be lifted to an epi.
"""
function lift(C::GModule, mp::Map)
  #m: G->group(C)
  #compute all(?) of H^2 that will describe groups s.th. m can be lifted to

  G = domain(mp)
  N = group(C)
  @assert isa(N, PcGroup)
  @assert codomain(mp) == N

  H2, z, _ = Oscar.GrpCoh.H_two(C; lazy = true)
  if order(H2) > 1
    global last_in = (C, mp)
  end
  R = relators(G)
  M = C.M
  D, pro, inj = direct_product([M for i=1:ngens(G)]..., task = :both)

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

  K, pK, iK = direct_product([M for i=1:length(R)]..., task = :both)
  S = relators(N)
  if length(S) != 0
    X, pX, iX = direct_product([M for i=1:length(S)]..., task = :both)
  end

  allG = []

  for h = H2
    GG, GGinj, GGpro, GMtoGG = Oscar.GrpCoh.extension(PcGroup, z(h))

    s = hom(D, K, [zero(K) for i=1:ngens(D)])
    gns = [GMtoGG([x for x = GAP.Globals.ExtRepOfObj(GapObj(h))], zero(M)) for h = gens(N)]
    gns = [map_word(mp(g), gns, init = one(GG)) for g = gens(G)]

    rel = [map_word(r, gns, init = one(GG)) for r = relators(G)]
    @assert all(x->isone(GGpro(x)), rel)
    rhs = [preimage(GGinj, x) for x = rel]
    s = hom(D, K, [K([preimage(GGinj, map_word(r, [gns[i] * GGinj(pro[i](h)) for i=1:ngens(G)])) for r = relators(G)] .- rhs) for h = gens(D)])

    k, mk = kernel(s)

    if is_zero(h) && length(S) > 0 #Satz 15 in Brueckner, parts a and b
      u = [preimage(mp, x) for x = gens(N)]
      ss = elem_type(X)[]
      sss = elem_type(D)[]
      for h = gens(D)
        uu = [map_word(x, [gns[i] * GGinj(pro[i](h)) for i=1:ngens(G)]) for x = u]
        push!(ss, X([preimage(GGinj, map_word(r, uu, init = one(GG))) for r = relators(N)]))
        push!(sss, D([preimage(GGinj, map_word(mp(G[i]), uu, init = one(GG))*inv(gns[i] * GGinj(pro[i](h)))) for i = 1:ngens(G)]))
      end
      ss = hom(D, X, ss)
      sss = hom(D, D, sss)
      kk, mkk = intersect(kernel(ss)[1], kernel(sss)[1])
      q, mq = quo(k, kk)
      k = q
      mk = pseudo_inv(mq)*mk
    end

    fl, pe = try
      true, preimage(s, K(rhs))
    catch
      false, zero(D)
    end
    if !fl
      @show :no_sol
      continue
    end

    for x = k
      if is_zero(h) && is_zero(x)
        @show :skip, dim(k)
        continue
      end
      hm = hom(G, GG, [gns[i] * GGinj(pro[i](-pe +  mk(x))) for i=1:ngens(G)])
      if is_surjective(hm)
        push!(allG, hm)
      else #should not be possible any more
        @show :not_sur
      end
    end
  end

  return allG
end

function solvable_quotient(G::Oscar.GAPGroup)
  q = cyclic_group(1)
  mp = hom(G, q, [one(q) for g in gens(G)])
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
=#

#= EXAMPLE
F = free_group("a", "b"); a,b = gens(F);
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

then ones wants the orbits of "elements" in G...

gset(matrix_group([matrix(x) for x= action(G)]))
@time orbits(ans)

this works, but constracts all elements...

G = gset(T, (a, g) -> preimage(mH2, mT(g, mH2(a))), collect(H2), closed = true)
orbits(G)

=#

end #module RepPc

using .RepPc
