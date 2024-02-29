module RepPc
using Oscar

export coimage

Base.pairs(M::MatElem) = Base.pairs(IndexCartesian(), M)
Base.pairs(::IndexCartesian, M::MatElem) = Base.Iterators.Pairs(M, CartesianIndices(axes(M)))

function Hecke.roots(a::FinFieldElem, i::Int)
  kx, x = polynomial_ring(parent(a), cached = false)
  return roots(x^i-a)
end

Oscar.matrix(phi::Generic.IdentityMap{<:AbstractAlgebra.FPModule}) = identity_matrix(base_ring(domain(phi)), dim(domain(phi)))


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

Implements: Brueckner, Chap 1.2.3
"""
function reps(K, G::Oscar.GAPGroup)
  @req is_finite(G) "the group is not finite"
  if order(G) == 1
    F = free_module(K, 1)
    h = hom(F, F, [F[1]])
    return [gmodule(F, G, typeof(h)[h for i = gens(G)])]
  end

  pcgs = GAP.Globals.Pcgs(G.X)
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
          # appliying det:
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
             needs to be lower left, and ubbber right the inverse.

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

Oscar.gen(M::AbstractAlgebra.FPModule, i::Int) = M[i]

Oscar.is_free(M::Generic.FreeModule) = true
Oscar.is_free(M::Generic.DirectSumModule) = all(is_free, M.m)

function Oscar.dual(h::Map{FinGenAbGroup, FinGenAbGroup})
  A = domain(h)
  B = codomain(h)
  @assert is_free(A) && is_free(B)
  return hom(B, A, transpose(h.map))
end

function Oscar.dual(h::Map{<:AbstractAlgebra.FPModule{ZZRingElem}, <:AbstractAlgebra.FPModule{ZZRingElem}})
  A = domain(h)
  B = codomain(h)
  @assert is_free(A) && is_free(B)
  return hom(B, A, transpose(matrix(h)))
end

function coimage(h::Map)
  return quo(domain(h), kernel(h)[1])
end

function Oscar.cokernel(h::Map)
  return quo(codomain(h), image(h)[1])
end

function Base.iterate(M::AbstractAlgebra.FPModule{T}) where T <: FinFieldElem
  k = base_ring(M)
  if dim(M) == 0
    return zero(M), iterate([1])
  end
  p = Base.Iterators.ProductIterator(Tuple([k for i=1:dim(M)]))
  f = iterate(p)
  return M(elem_type(k)[f[1][i] for i=1:dim(M)]), (f[2], p)
end

function Base.iterate(::AbstractAlgebra.FPModule{<:FinFieldElem}, ::Tuple{Int64, Int64})
  return nothing
end

function Base.iterate(M::AbstractAlgebra.FPModule{T}, st::Tuple{<:Tuple, <:Base.Iterators.ProductIterator}) where T <: FinFieldElem
  n = iterate(st[2], st[1])
  if n === nothing
    return n
  end
  return M(elem_type(base_ring(M))[n[1][i] for i=1:dim(M)]), (n[2], st[2])
end

function Base.length(M::AbstractAlgebra.FPModule{T}) where T <: FinFieldElem
  return Int(order(base_ring(M))^dim(M))
end

function Base.eltype(M::AbstractAlgebra.FPModule{T}) where T <: FinFieldElem
  return elem_type(M)
end

function Oscar.dim(M::AbstractAlgebra.Generic.DirectSumModule{<:FieldElem})
  return sum(dim(x) for x = M.m)
end

Oscar.is_finite(M::AbstractAlgebra.FPModule{<:FinFieldElem}) = true

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

  _ = Oscar.GrpCoh.H_two(C)
  ssc, mH2 = get_attribute(C, :H_two_symbolic_chain)
  sc = (x,y) -> ssc(x, y)[1]
  R = relators(G)
  M = C.M
  D, pro, inj = direct_product([M for i=1:ngens(G)]..., task = :both)
  a = sc(one(N), one(N))
  E = domain(a)
  DE, pDE, iDE = direct_product(D, E, task = :both)

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

   sc(a, b) yields sigma(a, b): E -> M, the value of the cohain as a map
   (ie. the once e in E is chosen, sc(a, b)(e) is THE cochain
  =#

  K, pK, iK = direct_product([M for i=1:length(R)]..., task = :both)
  s = hom(DE, K, [zero(K) for i=1:ngens(DE)])
  j = 1
  for r = R
    a = (one(N), hom(DE, M, [zero(M) for i=1:ngens(DE)]))
    for i = Oscar.GrpCoh.word(r)
      if i<0
        h = inv(mp(G[-i]))
        m = -pDE[1]*pro[-i]*action(C, h)  - pDE[2]*sc(inv(h), h)
      else
        h = mp(G[i])
        m = pDE[1]*pro[i]
      end
      # a *(h, m) = (x, y)(h, m) = (xh, m + y^h + si(x, h))
      a = (a[1]*h, m + a[2]*action(C, h) + pDE[2]*sc(a[1], h))
    end
    @assert isone(a[1])
    s += a[2]*iK[j]
    j += 1
  end
  #so kern(s) should be exactly all possible quotients that allow a
  #projection of G. They are not all surjective. However, lets try:
  k, mk = kernel(s)
  allG = []
  z = get_attribute(C, :H_two)[1]  #tail (H2) -> cochain

  seen = Set{Tuple{elem_type(D), elem_type(codomain(mH2))}}()
  #TODO: the projection maps seem to be rather slow - in particular
  #      as they SHOULD be trivial...
  for x = k
    epi = pDE[1](mk(x)) #the map
    chn = mH2(pDE[2](mk(x))) #the tail data
    if (epi,chn) in seen
      continue
    else
      push!(seen, (epi, chn))
    end
    #TODO: not all "chn" yield distinct groups - the factoring by the
    #      co-boundaries is missing
    #      not all "epi" are epi, ie. surjective. The part of the thm
    #      is missing...
    # (Thm 15, part b & c) (and the weird lemma)
#    @hassert :BruecknerSQ 2 all(x->all(y->sc(x, y)(chn) == last_c(x, y), gens(N)), gens(N))


    @hassert :BruecknerSQ 2 preimage(z, z(chn)) == chn
    GG, GGinj, GGpro, GMtoGG = Oscar.GrpCoh.extension(PcGroup, z(chn))
    @assert is_surjective(GGpro)
    if get_assertion_level(:BruecknerSQ) > 1
      _GG, _ = Oscar.GrpCoh.extension(z(chn))
      @assert is_isomorphic(GG, _GG)
    end

    function reduce(g) #in G
      h = mp(g)
      c = ssc(h, one(N))[2]
      if length(c) == 0
        return c
      end
      d = Int[abs(c[1]), sign(c[1])]
      for i=c[2:end]
        if abs(i) == d[end-1]
          d[end] += sign(i)
        else
          push!(d, abs(i), sign(i))
        end
      end
      return d
    end
    l= [GMtoGG(reduce(gen(G, i)), pro[i](epi)) for i=1:ngens(G)]
#    @show map(order, l), order(prod(l))
#    @show map(order, gens(G)), order(prod(gens(G)))

    h = try
          hom(G, GG, l)
        catch
          @show :crash
          continue
        end
    if !is_surjective(h)
#      @show :darn
      continue
    else
#      @show :bingo
    end
    push!(allG, h)
  end
  return allG
end

function solvable_quotient(G::Oscar.GAPGroup)
  q = cyclic_group(1)
  mp = hom(G, q, [one(q) for g in gens(G)])
end

end #module RepPc

using .RepPc

export coimage
