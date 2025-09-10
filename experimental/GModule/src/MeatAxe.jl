module MeatAxe
using Oscar
using Oscar.GModuleFromGap
import Oscar: GAPGroup
import Oscar.GModuleFromGap: is_regular_gmodule, restrict_endo


#very basic...
function gmodule_new(chi::Oscar.GAPGroupClassFunction)
  G = group(chi)
  s = [representative(x) for x = subgroup_classes(group(chi))]
  sort!(s, lt = (a,b) -> isless(order(b), order(a)))
  for i = s
    if !is_zero(scalar_product(permutation_character(G, i), chi))
      M = permutation_gmodule(G, i, QQ)
      t = split_into_homogenous(M)
      tt = [character(domain(x)) for x = t]
      return [domain(t[i]) for i=1:length(t) if !is_zero(scalar_product(tt[i], chi))]
    end
  end
  error("no plan yet")
end

#tests if any component of M gives a new representation
# res is the list of known reps
# t[i] is true if the rep is nor known. Index as in the character table
function _do_one!(res, t, M; is_irr::Bool = false)
  if !is_irr
    _t = split_into_homogenous(M)
  else
    _t = [hom(M, M, hom(M.M, M.M, gens(M.M)))]
  end
  for _x = _t
    x = domain(_x)
    _c = map(iszero, coordinates(character(x)))
    if t .& _c == t
#      @show :bad, "module known", _c
      continue
    end
#    @show "found", _c
    t .&= map(iszero, coordinates(character(x)))
    push!(res, x)
    if !any(t)
      return true
    end
  end
  return false
end

#for tensor products where the combined degree is bounded my limit
#try to decompose...
#as usual: do it with characters first
#TODO: if they are not G-modules...
function _try_tensor_products!(res::Vector{GModule}, t::Vector{Bool}; limit::Int = 100)
  for i = 1:length(res)
    for j = 2:length(res)
      if dim(res[i])*dim(res[j]) <= limit
        c = coordinates(tensor_product(character(res[i]), character(res[j])))
        cc = map(iszero, c) .& t
        if cc != t
#          @show :bingo
          M = tensor_product(res[i], res[j]; task = :none)
          _do_one!(res, t, M) && return true
        end
      end
    end
  end
  return false
end

function _try_squares!(res::Vector{GModule}, t::Vector{Bool}; limit::Int = 100)
  for i = 1:length(res)
    if 0 < dim(res[i])*(dim(res[i])-1)/2 <= limit
      c = coordinates(alternating_square(character(res[i])))
      cc = map(iszero, c) .& t
      M = alternating_square(res[i])
      _do_one!(res, t, M) && return true
    end
    if dim(res[i])*(dim(res[i])+1)/2 <= limit
      c = coordinates(symmetric_square(character(res[i])))
      cc = map(iszero, c) .& t
      M = symmetric_square(res[i])
      _do_one!(res, t, M) && return true
    end
  end
  return false
end

function _try_perm_character!(res::Vector{GModule}, t::Vector{Bool}, G::GAPGroup, U::GAPGroup; limit::Int = 100)
  c = coordinates(permutation_character(G, U))
  cc = map(iszero, c) .& t
  if cc != t #will give new rep.
#    @show :bingo
    M = permutation_gmodule(G, U, QQ)
    _do_one!(res, t, M) && return true
  end
  return false
end


#"main" function: try to find all irreducible rationals representations of G
# limit "limits" the size of the modules constructed on the way
#   (endomorphism rings are too expensive)
function gmodule_irred_rational(G::Group; limit::Int = 50, t::Union{Nothing, Vector{Bool}} = nothing)
  if is_abelian(G)
    return Oscar.GModuleFromGap._irred_abelian(G)
  end
  #linear characters are handled elsewhere. Here the expectation is
  #that in "t" the linear characters are marked as known
  #the code works for linear characters, it's just unneccessarily slow)
  X = character_table(G)
  s = [representative(x) for x = low_index_subgroup_classes(G, limit)]
  s = [x for x = s if order(x) < order(G)]
  sort!(s, lt = (a,b) -> isless(order(b), order(a)))
  if t !== nothing
    @assert length(t) == length(X)
  else
    t = [degree(x) <= limit for x = X]
  end
  res = GModule[]
  #permutation gmodules are good as the endomorphism ring
  #is "for free", the endo is used for the splitting
  #poss. the limits should be adjusted as this is easier...
  for i = s
    _try_perm_character!(res, t, G, i) && return res
  end
#  @show "after permutation modules", t

  _try_tensor_products!(res, t; limit) && return res
  _try_squares!(res, t; limit) && return res

  all_K = [x for x = subgroup_reps(G) if order(x) <= 20] #how to do this?

  for i=s
#    @show "order $(order(i))"
    XX = character_table(i)
    for chi = XX
      #careful with multiplicity
      deg = degree_of_character_field(chi)*degree(chi)*index(G, i) 
      if deg > limit^2
#        @show :too_large, degree_of_character_field(chi), degree(chi), index(G, i)
        continue
      end
      chi_G = induce(chi, G) #could be irreducible...
      cc = map(iszero, coordinates(chi_G)) .& t
      if cc != t
#        @show "could help"
        v = gmodule_new(chi)
        @assert length(v) == 1
        if is_irreducible(chi_G)
          V = induce(v[1], embedding(i, G))[1]
          _do_one!(res, t, V; is_irr = true) && return res
          continue
        end
        v = v[1]
        if deg <= limit
#          @show "direct induce"
          for j=1:length(t)
            if cc[j] != t[j] 
#              @show :induce
              V = induce(v, embedding(i, G))[1]
              _do_one!(res, t, V) && return res
            else
#              @show "bad"
            end
          end
        else
          local V
          have_V = false
          good_K = []
          for K = all_K
            dim_C = scalar_product(restrict(chi_G, K), trivial_character(K))
            if dim_C <= limit
              #test if K makes sense
              tt = [!t[i] || is_zero(scalar_product(restrict(X[i], K), trivial_character(K))) for i= 1:length(t)]
              if all(tt)
#                @show "condense would give no info"
                continue
              end
              push!(good_K, (K, dim_C, tt))
            end
          end
          if length(good_K) == 0
            continue
          end
          for ii = findall(!, reduce(.&, [x[3] for x = good_K]))
            t[ii] || continue #due to orbits
            char = galois_orbit_sum(X[ii])
            if degree(char) > limit
#              @show "too large, abandoning"
              continue
            end
#            @show "aiming for $i"
            for_i = [x for x = good_K if !x[3][ii]]
            sort!(for_i, lt = (a,b) -> a[2] < b[2])
            rand_G = [rand(G) for i=1:20]
#            @show length(for_i)
            K = for_i[1][1]
            if for_i[1][2] == 0
#              @show "oh shit"
              continue
            end
            if !have_V
              V = induce(v, embedding(i, G))[1]
              have_V = true
            end
#            @show "condense", V, K
            c, mc, phi = condense(V, K)
            @assert dim(c) <= limit
            sp = split_into_homogenous(c)
#            @show "split done"
            #TODO: use the traces in some form to decide which element in
            #      sp needs spinning - generically, there should be only one
            #      this is via Allan's trace vectors
            #      not sure what happens if split_into_homo has multiplicities
            #TODO: trace vectors: Allan adds more "generators" to the condensed
            #      module to have traces and to avoid problems from
            #      condensing Q[G] wrongly. Not sure if this is the correct
            #      way...
            tr = [trace_condensed(char, K, g) for g = rand_G]
            for _t = sp
              ttr = Any[]
              for g = rand_G
                _c = trace_condensed(V, _t.module_map*mc, phi, g)
                if isnothing(_c)
                  @show "rand_G too small..."
                  break
                end
                push!(ttr, _c)
                if ttr[end] != tr[length(ttr)]
#                  @show "stop as trace is wrong"
                  break
                end
              end
              @show ttr
              if tr == ttr
#                @show "bingo!!", degree(char)
                d = spin2(V, [mc(_t(x).data) for x = gens(domain(_t))]; limit = Int(degree(char)))
                if dim(d) > degree(char) || d == V
#                  @show coordinates(character(d))
#                  @show coordinates(char)
#                  @show coordinates(character(V))
#                  @show "spin too large", d, dim(d), limit, dim(V)
                  #OK: rand_G is too small
                  break
                  error("should not happen")
                end
                @assert dim(d) < dim(V)
#                @show dim(d), dim(V), coordinates(character(d))
                if dim(d) <= limit
                  _do_one!(res, t, d; is_irr = true) && return res
                end
              else
#                @show "wrong module lifted"
              end
              if !any(cc .& t)
#                @show "OK, leaving this level"
                #can't get any more out of V
                break
              end
            end
            if !any(cc .& t)
              #can't get any more out of V
              break
            end
          end
        end
      end
    end
  end
  return res
  error("no plan yet")
end

#a GModule over a (normal) number field, find the inequivalent 
#conjugate modules (should give a Galois-orbit if the characters)
function galois_orbit(C::GModule)
  k = base_ring(C)
  if k == QQ
    return [C]
  end
  a = automorphism_list(k) #can go wrong for Schur index
                           #if the field is not normal over Q
  z = Any[]
  cz = []
  for x = a
    t = C^x
    @assert Oscar.GrpCoh.is_consistent(t)
    ct = character(t)
    if ct in cz
      continue
    else
      push!(z, t)
      push!(cz, ct)
    end
  end
  return z
end

#C should be irreducible, rational. Find a unique
#abs. irred. in there - the others should be Galois conjugate
function abs_irred(C::GModule)
  if is_irreducible(character(C))
    return C
  end
  s = split_homogeneous(C)
  if length(s) == 0
    s = [C]
  else
    s = [domain(x) for x = s]
  end
  for y = s #add the Galois orbits...
    if is_irreducible(character(y))
      return y
    else
      t = split_homogeneous2(y)
      return t
    end
  end
end

function gmodule_new(G::Group; limit::Int = 50, t::Union{Nothing, Vector{Bool}} = nothing)
  X = character_table(G)
  if is_abelian(G)
    return Oscar.GModuleFromGap._abs_irred_abelian(G; cyclo = true)
  end
  #linear characters come from the max. ab. quot.
  Gab, mGab = maximal_abelian_quotient(G)
  zz = []
  z = Oscar.GModuleFromGap._abs_irred_abelian(Gab; cyclo = true)
  for x = z
    push!(zz, inflate(x, mGab))
  end

  if isnothing(t)
    t = [degree(x) > 1 for x = X]
  else
    t .&= [degree(x) > 1 for x = X]
  end

  z = gmodule_irred_rational(G; limit, t)
  for x = z
    append!(zz, galois_orbit(abs_irred(x)))
  end
  return zz
end

#if the Q[G]-module C is fixpoint condensed with K, the dim will be:
function dim_condensed(C::GModule, K::Oscar.GAPGroup)
  #Allan, Cor. 3.5.5.
  return Int(scalar_product(restrict(character(C), K), trivial_character(K)))
end

function trace_condensed(chi::Oscar.GAPGroupClassFunction, K::Oscar.GAPGroup, g::GAPGroupElem)
  return QQ(ZZ(1), order(K))*sum(chi(g*k) for k = K)
end

function trace_condensed(C::GModule, ms::Map, phi, g::GAPGroupElem)
  s = domain(ms)
  has, pre = has_preimage_with_preimage(ms, map(phi, action(C, g, map(ms, gens(s)))))
  has || return nothing
  return trace(hom(s, s, pre))
end

function idempotent(C::GModule, iKG::Map, v::AbstractAlgebra.FPModuleElem)
  #use a pc-rep tp reduce the calls to action and the number of terms
  #if C is huge and K small, using action on elements will be better
  #once action without matrix is implemented...
  return QQ(1, order(domain(iKG)))*sum(action(C, iKG(k), v) for k = domain(iKG))
end

function idempotent(C::GModule, iKG::Map, v::Vector{<:AbstractAlgebra.FPModuleElem})
  #use a pc-rep tp reduce the calls to action and the number of terms
  #if C is huge and K small, using action on elements will be better
  #once action without matrix is implemented...
  s = copy(v)
  for k = domain(iKG)
    isone(k) && continue
    s .+= action(C, iKG(k), v)
  end
  return QQ(1, order(domain(iKG))).*s
end

function idempotent(C::GModule, iKG::Map)
  #use a pc-rep tp reduce the calls to action and the number of terms
  K = domain(iKG)
  if is_solvable(K)
    hom = isomorphism(PcGroup, K) #K -> pc
    g = gens(codomain(hom))
    r = [relative_order(x) for x = g]
    a = [action(C, iKG(preimage(hom, x))) for x = g]
    I = id_hom(C.M)
    for i=length(r):-1:1
      a = action(C, iKG(preimage(hom, g[i])))
      J = id_hom(C.M) + a
      b = a
      for j = 2:relative_order(g[i])-1
        b *= a
#        @assert b == action(C, iKG(preimage(hom, g[i]^j)))
        J += b
      end
#      @assert J == sum([action(C, iKG(preimage(hom, g[i]^l))) for l=0:relative_order(g[i])-1])
      I = J*I
    end
    return QQ(1, order(domain(iKG)))*I
  end
  return QQ(1, order(domain(iKG)))*sum([action(C, iKG(k)) for k = domain(iKG)])
end

#ms: S -> M Q-modules
#try to find a nicer version by
#  - making it integral (clear denominators)
#  - saturating
#  - LLL-reducing.
#the hope is that the action matrices will be nice
function _simplify(ms::Map) #embedding map s -> M
  s = domain(ms)
  M = codomain(ms)
  x = matrix(ms)
  n = numerator(x)
  @show maximum(nbits, n)
  n = saturate(n)
  n = lll(n)
  @show maximum(nbits, n)
  return hom(s, M, matrix(QQ, n))
end

#condese the Q[G]-module C for the subgroup K
#returns 
# the A-module (without a group!!!, A = condensed Q[G])
# the vector-space embedding
# the idempotent
function condense(C::GModule, K::Oscar.GAPGroup; extra::Int = 5)
  #K should be a subgroup of G for the G-Module C
  G = group(C)
  iKG = embedding(K, G)
  #phi = 1/|K| sum_K k in Z[G] should be an idempotent
  #we want to get C phi as a module - without a group, this is algebra
  #only. The operation should be 
  #  phi g phi
  d = dim_condensed(C, K)
  M = C.M
  i = 1
  #TODO: spinning properly - or different strategy
  #      here I "condense" along a basis, alternatively, or better
  #      spin the elements to get operation and larger images
  #      there is/ might be the problem that gens(G) is too small,
  #      ie. the condensation is not the condenstation of Q[G]
  #      we want the image and the operation...
  phi = idempotent(C, iKG)
  ge = elem_type(C.M)[]
  s, ms = sub(C.M, ge)
  i = 0
  while length(ge) < d
    i += 1
    x = phi(C.M[i])
    if has_preimage_with_preimage(ms, x)[1]
      continue
    end
    push!(ge, x)
    s, ms = sub(C.M, ge)
  end
  #try to reduce/ improve
  ms = _simplify(ms)
  
  #action and preimage allows vectors
  #the "condense" problem for Q[G]: we don't know how many 
  # (and which) elements of G generate phi(QQ[G]).
  # Not having enough will meant preimage will fail to work
  # in the trace_condense...(where we try to compute more action)
  #TODO: sanity: for everything in rG, the traces will be "for free"
  #      as the action is already computed. so use it?
  rG = vcat(gens(G), [rand(G) for i=1:extra])
  return gmodule(nothing, [hom(s, s, preimage(ms, map(phi, action(C, g, map(ms, gens(s)))))) for g = rG]), ms, phi
end

#plain, vanilla spin, slow
function spin(C::GModule, v::Vector)
  s, ms = sub(C.M, v)
  v = map(ms, gens(s))
  done = false
  while !done
    done = true
    for m = v, g = gens(group(C))
      x = action(C, g, m)
      if has_preimage_with_preimage(ms, x)[1]
        continue
      end
      done = false
      push!(v, x)
      s, ms = sub(C.M, v)
      v = map(ms, gens(s))
      break
    end
  end
  return gmodule(group(C), [hom(s, s, preimage(ms, action(C, g, v))) for g = gens(group(C))])
end

#same as above, but bypass modules and work with matrices directly
function spin2(C::GModule, v::Vector; limit::Int = 50)
  #TODO: not use action(C, g, ...) but C.ac[1] in some form
  #      this might also work for non groups!
  #TODO: use similar for condensation?
  G = group(C)
  s = vcat([x.v for x = v]...)
  r = rref!(s)  
  piv = Int[]
  i = 1
  for j=1:r
    while is_zero_entry(s, j, i)
      i += 1
    end
    push!(piv, i)
  end
  new = []
  if base_ring(C) == QQ
    for i = 1:r
      t = s[i:i, :]
      mul!(t, inv(content(t)))
      st = maximum(nbits, t)
      for g = gens(G)
        push!(new, (t, g, st))
      end
    end
  end

  x = zero_matrix(QQ, 1, ncols(s))
  while length(new) > 0
    sort!(new, lt = (a,b) -> a[3] > b[3])
    m, g, sm = pop!(new) #order matters! try to use small elements...
    mul!(x, m, matrix(action(C, g)))  #if action is matrix free...
    reduce_mod!(x, view(s, 1:r, :))
    if is_zero(x)
      continue
    end
    i = 1
    while is_zero_entry(x, 1, i)
      i += 1
    end
    mul!(x, inv(x[1, i]))
    ta = searchsorted(piv, i)
    @assert length(ta) == 0
    if first(ta) == 1
      pushfirst!(piv, i)
      s = vcat(x, view(s, 1:r, :))
    elseif first(ta) <= r
      insert!(piv, first(ta), i)
      _s = view(s, 1:first(ta)-1, :)
      reduce_mod!(_s, x)
      s = vcat(_s, x, view(s, first(ta):r, :))
    else
      push!(piv, i)
      _s = view(s, 1:r, :)
      reduce_mod!(_s, x)
      s = vcat(_s, x)
    end
    if nrows(s) > limit
      return C
    end
    xx = x * inv(content(x))
    sx = maximum(nbits, xx)
    for g = gens(G)
      push!(new, (xx, g, sx))
    end
    r += 1
#      _x = rref(s)
#      @assert _x[1] == r
#      @assert _x[2] == s
  end
  v = [C.M(collect(view(s, i, :))) for i=1:r]
  ss, mss = sub(C.M, v)
  mss = _simplify(mss)
  v = map(mss, gens(ss))
  #preimage can handle vectors of elements, this means only one rref
  #possibly eventually homs should use the solve_ctx?
  #similarly, action can handle vectors, removing identical group theory
  return gmodule(group(C), [hom(ss, ss, preimage(mss, action(C, g, v))) for g = gens(group(C))])
end

#b is an endomorphism of M
# the kernel of a factor of the minpoly will give a submodule
# - if the minpoly splits non-trivially
# the quotient by the kernel would also work, but as this is
# also used with spinning, we need submodules.
function split_via_endo(b, M::GModule{<:Any, <:AbstractAlgebra.FPModule{QQFieldElem}})
  H = []
  iszero(b) && error("b is zero")
  f = minpoly(b)
  lf = factor(f)
  if length(lf) == 1
    return []
  end
  for (p, k) = lf
    x = (p^k)(b)
    h = hom(M, M, hom(M.M, M.M, matrix(x)))
    k, mk = kernel(h)
    @assert dim(k) > 0
#    q, mq = quo(M, mk.module_map) #TODO: for uncondense we cannot use
          #quotients
          #is quotient equivalent to kernel of other factor?
#    @assert dim(q) > 0
    append!(H, split_into_homogenous(k))
#    append!(H, split_into_homogenous(q))
#    break
  end
  return H
end

function Oscar.lll(M::QQMatrix)
  m, d = integral_split(M, ZZ)
  return lll(m)*QQ(1, d)
end

function Oscar.lll_basis(M::Hecke.AlgAssAbsOrd{MatAlgebra{QQFieldElem, QQMatrix}, ZZRing})
  A = algebra(M)
  b = basis(M, A)
  m = matrix(QQ, transpose(hcat([coefficients(x) for x = b]...)))
  m = lll(m)
  return [M(A(m[i, :])) for i=1:nrows(m)]
end

function Oscar.lll_basis(M::Hecke.AlgAssAbsOrd{ZZRing, MatAlgebra{QQFieldElem, QQMatrix}})
  A = algebra(M)
  b = basis(M, A)
  m = matrix(QQ, transpose(hcat([coefficients(x) for x = b]...)))
  m = lll(m)
  return [M(A(m[i, :])) for i=1:nrows(m)]
end

function Oscar.lll_basis(I::Hecke.AlgAssAbsOrdIdl{MatAlgebra{QQFieldElem, QQMatrix}, ZZRing})
  M = order(I)
  A = algebra(M)
  b = basis(I) # is in A
  m = matrix(QQ, transpose(hcat([coefficients(x) for x = b]...)))
  m = lll(m)
  return [M(A(m[i, :])) for i=1:nrows(m)]
end

#tries to split the homogeous module M without changing the field
function split_homogeneous(M::GModule{<:Any, <:AbstractAlgebra.FPModule{QQFieldElem}})
  #Steel, p31: MaximalOrderBasisSearch
  #            well, Step 1
  #            need to look for more elements - but how many?
  E, mE = endo(M)
  Z_M = maximal_order(E)

  if group(M) !== nothing
    chi = character(M)
    chi, k, m = galois_representative_and_multiplicity(chi)

    if m == 1
#      @show :is_irr
      return [hom(M, M, id_hom(M.M))]
    end
  end

  S = []
  B = lll_basis(Z_M)
  @assert all(!iszero, B)
  seen = Set{elem_type(E)}()
 
  cnt = 0
  for b = elem_gen_ctx(B)
    cnt += 1
    if cnt > length(B)^2 #TODO: have a sane stopping condition
      break
    end
    x = b.elem_in_algebra
    if x in seen
      continue
    end
    push!(seen, x)
    z = split_via_endo(x, M)
    if length(z) > 0
      return z
      append!(S, z)
    end
  end

  return S

  p = 20
  for i=1:10
    p = next_prime(p)

    P = prime_ideals_over(Z_M, p)[1]
    B = lll_basis(P)
    for b = elem_gen_ctx(B)
      x = b.elem_in_algebra
      if x in seen
        continue
      end
      push!(seen, x)
      z = split_via_endo(x, M)
      if length(z) > 0
        return z
        append!(S, z)
      end
    end
  end
  return S
end

function split_homogeneous2(M::GModule{<:Any, <:AbstractAlgebra.FPModule{QQFieldElem}})
  #Steel, p36: SplitHomogeneousByMinimalField(M)
  #            ... up to Step 2
  #TODO:-use this to reduce to the m=1 part - ignoring the Schur stuff
  #      (don't know if possible: rho_F is abs. irr, but as field is too large,
  #      the restriction of scalars has multiplicity again)
  #     -use the Amitsur paper instead of Fieker to reduce field? (Does not need
  #      normality, possibly)
  #TODO:-does e in Step 2 exists for A-modules (as opposed to G-modules)? Tommy
  #      spontaneously said yes...
  #XXX:  Step 3 exists - but requires the field to be normal, hence sucks in general
  #      the search for e need to be more intelligent and possibly try  more elements
  #      to find small normal closure
  #     -or try smaller degree to get some split? Any reduction in m helps
  #     -use characters to see what we do not want to split: the same module
  #      might be in a mult 8 or mult 2 component...
  #     -also do the conic case (for ms = 2)
  #     -for m=1, s=3: Jessica Cologna reduces to conics as well
  #
  #TODO:-make sure the interaction between Hecke (Tommy) and Oscar (Claus)
  #      is efficient and uses appropriate caching...
  #     -finally get some examples going and see where the infrastructure can be 
  #      improved
  #     - e.g. write & use factored_minpoly, factored_charpoly ...
  E, mE = endo(M)
  Z_M = maximal_order(E)


  chi = character(M)
  #homogeouns can have 2 problems
  # - schur index
  # - multiplicity
  # eg X is rational, but has schur index 2 => there is one rep above over
  # an extension
  # or 
  #   X is rational, but char field is quadratic, there are 2 reps over
  # the extension. 
  # in both basis is mult. 2
  # they can happen together...
  chi, k, m = galois_representative_and_multiplicity(chi)

  #we should have k*(m*s)^2 = dim(E), so
  #s = sqrt(dim(E)/k)/m
#  @show s = divexact(root(divexact(dim(E), k), 2), m)
  s = schur_index(chi)
  m = divexact(m, s)

  first = true
  local best_f::QQPolyRingElem
  local best_i::elem_type(E)
  local last_cnt::Int

  cnt = 0

  for _i=elem_gen_ctx(lll_basis(Z_M))
    cnt += 1
    if !first && cnt - last_cnt > degree(Z_M)
      break
    end
    i = _i.elem_in_algebra
    f = minpoly(i)
    lf = factor(f)
    if length(lf) == 1 && degree(f) == s*m*k &&  haskey(lf.fac, f)
      if first
        best_f = f
        best_i = i
        first = false
        last_cnt = cnt
      else
        if length(string(f)) < length(best_f) ||
          discriminant(f) < discriminant(best_f)
          best_f = f
          best_i = i
          last_cnt = cnt
        end
      end
    end
  end
  if !first
    Ka, a = number_field(best_f)
    k = eigenspace(mE(best_i), a)
#      Mf = extension_of_scalars(M, Ka)
#      h = hom(Mf, Mf, hom(Mf.M, Mf.M, map_entries(Ka, i.matrix) - a*identity_matrix(Ka, dim(M))))
#      k = kernel(h)
    #k is abs. irr. but field is too large
    return k
  end
  return :nothing_found
end

#XXX: is this the correct interface?
function Oscar.eigenspace(f::Oscar.GModuleHom{<:Any, T, T}, a::AbsSimpleNumFieldElem) where T <: AbstractAlgebra.FPModule{QQFieldElem}
  Ka = parent(a)
  M = domain(f)
  @req codomain(f) == M "1st argument must be an endomorphism"
  Mf = extension_of_scalars(M, Ka)
  h = hom(Mf, Mf, hom(Mf.M, Mf.M, map_entries(Ka, f.module_map.matrix) - a*identity_matrix(Ka, dim(M))))
  return kernel(h)[1]
end

"""
  Given a vector of elements, return an iterator for small random
  combinations. The first elments returned are the input, then sums of 
  random pairs, triples, ...
"""
mutable struct elem_gen_ctx{T <: Any} 
  bas::Vector{T}
  function elem_gen_ctx(a::Vector{T}) where T
    return new{T}(a)
  end
end

function Base.iterate(a::elem_gen_ctx)
  return a.bas[1], (1, length(a.bas)-1)
end

function Base.iterate(a::elem_gen_ctx, st::Tuple{Int, Int})
  lev, len = st
  if lev == 1 && len > 0
    len -= 1
    return a.bas[end-len], (lev, len)
  end
  if len <= 1
    lev += 1
    if lev > length(a.bas)
      return nothing
    end
    len = Int(min(2^30, div(binomial(ZZ(length(a.bas)), ZZ(lev)), 2)))
    # div(, 2) to avoid duplicates by birthday paradox
  end
  return sum(rand(a.bas) for i = 1:lev), (lev, len-1)
end
last_M = []
function split_into_homogenous(M::GModule{<:Any, <:AbstractAlgebra.FPModule{QQFieldElem}})
  #Steel, p28: HomogeneousComponents(M)
  #Careful: the list in the end contains homogenous components - but

  #         with lots of repetition
  #TODO: write and use CentreOfEndomorphismRing
  #      have a sane overall strategy
  #TODO: center_of_endo for sub and quo modules
  #      irreducibles for abelian groups
  if is_regular_gmodule(M) || has_attribute(M, :center_endo)
    C, mC = Oscar.GModuleFromGap.center_of_endo(M)
  else
    E, mE = endo(M)
    C, mC = center(E)
    mC = mC*mE
  end
  H = []
  for b = elem_gen_ctx(basis(C))
    f = minpoly(b)
    lf = factor(f)
    if length(lf) == 1
      if degree(f) == dim(C)
        return [hom(M, M, id_hom(M.M))]
      end
      continue
    end
#    @show lf, length(basis(C))
    if degree(f) < length(basis(C))
#      @show :try_again
      continue
    end
    #Q: under this conditions is the following true:
    # - k is irreducible over Q
    # - it will split (one way or other) over the field created by p
    for (p, k) = lf
      x = (p^k)(b)
      h = mC(x)
      k, mk = kernel(h)
      set_attribute!(k, :get_center_endo => x->restrict_endo(mC, mk)[2])
      if has_attribute(M, :endo)
        set_attribute!(k, :get_endo => x -> restrict_endo(get_attribute(M, :endo), mk)[2])
      end
#      if degree(p) > 1 #can we use it to split over this extension?
#        set_attribute!(k, :split =>(p, b, mC))
#      end
#      q, mq = quo(M, mk.module_map)
#      @assert dim(q) > 0 && dim(k) > 0
#      _, mF = restrict_endo(mC, mk)
#      set_attribute!(k, :center_endo => mF)
      push!(H, mk)
#      append!(H, split_into_homogenous(k))
#      _, mF = proj_endo(mC, hom(M, q, mq; check = false))
#      set_attribute!(q, :center_endo => mF)
#      append!(H, split_into_homogenous(q))
#      return H
    end
    return H
  end
  return H
end

#= Example
  G = transitive_group(11, 6)
  X = character_table(G)
  c = Oscar.GModuleFromGap.gmodule_new(X[6])

  depending on strategy:
    there are elements in the center where the minpoly has too small degree
    then the resulting modules are not irreducible
  if we search for elements with maximal degree minpoly, then indeed the 
      resulting kernels are irr.
  but the search might be longer...


G = transitive_group(11, 6)
s = subgroup_reps(G)

Xs = character_table(s[end-2])
Oscar.MeatAxe.gmodule_new(Xs[4])
v = ans[1]
V = induce(v, embedding(group(v), G))[1]

[ Oscar.MeatAxe.dim_condensed(V, u) for u = s if order(u) < 20]


c = Oscar.MeatAxe.condense(V, s[11])
Oscar.MeatAxe.trace_condensed(V, c[2], c[3], G[1])
Oscar.MeatAxe.trace_condensed(V, c[2], c[3], G[2])
z = [rand(G) for i=1:10]

[Oscar.MeatAxe.trace_condensed(X[3]+X[4], s[11], g) for g = z]
vv = Oscar.MeatAxe.split_into_homogenous(c[1])
Oscar.MeatAxe.trace_condensed(V, vv[1].module_map*c[2], c[3], G[1])
Oscar.MeatAxe.spin2(V, map(vv[1].module_map*c[2], gens(domain(vv[1]).M)))
coordinates(character(ans))


=#  


end #module
