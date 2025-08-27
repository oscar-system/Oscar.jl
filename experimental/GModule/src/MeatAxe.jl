module MeatAxe
using Oscar
using Oscar.GModuleFromGap
import Oscar: GAPGroup
import Oscar.GModuleFromGap: is_regular_gmodule


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

function _do_one!(res, t, M)
  _t = split_into_homogenous(M)
  for _x = _t
    x = domain(_x)
    _c = map(iszero, coordinates(character(x)))
    if t .& _c == t
      @show :bad, "module known", _c
      continue
    end
    @show "found", _c
    t .&= map(iszero, coordinates(character(x)))
    push!(res, x)
    if !any(t)
      return true
    end
  end
  return false
end

function _try_tensor_products!(res::Vector{GModule}, t::Vector{Bool}; limit::Int = 100)
  for i = 1:length(res)
    for j = 2:length(res)
      if dim(res[i])*dim(res[j]) <= limit
        c = coordinates(tensor_product(character(res[i]), character(res[j])))
        cc = map(iszero, c) .& t
        if cc != t
          @show :bingo
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
    @show :bingo
    M = permutation_gmodule(G, U, QQ)
    _do_one!(res, t, M) && return true
  end
  return false
end


function gmodule_new(G::Group; limit::Int = 100, t::Union{Nothing, Vector{Bool}} = nothing)
  if is_abelian(G)
    return Oscar.GModuleFromGap._irred_abelian(G)
  end
  X = character_table(G)
  s = [representative(x) for x = low_index_subgroup_classes(G, limit)]
  s = [x for x = s if order(x) < order(G)]
  sort!(s, lt = (a,b) -> isless(order(b), order(a)))
  if t !== nothing
    @assert length(t) == length(X)
  else
    t = [true for x = X]
  end
  res = GModule[]
  #permutation gmodules are good as the endomorphism ring
  #is "for free", the endo is used for the splitting
  #poss. the limits should be adjusted as this is easier...
  for i = s
    _try_perm_character!(res, t, G, i) && return res
  end
  @show "after permutation modules", t

  _try_tensor_products!(res, t; limit) && return res
  _try_squares!(res, t; limit) && return res

  all_K = [x for x = subgroup_reps(G) if order(x) <= 20] #how to do this?

  for i=s
    @show "order $(order(i))"
    XX = character_table(i)
    for chi = XX
      #careful with multiplicity
      deg = degree_of_character_field(chi)*degree(chi)*index(G, i) 
      if deg > limit^2
        @show :too_large, degree_of_character_field(chi), degree(chi), index(G, i)
        continue
      end
      chi_G = induce(chi, G)
      cc = map(iszero, coordinates(chi_G)) .& t
      if cc != t
        @show "could help"
        v = gmodule_new(chi)
        @assert length(v) == 1
        v = v[1]
        if deg <= limit
          @show "direct induce"
          for j=1:length(t)
            if cc[j] != t[j] 
              @show :induce
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
              @show t 
              @show [is_zero(scalar_product(restrict(X[i], K), trivial_character(K))) for i = 1:length(t)]
              if all(i->!t[i] || is_zero(scalar_product(restrict(X[i], K), trivial_character(K))), 1:length(t))
                @show "condense would give no info"
                continue
              end
              push!(good_K, (K, dim_C))
            end
          end
          sort!(good_K, lt = (a,b) -> a[2] < b[2])
          randG = [rand(G) for i=1:20]
          @show good_K
          for _K = good_K[1:1]
            K = _K[1]
            if !have_V
              V = induce(v, embedding(i, G))[1]
              have_V = true
            end
            @show "condense", V, K
            c, mc = condense(V, K)
            @assert dim(c) <= limit
            sp = split_into_homogenous(c)
            @show "split done"
            #TODO: use the traces in some form to decide which element in
            #      sp needs spinning - generically, there should be only one
            #      this is via Allan's trace vectors
            #      not sure what happens if split_into_homo has multiplicities
            #TODO: trace vectors: Allan adds more "generators" to the condensed
            #      module to have traces and to avoid problems from
            #      condensing Q[G] wrongly. Not sure if this is the correct
            #      way...
            for _t = sp
              d = spin2(V, [mc(_t(x).data) for x = gens(domain(_t))]; limit)
              if dim(d) > limit || d == V
                @show "spin too large", d, limit
              end
              dim(d) > limit && continue
              @assert dim(d) < dim(V)
              @show dim(d), dim(V), coordinates(character(d))
              if dim(d) <= limit
                _do_one!(res, t, d) && return res
              end
              if !any(cc .& t)
                @show "OK, leaving this level"
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

function dim_condensed(C::GModule, K::Oscar.GAPGroup)
  #Allan, Cor. 3.5.5.
  return Int(scalar_product(restrict(character(C), K), trivial_character(K)))
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
  return QQ(1, order(domain(iKG)))*sum([action(C, iKG(k)) for k = domain(iKG)])
end

function condense(C::GModule, K::Oscar.GAPGroup)
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
  x = matrix(ms)
  n = numerator(x)
  @show maximum(nbits, n)
  n = saturate(n)
  n = lll(n)
  @show maximum(nbits, n)
  ms = hom(s, C.M, matrix(QQ, n))
  #action and preimage allows vectors
  return gmodule(nothing, [hom(s, s, preimage(ms, map(phi, action(C, g, map(ms, gens(s)))))) for g = gens(G)]), ms
#  return gmodule(nothing, [hom(s, s, preimage(ms, idempotent(C, iKG, action(C, g, map(ms, gens(s)))))) for g = gens(G)]), ms
end

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

  done = false
  x = zero_matrix(QQ, 1, ncols(s))
  while !done
    done = true
    for m = 1:r, g = gens(group(C))
      #(bad) strategy: for each basis element (row in s) test if the image
      #is also in.
      #non-trivial part: use reduce_mod! as a membership test 
      #  ... and insert the reduced vector into the correct row
      #  of s - to not redo the rref
      #possibly not test the same row/ gen combination over and over?
      #might be tricky as elements change due to reduction to the top
      #maybe not do this? (the reduce_mod! below). For the spin we don't need
      #it, for the upper reduce_mod! (membership) neither ...
      mul!(x, view(s, m:m, :), matrix(action(C, g)))  #if action is matrix free...
      reduce_mod!(x, view(s, 1:r, :))
      if is_zero(x)
        continue
      end
      done = false
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
      r += 1
#      _x = rref(s)
#      @show x[2] - s
#      @assert _x[1] == r
#      @assert _x[2] == s
      break
    end
  end
  v = [C.M(collect(view(s, i, :))) for i=1:r]
  ss, mss = sub(C.M, v)
  #preimage can handle vectors of elements, this means only one rref
  #possibly eventually homs should use the solve_ctx?
  #similarly, action can handle vectors, removing identical group theory
  return gmodule(group(C), [hom(ss, ss, preimage(mss, action(C, g, v))) for g = gens(group(C))])
#  return gmodule(group(C), [hom(ss, ss, [preimage(mss, action(C, g, m)) for m = v]) for g = gens(group(C))])
end




function split_via_endo(b, M::GModule{<:Any, <:AbstractAlgebra.FPModule{QQFieldElem}})
  H = []
  iszero(b) && error("b is zero")
  f = minpoly(b)
  @show lf = factor(f)
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
      @show :is_irr
      return [M]
    end
  end

  S = []
  B = lll_basis(Z_M)
  @assert all(!iszero, B)
  seen = Set{elem_type(E)}()
  
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

  for _i=elem_gen_ctx(lll_basis(Z_M))
    i = _i.elem_in_algebra
    f = minpoly(i)
    lf = factor(f)
    if length(lf) == 1 && degree(f) == s*m*k &&  haskey(lf.fac, f)
      if first
        best_f = f
        best_i = i
        first = false
      else
        if length(string(f)) < length(best_f) ||
          discriminant(f) < discriminant(best_f)
          best_f = f
          best_i = i
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

function split_into_homogenous(M::GModule{<:Any, <:AbstractAlgebra.FPModule{QQFieldElem}})
  #Steel, p28: HomogeneousComponents(M)
  #Careful: the list in the end contains homogenous components - but
  randG = [rand(G) for i=1:20]
  #         with lots of repetition
  #TODO: write and use CentreOfEndomorphismRing
  #      have a sane overall strategy
  #TODO: center_of_endo for sub and quo modules
  #      irreducibles for abelian groups
  if is_regular_gmodule(M) || has_attribute(M, :center_endo)
    C, mC = center_of_endo(M)
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
    @show lf, length(basis(C))
    if degree(f) < length(basis(C))
      @show :try_again
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

=#  


end #module
