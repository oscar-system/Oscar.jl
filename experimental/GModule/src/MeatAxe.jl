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
      tt = map(character, t)
      return [t[i] for i=1:length(t) if !is_zero(scalar_product(tt[i], chi))]
    end
  end
  error("no plan yet")
end

function _do_one!(res, t, M)
  _t = split_into_homogenous(M)
  for x = _t
    _c = map(iszero, coordinates(character(x)))
    if t .& _c == t
      @show :bad
      continue
    end
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

  for i=s
    @show "order $(order(i))"
    XX = character_table(i)
    for chi = XX
      #careful with multiplicity
      if degree_of_character_field(chi)*degree(chi)*index(G, i) > limit
        @show :too_large, degree_of_character_field(chi), degree(chi), index(G, i)
        continue
      end
      cc = map(iszero, coordinates(induce(chi, G))) .& t
      if cc != t
        @show "should help"
        v = gmodule_new(chi)
        @assert length(v) == 1
        v = v[1]
        for j=1:length(t)
          if cc[j] != t[j] 
            @show :induce
            V = induce(v, embedding(i, G))[1]
            _do_one!(res, t, V) && return res
          else
            @show "bad"
          end
        end
      end
    end
  end
  return res
  error("no plan yet")
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
    q, mq = quo(M, mk.module_map)
    @assert dim(q) > 0
    append!(H, split_into_homogenous(k))
    append!(H, split_into_homogenous(q))
    break
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
  chi = character(M)
  chi, k, m = galois_representative_and_multiplicity(chi)

  if m == 1
    @show :is_irr
    return [M]
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
        return [M]
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
      push!(H, k)
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
