module GModuleFromGap
using Oscar
using Hecke

import Oscar:gmodule, GAPWrap

import AbstractAlgebra: Group, Module
import Base: parent

function irreducible_modules(G::Oscar.GAPGroup)
  im = GAP.Globals.IrreducibleRepresentations(G.X)
  IM = GModule[] 
  K = abelian_closure(QQ)[1]
  for m in im
    z = map(x->matrix(map(y->map(K, y), m(x.X))), gens(G))
    if ngens(G) == 0
      F = free_module(K, 0)
      zz = typeof(hom(F, F, elem_type(F)[]))[]
    else
      F = free_module(K, nrows(z[1]))
      zz = map(x->hom(F, F, x), z)
    end
    push!(IM, gmodule(F, G, zz))
  end
  return IM
end

function minimize(::typeof(CyclotomicField), a::AbstractArray{nf_elem})
  fl, c = Hecke.iscyclotomic_type(parent(a[1]))
  @assert all(x->parent(x) == parent(a[1]), a)
  @assert fl
  for p = keys(factor(c).fac)
    while c % p == 0
      K, _ = cyclotomic_field(Int(div(c, p)), cached = false)
      b = similar(a)
      OK = true
      for x = eachindex(a)
        y = Hecke.force_coerce_cyclo(K, a[x], Val{false})
        if y == false
          OK = false
        else
          b[x] = y
        end
      end
      if OK
        a = b
        c = div(c, p)
      else
        break
      end
    end
  end
  return a
end

function minimize(::typeof(CyclotomicField), a::MatElem{nf_elem})
  return matrix(minimize(CyclotomicField, a.entries))
end

function minimize(::typeof(CyclotomicField), a::nf_elem)
  return minimize(CyclotomicField, [a])[1]
end

function Oscar.conductor(a::nf_elem)
  return conductor(parent(minimize(CyclotomicField, a)))
end

function Oscar.conductor(a::QabElem)
  return conductor(data(a))
end

function irreducible_modules(::Type{AnticNumberField}, G::Oscar.GAPGroup)
  z = irreducible_modules(G)
  Z = GModule[]
  for m in z
    a = gmodule(CyclotomicField, m)
    k, mk = Hecke.subfield(base_ring(a), vec(collect(vcat(map(mat, a.ac)...))))
    if k != base_ring(a)
      F = free_module(k, dim(m))
      push!(Z, gmodule(group(m), [hom(F, F, map_entries(inv(mk), mat(x))) for x = a.ac]))
    else
      push!(Z, a)
    end
  end
  return Z
end

function irreducible_modules(::typeof(CyclotomicField), G::Oscar.GAPGroup)
  z = irreducible_modules(G)
  return [gmodule(CyclotomicField, m) for m in z]
end

function irreducible_modules(::FlintRationalField, G::Oscar.GAPGroup)
  z = irreducible_modules(CyclotomicField, G)
  return [gmodule(QQ, m) for m in z]
end

function irreducible_modules(::FlintIntegerRing, G::Oscar.GAPGroup)
  z = irreducible_modules(QQ, G)
  return [gmodule(ZZ, m) for m in z]
end

function gmodule(::typeof(CyclotomicField), C::GModule)
  @assert isa(base_ring(C), QabField)
  d = dim(C)
  l = 1
  for g = C.ac
    l = lcm(l, lcm(collect(map_entries(x->Hecke.iscyclotomic_type(parent(x.data))[2], mat(g)))))
  end
  K = cyclotomic_field(base_ring(C), l)[1]
  F = free_module(K, dim(C))
  if d == 0 
    h = hom(F, F, elem_type(F)[])
    return gmodule(F, group(C), typeof(h)[hom(F, F, map_entries(x->K(x.data), mat(x))) for x = C.ac])
  end
  return gmodule(F, group(C), [hom(F, F, map_entries(x->K(x.data), mat(x))) for x = C.ac])
end

import Base: ^
function ^(C::GModule{<:Any, Generic.FreeModule{nf_elem}}, phi::Map{AnticNumberField, AnticNumberField})
  F = free_module(codomain(phi), dim(C))
  return GModule(group(C), [hom(F, F, map_entries(phi, mat(x))) for x = C.ac])
end

function ^(C::GModule{<:Any, T}, h::Map{S, S}) where T <: S where S
  return GModule(group(C), [inv(h)*x*h for x = C.ac])
end

function ^(C::GModule{<:Any, Generic.FreeModule{QabElem}}, phi::Map{QabField, QabField})
  F = free_module(codomain(phi), dim(C))
  return GModule(F, group(C), [hom(F, F, map_entries(phi, mat(x))) for x = C.ac])
end

function gmodule(::FlintRationalField, C::GModule{<:Any, Generic.FreeModule{nf_elem}})
  F = free_module(QQ, dim(C)*degree(base_ring(C)))
  return GModule(F, group(C), [hom(F, F, hvcat(dim(C), [representation_matrix(x) for x = transpose(mat(y))]...)) for y = C.ac])
end

function gmodule(k::Nemo.GaloisField, C::GModule{<:Any, Generic.FreeModule{fmpq}})
  F = free_module(k, dim(C))
  return GModule(group(C), [hom(F, F, map_entries(k, mat(x))) for x=C.ac])
end

function gmodule(mk::Map{AnticNumberField, <:FinField}, C::GModule{<:Any, Generic.FreeModule{nf_elem}})
  k = codomain(mk)
  @assert domain(mk) == base_ring(C)
  F = free_module(k, dim(C))
  return GModule(group(C), [hom(F, F, map_entries(mk, mat(x))) for x=C.ac])
end

function Hecke.modular_proj(C::GModule{T, Generic.FreeModule{nf_elem}}, me::Hecke.modular_env) where T
  R = []
  z = map(x->(Hecke.modular_proj(x.matrix, me)), C.ac)
  for i=1:length(z[1])
    F = free_module(base_ring(z[1][i]), dim(C))
    @assert all(j->base_ring(z[j][i]) == base_ring(z[1][i]), 1:length(z))
    push!(R, GModule(group(C), [hom(F, F, t[i]) for t = z]))
    @assert all(i->base_ring(mat(R[end].ac[i])) == base_ring(R[end]), 1:length(R[end].ac))
  end
  return R
end

function Gap(C::GModule{<:Any, <:Generic.FreeModule{<:FinFieldElem}}, h=Oscar.ring_iso_oscar_gap(base_ring(C)))
  z = get_attribute(C, :Gap)
  if z !== nothing
    return z
  end
  z = GAP.Globals.GModuleByMats(GAP.julia_to_gap([GAP.julia_to_gap(map(h, Matrix(mat(x)))) for x = C.ac]), codomain(h))
  set_attribute!(C, :Gap=>z)
  return z
end

function Oscar.isirreducible(C::GModule{<:Any, <:Generic.FreeModule{<:FinFieldElem}})
  G = Gap(C)
  return GAP.Globals.MTX.IsIrreducible(G)
end

function isabsolutely_irreducible(C::GModule{<:Any, <:Generic.FreeModule{<:FinFieldElem}})
  G = Gap(C)
  return GAP.Globals.MTX.IsAbsolutelyIrreducible(G)
end

function isdecomposable(C::GModule{<:Any, <:Generic.FreeModule{<:FinFieldElem}})
  G = Gap(C)
  return !GAP.Globals.MTX.IsIndecomposable(G)
end

function Oscar.hom(C::T, D::T) where T <: GModule{<:Any, <:Generic.FreeModule{<:FieldElem}}
  b = hom_base(C, D)
  H, mH = hom(C.M, D.M)
  s, ms = sub(H, [H(vec(collect(x))) for x = b])
  return GModule(group(C), [hom(s, s, [preimage(ms, H(vec(collect(inv(mat(C.ac[i]))*g*mat(D.ac[i]))))) for g = b]) for i=1:ngens(group(C))]), ms, mH
end

function Oscar.hom(F::Generic.FreeModule{T}, G::Generic.FreeModule{T}) where T
  k = base_ring(F)
  @assert base_ring(G) == k
  H = free_module(k, dim(F)*dim(G))
  return H, MapFromFunc(x->hom(F, G, matrix(k, dim(F), dim(G), vec(collect(x.v)))), y->H(vec(collect(transpose(mat(y))))), H, Hecke.MapParent(F, G, "homomorphisms"))
end

function hom_base(C::T, D::T) where T <: GModule{<:Any, <:Generic.FreeModule{<:FinFieldElem}}
  @assert base_ring(C) == base_ring(D)
  h = Oscar.ring_iso_oscar_gap(base_ring(C))
  hb = GAP.Globals.MTX.BasisModuleHomomorphisms(Gap(C, h), Gap(D, h))
  n = length(hb)
  b = [matrix([preimage(h, x[i, j]) for i in 1:GAPWrap.NrRows(x), j in 1:GAPWrap.NrCols(x)]) for x in hb]
#  b = map(x->matrix(map(y->preimage(h, y), Matrix{Any}(x))), hb)
#  @show [mat(C.ac[i])*b[1] == b[1]*mat(D.ac[i]) for i=1:length(C.ac)]
  return b
end

"""
  C*T[i] = T[i]*D
on return
"""
function hom_base(C::_T, D::_T) where _T <: GModule{<:Any, <:Generic.FreeModule{nf_elem}}
  @assert base_ring(C) == base_ring(D)

  p = Hecke.p_start
  p = 2^10
  p = 127
  m_in = map(mat, C.ac)
  m_out = map(mat, D.ac)
  local T
  pp = fmpz(1)
  k = base_ring(C)
  @assert base_ring(m_in[1]) == k
  @assert base_ring(m_in[1]) == k
  while true
    p = next_prime(p)
    me = modular_init(k, p)
    z1 = Hecke.modular_proj(C, me)
    if C === D
      z2 = z1
    else
      z2 = Hecke.modular_proj(D, me)
    end
    t = []
    for i=1:length(z1)
      push!(t, hom_base(z1[i], z2[i]))
    end
    tt = [Hecke.modular_lift([t[i][j] for i=1:length(z1)], me) for j=1:length(t[1])]
    if length(tt) == 0
      return []
    end
    @assert base_ring(tt[1]) == k
    if isone(pp)
      pp = fmpz(p)
      T = tt
    else
      T = [induce_crt(tt[i], T[i], fmpz(p), pp) for i=1:length(T)]
      @assert base_ring(T[1]) == k
      pp *= p
      S = []
      for t = T
        fl, s = induce_rational_reconstruction(t, pp)
        fl || break
        push!(S, s)
      end
      @assert base_ring(S[1]) == k
      s = S[1]
      if length(S) == length(T)
        if all(s->all(i->m_in[i]*s ==  s*m_out[i], 1:length(m_in)), S)
          return S
        end
      end
    end
  end
end

function hom_base(C::_T, D::_T) where _T <: GModule{<:Any, <:Generic.FreeModule{fmpq}}
  @assert base_ring(C) == base_ring(D)

  p = Hecke.p_start
  p = 2^10
  p = 127
  m_in = map(mat, C.ac)
  m_out = map(mat, D.ac)
  local T
  pp = fmpz(1)
  k = base_ring(C)
  @assert base_ring(m_in[1]) == k
  @assert base_ring(m_in[1]) == k
  @assert k == QQ
  while true
    p = next_prime(p)
    z1 = gmodule(GF(p), C)
    if C === D
      z2 = z1
    else
      z2 = gmodule(GF(p), D)
    end
    
    t = hom_base(z1, z2)
    tt = [lift(s)  for s=t]
    @assert base_ring(tt[1]) == ZZ
    if isone(pp)
      pp = fmpz(p)
      T = tt
    else
      T = [induce_crt(tt[i], fmpz(p), T[i], pp)[1] for i=1:length(T)]
      @assert base_ring(T[1]) == ZZ
      pp *= p
      S = []
      for t = T
        fl, s = induce_rational_reconstruction(t, pp)
        fl || break
        push!(S, s)
      end
      if nbits(pp) > 1000
        error("ndw")
      end
      if length(S) == length(T)
        if all(s->all(i->m_in[i]*s ==  s*m_out[i], 1:length(m_in)), S)
          return S
        end
      end
    end
  end
end

function gmodule(K::AnticNumberField, M::GModule{<:Any, <:Generic.FreeModule{nf_elem}})
  F = free_module(K, dim(M))
  return gmodule(F, group(M), [hom(F, F, map_entries(K, mat(x))) for x = M.ac])
end

function (K::QabField)(a::nf_elem)
  fl, f = Hecke.iscyclotomic_type(parent(a))
  @assert fl
  return QabElem(a, f)
end

function hom_base(C::_T, D::_T) where _T <: GModule{<:Any, <:Generic.FreeModule{QabElem}}
  C1 = gmodule(CyclotomicField, C)
  D1 = gmodule(CyclotomicField, D)
  fl, Cf = Hecke.iscyclotomic_type(base_ring(C1))
  @assert fl
  fl, Df = Hecke.iscyclotomic_type(base_ring(D1))
  @assert fl
  l = lcm(Cf, Df)
  K, _ = cyclotomic_field(base_ring(C), l)
  if l != Cf
    C1 = gmodule(K, C1)
  end
  if l != Df
    D1 = gmodule(K, D1)
  end
  h = hom_base(C1, D1)
  if length(h) == 0
    return h
  end
  return map(x->map_entries(base_ring(C), x), h)
end

function gmodule(::FlintRationalField, C::GModule{<:Any, <:Generic.FreeModule{fmpz}})
  F = free_module(QQ, dim(C))
  return GModule(group(C), [hom(F, F, map_entries(QQ, mat(x))) for x = C.ac])
end

function hom_base(C::_T, D::_T) where _T <: GModule{<:Any, <:Generic.FreeModule{fmpz}}

  h = hom_base(gmodule(QQ, C), gmodule(QQ, D))
  H = vcat([integral_split(matrix(QQ, 1, dim(C)^2, vec(collect(x))), ZZ)[1] for x = h]...)
  H = Hecke.saturate(H)
  return [matrix(ZZ, dim(C), dim(C), vec(collect(H[i, :]))) for i=1:nrows(H)]
end

function gmodule(::FlintIntegerRing, C::GModule{<:Any, <:Generic.FreeModule{fmpq}})
  ma = map(mat, C.ac)
  M = identity_matrix(QQ, dim(C))
  while true
    N = reduce(vcat, [M*x for x = ma])
    H = hnf(integral_split(N, ZZ)[1])[1:dim(C), :]
    if H == M
      break
    end
    M = map_entries(QQ, H)
  end
  M = inv(M)
  h = hom(C.M, C.M, M)
  D = C^h
  return gmodule(group(C), [integral_split(x, ZZ)[1] for x = action_matrices(D)])
end

function Base.transpose(C::GModule{<:Any, <:Generic.FreeModule})
  return gmodule(group(C), [transpose(x) for x = action_matrices(C)])
end

function Oscar.dual(C::GModule{<:Any, <:Generic.FreeModule})
  D = gmodule(group(C), [inv(transpose(x)) for x = action_matrices(C)])
  return D
end

#if C is abs. irr <=> hom is 1 dim => this always works
#otherwise, this gives a basis of the symmetric matrices, stabilized
#by the group - but possibly not pos. definite
#
#to always get pos. def. one needs a different algo
# - sum over group
# - approximate reynolds as in invar thy (for number fields and Q/ Z)
function invariant_forms(C::GModule{<:Any, <:Generic.FreeModule})
  D = Oscar.dual(C)
  h = hom_base(C, D)
  r, k = kernel(transpose(vcat([matrix(base_ring(C), 1, dim(C)^2, vec(x-transpose(x))) for x = h]...)))
  return [sum(h[i]*k[i, j] for i=1:length(h)) for j=1:r]
end

function Oscar.gmodule(G::Oscar.GAPGroup, v::Vector{<:MatElem})
  @assert length(v) == ngens(G)
  R = base_ring(v[1])
  @assert all(x->R == base_ring(x), v)
  @assert nrows(v[1]) == ncols(v[1])
  @assert all(x->size(v[1]) == size(x), v)
  F = free_module(R, nrows(v[1]))
  return gmodule(G, [hom(F, F, x) for x = v])
end


function Oscar.gmodule(::Type{GrpAbFinGen}, C::GModule{T, Generic.FreeModule{fmpz}}) where {T <: Oscar.GAPGroup}
  A = free_abelian_group(rank(C.M))
  return Oscar.gmodule(Group(C), [hom(A, A, mat(x)) for x = C.ac])
end

function Oscar.gmodule(::Type{GrpAbFinGen}, C::GModule{T, Generic.FreeModule{gfp_fmpz_elem}}) where {T <: Oscar.GAPGroup}
  A = abelian_group([characteristic(base_ring(C)) for i=1:rank(C.M)])
  return Oscar.gmodule(A, Group(C), [hom(A, A, map_entries(lift, mat(x))) for x = C.ac])
end

#to bypass the vec(collect(M)) which copies twice
function Base.vec(M::Generic.Mat)
  return vec(M.entries)
end

function Base.vec(M::MatElem)
  r = elem_type(base_ring(M))[]
  for j=1:ncols(M)
    for i=1:nrows(M)
      push!(r, M[i, j])
    end
  end
  return r
end

function Oscar.simplify(C::GModule{<:Any, <:Generic.FreeModule{fmpq}})
  return gmodule(QQ, Oscar.simplify(gmodule(ZZ, C))[1])
end

function action_matrices(C::GModule{<:Any, <:Generic.FreeModule})
  return map(mat, action(C))
end

function Oscar.simplify(C::GModule{<:Any, <:Generic.FreeModule{fmpz}})
 f = invariant_forms(C)[1]
 @assert all(i->det(f[1:i, 1:i])>0, 1:nrows(f))
 m = map(mat, C.ac)
 S = identity_matrix(ZZ, dim(C))
 while true
   L, T = lll_gram_with_transform(f)
   Ti = inv(T)
   n = [T*x*Ti for x = m]
   if length(string(n)) >= length(string(m))
     return C, S
   end
   S = T*S
   C = gmodule(group(C), n)
   f = invariant_forms(C)[1]
   M = n
 end
end

function Hecke.induce_crt(a::Generic.MatSpaceElem{nf_elem}, b::Generic.MatSpaceElem{nf_elem}, p::fmpz, q::fmpz)
  c = parent(a)()
  pi = invmod(p, q)
  mul!(pi, pi, p)
  pq = p*q
  z = fmpz(0)

  for i=1:nrows(a)
    for j=1:ncols(a)
      c[i,j] = Hecke.induce_inner_crt(a[i,j], b[i,j], pi, pq, z)
    end
  end
  return c
end

function Hecke.induce_rational_reconstruction(a::Generic.MatSpaceElem{nf_elem}, pg::fmpz)
  c = parent(a)()
  for i=1:nrows(a)
    for j=1:ncols(a)
      fl, c[i,j] = rational_reconstruction(a[i,j], pg)
      fl || return fl, c
    end
  end
  return true, c
end

function Hecke.induce_rational_reconstruction(a::fmpz_mat, pg::fmpz)
  c = zero_matrix(QQ, nrows(a), ncols(a))
  for i=1:nrows(a)
    for j=1:ncols(a)
      fl, n, d = rational_reconstruction(a[i,j], pg)
      fl || return fl, c
      c[i,j] = n//d
    end
  end
  return true, c
end


export irreducible_modules, isabsolutely_irreducible, isdecomposable

## Fill in some stubs for Hecke

function _to_gap(h, x::Vector)
  return GAP.Globals.GModuleByMats(GAP.julia_to_gap([GAP.julia_to_gap(map(h, Matrix(y))) for y in x]), codomain(h))
end

function _gap_matrix_to_julia(h, g)
  return matrix(domain(h), [map(y -> preimage(h, y), gg) for gg in GAP.gap_to_julia(g)])
end

function _to_julia(h, C)
  return [ matrix(domain(h), [map(y -> preimage(h, y), gg) for gg in GAP.gap_to_julia(g)]) for g in GAP.Globals.MTX.Generators(C)]
end

if isdefined(Hecke, :stub_composition_factors)
  function Hecke.stub_composition_factors(x::Vector{T}) where {T}
    F = base_ring(x[1])
    h = Oscar.ring_iso_oscar_gap(F)
    V = _to_gap(h, x)
    Vcf = GAP.Globals.MTX.CompositionFactors(V)
    res = Vector{T}[]
    for C in Vcf
      push!(res, _to_julia(h, C))
    end
    return res
  end
end

if isdefined(Hecke, :stub_basis_hom_space)
  function Hecke.stub_basis_hom_space(x::Vector, y::Vector)
    F = base_ring(x[1])
    h = Oscar.ring_iso_oscar_gap(F)
    @assert base_ring(x[1]) == base_ring(y[1])
    @assert length(x) == length(y)
    hb = GAP.Globals.MTX.BasisModuleHomomorphisms(_to_gap(h, x), _to_gap(h, y))
    hbb = [_gap_matrix_to_julia(h, g) for g in GAP.gap_to_julia(hb)]
    return hbb
  end
end

end #module GModuleFromGap

using .GModuleFromGap

export irreducible_modules, isabsolutely_irreducible, isdecomposable

module RepPc
using Oscar

Base.pairs(M::MatElem) = Base.pairs(IndexCartesian(), M)
Base.pairs(::IndexCartesian, M::MatElem) = Base.Iterators.Pairs(M, CartesianIndices(axes(M)))

function Hecke.roots(a::fq_nmod, i::Int)
  kx, x = PolynomialRing(parent(a), cached = false)
  return roots(x^i-a)
end

#=TODO
 - construct characters along the way as well?
 - compare characters rather than the hom_base
 - maybe reason from theory what reps are going to be new?
 - conjugate to smallest field?
 - allow trivial stuff
=# 

function reps(K, G::PcGroup)
  s, ms = sub(G, [gens(G)[end]])
  o = Int(order(s))
  @assert isprime(o)
  z = roots(K(1), o)
  F = free_module(K, 1)
  R = [gmodule(F, s, [hom(F, F, [r*F[1]])]) for r = z]

  for i=ngens(G)-1:-1:1
    h = G[i]
    ns, mns = sub(G, gens(G)[i:end])
    p = Int(divexact(order(ns), order(s)))
    @assert isprime(p)
    new_R = []
    #TODO: use extend below
    for r = R
      F = r.M
      rh = gmodule(group(r), [action(r, preimage(ms, x^h)) for x = gens(s)])
      l = Oscar.GModuleFromGap.hom_base(r, rh)
      @assert length(l) <= 1
      nr = []
      Y = mat(action(r, preimage(ms, h^p)))
      if length(l) == 1
        X = l[1]
        Xp = X^p
        #Brueckner: C*Xp == Y for some scalar C
        i = findfirst(x->!iszero(x), Xp)
        @assert !iszero(Y[i])
        C = divexact(Y[i], Xp[i])
        @assert C*Xp == Y
        # I think they should always be roots of one here.
        rt = roots(C, p)
        Y = r.ac
        for x = rt
          nw = gmodule(F, ns,  vcat([hom(F, F, x*X)], Y))
          push!(nr, nw)
        end
      else #need to extend dim
        n = dim(r)
        F = free_module(K, dim(r)*p)
        z = zero_matrix(K, dim(F), dim(F))
        z[(p-1)*n+1:end, 1:n] = Y
        for i=1:p-1
          z[(i-1)*n+1:i*n, i*n+1:(i+1)*n] = identity_matrix(K, n)
        end
        md = [hom(F, F, z)]
        for g = gens(s)
          z = zero_matrix(K, dim(F), dim(F))
          for j=1:p
            Y = action(r, g)
            z[(j-1)*n+1:j*n, (j-1)*n+1:j*n] = mat(Y)
            g = preimage(ms, g^h)
          end
          push!(md, hom(F, F, z))
        end
        push!(nr, gmodule(F, ns, md))
      end
      for nw = nr
        if any(x->length(Oscar.GModuleFromGap.hom_base(x, nw))>0, new_R)
          continue
        else
          push!(new_R, nw)
        end
      end
    end
    s, ms = ns, mns
    R = new_R
  end
  return R
end

function extend(C::GModule, m::Map)
  #N acts and is normal in <h, N> 
  #m injects N into <h, N>  (at least h, maybe more)
  #h = gen(domain(m), 1)
  #h has order p in <h, N>/N
  #Satz 10 in Brueckner

  F = Module(C)
  N = group(C)
  Nh = codomain(m)
  @assert ngens(N) + 1 == ngens(Nh)
  @assert all(x->m(gen(N, i)) == gen(Nh, i+1), 1:ngens(N))

  h = gen(Nh, 1)
  p = divexact(order(Nh), order(N))
  @assert isprime(p)

  F = C.M
  K = base_ring(F)
  Ch = gmodule(N, [action(C, preimage(m, m(x)^h)) for x = gens(N)])
  l = Oscar.GModuleFromGap.hom_base(C, Ch)
  @assert length(l) <= 1
  nr = []
  Y = mat(action(C, preimage(m, h^p)))
  if length(l) == 1
    X = l[1]
    Xp = X^p
    #Brueckner: C*Xp == Y for some scalar C
    i = findfirst(x->!iszero(x), Xp)
    @assert !iszero(Y[i])
    C = divexact(Y[i], Xp[i])
    @assert C*Xp == Y
    # I think they should always be roots of one here.
    rt = roots(C, p)
    Y = r.ac
    for x = rt
      nw = gmodule(F, Nh,  vcat([hom(F, F, x*X)], Y))
      push!(nr, nw)
    end
  else #need to extend dim
    n = dim(r)
    F = free_module(K, dim(r)*p)
    z = zero_matrix(K, dim(F), dim(F))
    z[(p-1)*n+1:end, 1:n] = Y
    for i=1:p-1
      z[(i-1)*n+1:i*n, i*n+1:(i+1)*n] = identity_matrix(K, n)
    end
    md = [hom(F, F, z)]
    for g = gens(s)
      z = zero_matrix(K, dim(F), dim(F))
      for j=1:p
        Y = action(r, g)
        z[(j-1)*n+1:j*n, (j-1)*n+1:j*n] = mat(Y)
        g = preimage(m, g^h)
      end
      push!(md, hom(F, F, z))
    end
    push!(nr, gmodule(F, Nh, md))
  end
  return nr
end


function brueckner(G::FPGroup)
  #=
  actually, the initial prime list should be from
  the order of the maximal abelian quotient
  However, this is intended to be "pure"...
  =#
  F = free_module(QQ, 1)
  h = hom(F, F, [F[1]])
  Q = free_group(1)
  Q = quo(Q, [Q[1]])[1]
  mQ = hom(G, Q, gens(G), [Q[1] for i=1:ngens(G)])

  C = gmodule(F, G, [h for g = gens(G)])
  CZ = gmodule(GrpAbFinGen, gmodule(ZZ, C))
  a, b = Oscar.GrpCoh.H_one_maps(CZ)
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
  =#
  da = dual(a)
  db = dual(b)
  q = quo(kernel(da)[1], image(db)[1])[1]
  t = torsion_subgroup(q)[1]
  lp = collect(keys(factor(order(t)).fac))
  allR = Dict{Int, Vector{Any}}()
  @assert ngens(Q) == 1
  allR[0] = [gmodule(F, Q, typeof(h)[])]
  for p = lp
    F = free_module(GF(p), 1)
    h = hom(F, F, [F[1]])
    allR[p] = [gmodule(F, Q, [h])]
  end
  return allR
end
Base.getindex(M::AbstractAlgebra.FPModule, i::Int) = i==0 ? zero(M) : gens(M)[i]
Oscar.gen(M::AbstractAlgebra.FPModule, i::Int) = M[i]

function dual(h::Map{GrpAbFinGen, GrpAbFinGen})
  A = domain(h)
  B = codomain(h)
  @assert isfree(A) && isfree(B)
  return hom(B, A, transpose(h.map))
end

function coimage(h::Map{GrpAbFinGen, GrpAbFinGen})
  return quo(domain(h), kernel(h)[1])
end

function lift(C::GModule, mp::Map)
  #m: G->group(C)
  #compute all(?) of H^2 that will descibe groups s.th. m can be lifted to

  G = domain(mp)
  N = group(C)
  @assert codomain(mp) == N

  _ = Oscar.GrpCoh.H_two(C)
  sc, mH2 = get_attribute(C, :H_two_symbolic_chain)
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
        m = -pro[-i]
      else
        h = mp(G[i])
        m = pro[i]
      end
      # a *(h, m) = (x, y)(h, m) = (xh, m^h + y + si(x, h))
      a = (a[1]*h, pDE[1]*m*action(C, h) + a[2] + pDE[2]*sc(a[1], h))
    end
    @assert isone(a[1])
    s += a[2]*iK[j]
    j += 1
  end
  #so kern(s) should be exactly all possible quotients that allow a 
  #projection of G. They are not all surjective. However, lets try:
  k, mk = kernel(s)
  allG = []
  z = get_attribute(C, :H_two)[1]
  for x = k
    epi = pDE[1](mk(x)) #the map
    chn = pDE[2](mk(x)) #the tail data
    #TODO: not all "chn" yield distinct groups - the factoring by the 
    #      co-boundaries is missing
    #      not all "epi" are epi, ie. sujective. The part of the thm
    #      is missing...
    # (Thm 15, part b & c)
    GG, GGinj, GGpro, GMtoGG = Oscar.GrpCoh.extension(z(chn))
    #map G[i] -> <mp(G[i]), pro[i](epi)>
    push!(allG, (GG, [hom(G, GG, gens(G), [GMtoGG(mp(G[i]), pro[i](epi)) for i=1:ngens(G)])]))
  end
  return allG
end

end #module RepPc

using .RepPc
