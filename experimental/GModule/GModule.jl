module GModuleFromGap
using Oscar
using Hecke

import Oscar:gmodule

import AbstractAlgebra: Group, Module
import Base: parent

function GAP.gap_to_julia(::Type{QabElem}, a::GAP.GapObj) #which should be a Cyclotomic
  c = GAP.Globals.Conductor(a)
  E = abelian_closure(QQ)[2](c)
  z = parent(E)(0)
  co = GAP.Globals.CoeffsCyc(a, c)
  for i=1:c
    if !iszero(co[i])
      z += fmpq(co[i])*E^(i-1)
    end
  end
  return z
end

(::QabField)(a::GAP.GapObj) = GAP.gap_to_julia(QabElem, a)

function irreducible_modules(G::Oscar.GAPGroup)
  im = GAP.Globals.IrreducibleRepresentations(G.X)
  IM = GModule[] 
  K = abelian_closure(QQ)[1]
  for m in im
    z = map(x->matrix(map(y->map(K, y), m(x.X))), gens(G))
    F = free_module(K, nrows(z[1]))
    push!(IM, gmodule(G, map(x->hom(F, F, x), z)))
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

function gmodule(::typeof(CyclotomicField), C::GModule)
  @assert isa(base_ring(C), QabField)
  d = dim(C)
  l = 1
  for g = C.ac
    l = lcm(l, lcm(collect(map_entries(x->Hecke.iscyclotomic_type(parent(x.data))[2], mat(g)))))
  end
  K = cyclotomic_field(l, cached = false)[1]
  F = free_module(K, dim(C))
  return gmodule(group(C), [hom(F, F, map_entries(x->K(x.data), mat(x))) for x = C.ac])
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
  return GModule(group(C), [hom(F, F, map_entries(phi, mat(x))) for x = C.ac])
end

function gmodule(::FlintRationalField, C::GModule{<:Any, Generic.FreeModule{nf_elem}})
  F = free_module(QQ, dim(C)*degree(base_ring(C)))
  return GModule(group(C), [hom(F, F, hvcat(dim(C), [representation_matrix(x) for x = transpose(mat(y))]...)) for y = C.ac])
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
  z = AbstractAlgebra.get_special(C, :Gap)
  if z !== nothing
    return z
  end
  z = GAP.Globals.GModuleByMats(GAP.julia_to_gap([GAP.julia_to_gap(map(h, Matrix(mat(x)))) for x = C.ac]), codomain(h))
  AbstractAlgebra.set_special(C, :Gap=>z)
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
  b = map(x->matrix(map(y->preimage(h, y), Oscar.GAP.gap_to_julia(Matrix{Any}, x))), hb)
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
  F = free_module(ZZ, dim(C))
  return gmodule(group(C), [hom(F, F, integral_split(mat(x), ZZ)[1]) for x = D.ac])
end

function Base.transpose(C::GModule{<:Any, <:Generic.FreeModule})
  return gmodule(group(C), [hom(C.M, C.M, transpose(mat(x))) for x = C.ac])
end

function Oscar.dual(C::GModule{<:Any, <:Generic.FreeModule})
  D = gmodule(group(C), [hom(C.M, C.M, inv(transpose(mat(x)))) for x = C.ac])
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
  r, k = kernel(transpose(vcat([matrix(base_ring(C), 1, dim(C)^2, vec(collect(x-transpose(x)))) for x = h]...)))
  return [sum(h[i]*k[i, j] for i=1:length(h)) for j=1:r]
end

function simplify(C::GModule{<:Any, <:Generic.FreeModule{fmpz}})

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

end #module GModuleFromGap

using .GModuleFromGap

export irreducible_modules, isabsolutely_irreducible, isdecomposable

