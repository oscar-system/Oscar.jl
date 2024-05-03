include("Cohomology.jl")
include("GaloisCohomology.jl")
include("Misc.jl")

module GModuleFromGap
using Oscar
using Hecke
import Hecke: data

#XXX: clash of names!
#   gmodule(k, C) vs gmodule_ver(k, C)
#   first does a "restriction of scalars" or blow up with the rep mat
#   second tries to conjugate down to k

import Oscar: _vec, gmodule, GAPWrap
import Oscar.GrpCoh: MultGrp, MultGrpElem

import AbstractAlgebra: Group, Module
import Base: parent

"""
    restriction_of_scalars(M::GModule, phi::Map)

Return the `S`-module obtained by restricting the scalars of `M`
from `R` to `S`, where `phi` is an embedding of `S` into `R`.

If `R` has `S`-rank `d` and `M` has rank `n` then the returned module
has rank `d*n`.

# Examples
```jldoctest
julia> G = dihedral_group(20);

julia> T = character_table(G);

julia> C = gmodule(T[8]);

julia> C = gmodule(CyclotomicField, C);

julia> h = subfields(base_ring(C), degree = 2)[1][2];

julia> restriction_of_scalars(C, h)
G-module for G acting on vector space of dimension 4 over number field

julia> restriction_of_scalars(C, QQ)
G-module for G acting on vector space of dimension 8 over QQ

```
"""
function restriction_of_scalars(M::GModule{<:Oscar.GAPGroup, <:AbstractAlgebra.FPModule{<:FieldElem}}, phi::Map)
  #works iff relative_field above works. At least for AbsSimpleNumField and
  #finite fields
  @assert codomain(phi) == base_ring(M)
  d = divexact(degree(codomain(phi)), degree(domain(phi)))
  F = free_module(domain(phi), dim(M)*d)
  _, _, rep = relative_field(phi)

  return GModule(F, group(M), [hom(F, F, hvcat(dim(M), [rep(x) for x in transpose(matrix(y))]...)) for y in M.ac])
end

function restriction_of_scalars(C::GModule{<:Any, <:AbstractAlgebra.FPModule{AbsSimpleNumFieldElem}}, ::QQField)
  F = free_module(QQ, dim(C)*degree(base_ring(C)))
  return GModule(F, group(C), [hom(F, F, hvcat(dim(C), [representation_matrix(x) for x in transpose(matrix(y))]...)) for y in C.ac])
end


"""
    extension_of_scalars(M::GModule, phi::Map)

Return the `S`-module obtained by extending the scalars of `M`
from `R` to `S`, where `phi` is an embedding of `R` into `S`.

If `M` has `R`-rank `d` then the returned module has `S`-rank `d`.

The syntax `M ⊗ phi` is supported.

# Examples

(cases that `S` is a finite field or a number field, and `R` is a subfield;
case that `S` is a number field and `R` is the ring of integers in `S`,
for example `R = ZZ` and `S = QQ`)
"""
function extension_of_scalars(M::GModule, phi::Map)
  @assert domain(phi) == base_ring(M)

  d = dim(M)
  F = free_module(codomain(phi), d)

  return GModule(F, group(M), [hom(F, F, map_entries(phi, matrix(x))) for x in M.ac])
end


"""
    can_be_defined_over(M::GModule, phi::Map)

Return `true` if there is an `S`-module `N` such that
`extension_of_scalars(N, phi)` is isomorphic with `M`,
and `false` otherwise.

`phi` is an embedding of `S` into `R` and `M` is a module over `R`.

# Examples

(case that `R` is a number field and `S` is the ring of integers in `R`)
"""
function can_be_defined_over(M::GModule, phi::Map)
  error("not yet ...")
end

function can_be_defined_over(M::GModule{<:Any, <:AbstractAlgebra.FPModule{<:FinFieldElem}}, phi::Map)
  # Only works for irreducible modules
  k = domain(phi)
  K = base_ring(M)
  d = absolute_degree(k)
  @assert absolute_degree(K) != d
  s = absolute_frobenius(K, d)
  os = divexact(absolute_degree(K), d)
  hB = hom_base(M, gmodule(M.M, Group(M),
                      [hom(M.M, M.M, map_entries(s, matrix(x))) for x = M.ac]))
  if length(hB) != 1
    length(hB) > 1 && error("Module not irreducible")
    length(hB) == 0 && return false
  end
  return true
end


"""
    can_be_defined_over_with_data(M::GModule, phi::Map)

Is similar to [`can_be_defined_over`](@ref), but the return value is
a triple `(true, N, psi)` if there is an `S`-module `N` such that
`extension_of_scalars(N, phi)` is isomorphic with `M`
and `psi` is a map from `N` to `M`.

`phi` is an embedding of `S` into `R` and `M` is a module over `R`.

# Examples

(same as for `can_be_defined_over`)
"""
function can_be_defined_over_with_data(M::GModule, phi::Map)
  error("not yet ...")
end

function can_be_defined_over_with_data(M::GModule{<:Any, <:AbstractAlgebra.FPModule{<:FinFieldElem}}, phi::Map)
  # Only works for irreducible modules
  k = domain(phi)
  K = base_ring(M)
  d = absolute_degree(k)
  @assert absolute_degree(K) != d

  s = absolute_frobenius(K, d)
  os = divexact(absolute_degree(K), d)
  hB = hom_base(M, gmodule(M.M, Group(M),
                      [hom(M.M, M.M, map_entries(s, matrix(x))) for x = M.ac]))
  if length(hB) != 1
    length(hB) > 1 && error("Module not irreducible")
    length(hB) == 0 && return false
  end

  B = hB[1]
  D = norm(B, s, os)
  lambda = D[1,1]
  @hassert :MinField 2 D == lambda*identity_matrix(K, dim(C))
  alpha = norm_equation(K, preimage(phi, lambda))
  B *= inv(alpha)
  @hassert :MinField 2 isone(norm(B, s, os))
  D = hilbert90_cyclic(B, s, os)
  Di = inv(D)
  F = free_module(k, dim(M))
  N = gmodule(F, Group(M), [hom(F, F, map_entries(x -> preimage(phi, x), Di*matrix(x)*D)) for x = C.ac])
  # TODO: Get this map working
  # psi = GModuleHom(N, M, , phi)
  psi = nothing
  return (true, N, psi)
end


"""
    descent_to(M::GModule, phi::Map)

Return an `S`-module `N` such that
`extension_of_scalars(N, phi)` is isomorphic with `M` if such an `S`-module
exists, otherwise throw an exception.

`phi` is an embedding of `S` into `R` and `M` is a module over `R`.

Use [`can_be_defined_over`](@ref) in order to check whether `M` can be
written over `S`.
Use [`can_be_defined_over_with_data`](@ref) in order to check whether `M` can be
written over `S` and to get `N` if it exists.

# Examples

(case that `R` is a number field and `S` is the ring of integers in `R`)
"""
function descent_to(M::GModule, phi::Map)
  error("not yet ...")
end


function descent_to(M::GModule{<:Any, <:AbstractAlgebra.FPModule{<:FinFieldElem}}, phi::Map)
  # Only works for irreducible modules
  k = domain(phi)
  success, N, psi = can_be_defined_over_with_data(M, phi)
  success || error("Module cannot be written over $k")
  return N
end


"""
    descent_to_minimal_degree_field(M::GModule)

Return a module `N` over a field `S` of minimal degree such that
`extension_of_scalars(N, phi)` and `extension_of_scalars(M, psi)`
are isomorphic, where `phi` and `psi` are maps from `S` and `R`,
respectively, to the compositum of `S` and `R`,
where `M` is a module over `R`.

(In general, `S` need not be a subfield of `R`.)

# Examples

(modules over finite fields or number fields)
"""
function descent_to_minimal_degree_field(C::GModule{<:Any, <:AbstractAlgebra.FPModule{fpFieldElem}})
  return C
end

function descent_to_minimal_degree_field(C::GModule{<:Any, <:AbstractAlgebra.FPModule{<:FinFieldElem}})
  #always over char field
  K = base_ring(C)
  d = 0
  while d < absolute_degree(K)-1
    d += 1
    absolute_degree(K) % d == 0 || continue
    k = GF(characteristic(K), d)
    D = gmodule_over(k, C, do_error = false)
    D === nothing || return D
  end
  return C
end

function descent_to_minimal_degree_field(C::GModule{<:Any, <:AbstractAlgebra.FPModule{AbsSimpleNumFieldElem}})
  return _minimize(C)
end


"""
    invariant_lattice_classes(M::GModule, phi::Map)

Return representatives of the equivalence classes of `G`-invariant
`S`-lattices in the `R`-module `M`,
where `phi` is an embedding of `S` into `R`.

# Examples

(case that `R` is a number field and `S` is the ring of integers in `R`)
"""
function invariant_lattice_classes(M::GModule{<:Oscar.GAPGroup, <:AbstractAlgebra.FPModule{QQFieldElem}}, phi::Map{ZZRing, QQField})
  MZ = gmodule(ZZ, M)
  return invariant_lattice_classes(MZ)
end

function invariant_lattice_classes(M::GModule{<:Oscar.GAPGroup, <:AbstractAlgebra.FPModule{ZZRingElem}})
  res = Any[(M, sub(M.M, gens(M.M))[2])]
  sres = 1
  new = true
  lp = keys(factor(order(M.G)).fac)
  while new
    new  = false
    lres = length(res)
    for X in res[sres:end]
      for p in lp
        F = free_module(GF(p), dim(M))
        pM = gmodule(M.G, [hom(F, F, map_entries(base_ring(F), x)) for x in map(matrix, action(M))])
        S = maximal_submodule_bases(pM)
        pG = p.*gens(M.M)
        for s in S
          x, mx = sub(M.M, vcat(pG, [M.M(map_entries(x->lift(ZZ, x), s[i:i, :])) for i in 1:nrows(s)]))
          r = (sub(M, mx), mx)
          if any(x->is_isomorphic(r[1], x[1]), res)
            continue
          else
            new = true
            push!(res, r)
          end
        end
      end
    end
    lres = length(res)
  end
  return res
end

Oscar.rank(M::AbstractAlgebra.Generic.Submodule{ZZRingElem}) = ngens(M)

function maximal_submodule_bases(M::GModule{<:Oscar.GAPGroup, <:AbstractAlgebra.FPModule{<:FinFieldElem}})
  C = Gap(M)
  S = GAP.Globals.MTX.BasesMaximalSubmodules(C)
  res = dense_matrix_type(base_ring(M))[]
  for s = S
    if length(s) == 0
      m = matrix(base_ring(M), 0, dim(M), [])
    else
      m = matrix(base_ring(M), s)
    end
    push!(res, m)
  end
  return res
end

function maximal_submodules(M::GModule{<:Oscar.GAPGroup, <:AbstractAlgebra.FPModule{<:FinFieldElem}})
  return [sub(M, s) for s = maximal_submodule_bases(M)]
end


function __init__()
  add_verbosity_scope(:BruecknerSQ)
  set_verbosity_level(:BruecknerSQ, 0)

  add_assertion_scope(:BruecknerSQ)
  set_assertion_level(:BruecknerSQ, 0)

  add_assertion_scope(:MinField)
  set_assertion_level(:MinField, 0)

  add_verbosity_scope(:MinField)
  set_verbosity_level(:MinField, 0)
end

function irreducible_modules(k::FinField, G::Oscar.GAPGroup)
  h = Oscar.iso_oscar_gap(k)
  hi = inv(h)
  im = GAP.Globals.IrreducibleRepresentations(GapObj(G), codomain(h))
  IM = GModule[]
  for m in im
    z = map(x->matrix(map(y->map(hi, y), m(GapObj(x)))), gens(G))
    if ngens(G) == 0
      F = free_module(k, 0)
      zz = typeof(hom(F, F, elem_type(F)[]))[]
    else
      F = free_module(k, nrows(z[1]))
      zz = map(x->hom(F, F, x), z)
    end
    push!(IM, gmodule(F, G, zz))
  end
  return IM
end

function irreducible_modules(G::Oscar.GAPGroup)
  im = GAP.Globals.IrreducibleRepresentations(GapObj(G))
  IM = GModule[]
  K = abelian_closure(QQ)[1]
  for m in im
    z = map(x->matrix(map(y->map(K, y), m(GapObj(x)))), gens(G))
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

@doc raw"""
    trivial_gmodule(G::Group, M::FinGenAbGroup)
    trivial_gmodule(G::Group, M::AbstractAlgebra.FPModule)

Return the `G`-module over the underlying set `M` with trivial `G`-action,
i.e., `g(m) == m` for all $g\in G$ and $m\in M$.
"""
function trivial_gmodule(G::Oscar.GAPGroup, M::Union{FinGenAbGroup, AbstractAlgebra.FPModule})
  I = hom(M, M, gens(M))
  return Oscar.gmodule(M, G, typeof(I)[I for x = gens(G)])
end

function Oscar.gmodule(::Type{AbsSimpleNumField}, M::GModule{<:Oscar.GAPGroup, <:AbstractAlgebra.FPModule{AbsSimpleNumFieldElem}})
  k, mk = Hecke.subfield(base_ring(M), vec(collect(reduce(vcat, map(matrix, M.ac)))))
  if k != base_ring(M)
    F = free_module(k, dim(M))
    return gmodule(group(M), [hom(F, F, map_entries(pseudo_inv(mk), matrix(x))) for x = M.ac])
  end
  return M
end

function Oscar.gmodule(::Type{AbsSimpleNumField}, M::GModule{<:Oscar.GAPGroup, <:AbstractAlgebra.FPModule{QQAbElem{AbsSimpleNumFieldElem}}})
  return gmodule(AbsSimpleNumField, gmodule(CyclotomicField, M))
end

function irreducible_modules(::Type{AbsSimpleNumField}, G::Oscar.GAPGroup; minimal_degree::Bool = false)
  z = irreducible_modules(G)
  Z = GModule[]
  for m in z
    a = gmodule(CyclotomicField, m)
    b = gmodule(AbsSimpleNumField, a)
    push!(Z, b)
  end

  if !minimal_degree
    return Z
  else
    res = GModule[]
    return map(_minimize, Z)
  end
end

function _minimize(V::GModule{<:Oscar.GAPGroup, <:AbstractAlgebra.FPModule{AbsSimpleNumFieldElem}})
  k, m = _character_field(V)
  chi = character(V)
  d = schur_index(chi)
  if d != 1
    @vprint :MinField 1  "non-trivial Schur index $d found\n"
  end
  if d !== nothing && d*degree(k) == degree(base_ring(V))
    return V
  elseif d == 1
    @vprint :MinField 1 "Going from $(degree(base_ring(V))) to $(degree(k))"
    Vmin = gmodule_over(m, V)
    return Vmin
  else
    if d === nothing
      d = 1 # only a lower bound is known
    end
    s = subfields(base_ring(V))
    s = [x for x in s if degree(x[1]) >= d*degree(k)]
    sort!(s, lt = (a,b) -> degree(a[1]) < degree(b[1]))
    for (m, mm) in s
      if m == base_ring(V)
        @vprint :MinField 1 "no smaller field possible\n"
        return V
      end
      @vprint :MinField 1 "trying descent to $m...\n"
      Vmin = gmodule_over(mm, V, do_error = false)
      if Vmin !== nothing
        @vprint :MinField 1 "...success\n"
        return Vmin
      end
    end
  end
end

function irreducible_modules(::typeof(CyclotomicField), G::Oscar.GAPGroup)
  z = irreducible_modules(G)
  return [gmodule(CyclotomicField, m) for m in z]
end

function irreducible_modules(::QQField, G::Oscar.GAPGroup)
  #if cyclo is not minimal, this is not irreducible
  z = irreducible_modules(CyclotomicField, G)
  return [gmodule(QQ, descent_to_minimal_degree_field(m)) for m in z]
end

function irreducible_modules(::ZZRing, G::Oscar.GAPGroup)
  z = irreducible_modules(QQ, G)
  return [gmodule(ZZ, m) for m in z]
end

"""
    gmodule(k::Field, C::GModule)

TODO
"""
function gmodule(a::Type{CyclotomicField}, C::GModule)
  @assert isa(base_ring(C), QQAbField)
  d = dim(C)
  l = 1
  for g = C.ac
    l = lcm(l, lcm(collect(map_entries(x->Hecke.is_cyclotomic_type(parent(x.data))[2], matrix(g)))))
  end
  K = cyclotomic_field(base_ring(C), l)[1]
  F = free_module(K, dim(C))
  if d == 0
    h = hom(F, F, elem_type(F)[])
    return gmodule(F, group(C), typeof(h)[hom(F, F, map_entries(x->K(x.data), matrix(x))) for x = C.ac])
  end
  return gmodule(F, group(C), [hom(F, F, map_entries(x->K(x.data), matrix(x))) for x = C.ac])
end

function gmodule(k::Union{Nemo.fpField, FqField}, C::GModule{<:Oscar.GAPGroup, FinGenAbGroup})
  @assert absolute_degree(k) == 1
  q, mq = quo(C.M, characteristic(k))
  s, ms = snf(q)

  r = ngens(s)
  F = free_module(k, r)
  mp = [FinGenAbGroupHom(ms*pseudo_inv(mq)*x*mq*pseudo_inv(ms)) for x= C.ac]
  return gmodule(F, group(C), [hom(F, F, map_entries(k, x.map)) for x = mp])
end

function gmodule(k::Nemo.fpField, mC::Hecke.MapClassGrp)
  return gmodule(k, gmodule(ray_class_field(mC)))
end

function Base.:^(C::GModule{<:Any, <:AbstractAlgebra.FPModule{AbsSimpleNumFieldElem}}, phi::Map{AbsSimpleNumField, AbsSimpleNumField})
  F = free_module(codomain(phi), dim(C))
  return GModule(group(C), [hom(F, F, map_entries(phi, matrix(x))) for x = C.ac])
end

function Base.:^(C::GModule{<:Any, T}, h::Map{S, S}) where T <: S where S
  return GModule(group(C), [inv(h)*x*h for x = C.ac])
end

function Base.:^(C::GModule{<:Any, <:AbstractAlgebra.FPModule{QQAbElem}}, phi::Map{QQAbField, QQAbField})
  F = free_module(codomain(phi), dim(C))
  return GModule(F, group(C), [hom(F, F, map_entries(phi, matrix(x))) for x = C.ac])
end

function gmodule(::QQField, C::GModule{<:Any, <:AbstractAlgebra.FPModule{AbsSimpleNumFieldElem}})
  F = free_module(QQ, dim(C)*degree(base_ring(C)))
  return GModule(F, group(C), [hom(F, F, hvcat(dim(C), [representation_matrix(x) for x = transpose(matrix(y))]...)) for y = C.ac])
end

gmodule(k::Nemo.fpField, C::GModule{<:Any, <:AbstractAlgebra.FPModule{fpFieldElem}}) = C

function _character(C::GModule{<:Any, <:AbstractAlgebra.FPModule{<:AbstractAlgebra.FieldElem}})
  G = group(C)
  phi = epimorphism_from_free_group(G)
  ac = Oscar.GrpCoh.action(C)
  iac = Oscar.GrpCoh.inv_action(C)

  n = dim(C)
  K = base_ring(C)

  chr = []
  for c = conjugacy_classes(character_table(G))
    r = representative(c)
    if isone(r)
      push!(chr, (c, K(n)))
      continue
    end
    #use T = action(C, r) instead? Trying this, seems to work.
    # p = preimage(phi, r)
    # T = map_word(p, ac; genimgs_inv = iac)
    T = action(C,r)
    push!(chr, (c, trace(matrix(T))))
  end
  return chr
end

"""
Return Z[G] and a function f that, when applied to a G-module M will return
a map representing the action of Z[G] on M:

f(C) yields the extension of g -> action(C, g)

The third value returned is a Map between group elements and the index of the
corresponding module generator.
"""
function natural_gmodule(G::Oscar.GAPGroup, R::Ring)
  M = free_module(R, Int(order(G)))
  ge = collect(G)
  ZG = gmodule(G, [hom(M, M, [M[findfirst(isequal(ge[i]*g), ge)] for i=1:length(ge)]) for g = gens(G)])
  return ZG, C->(x -> sum(x[i]*action(C, ge[i]) for i=1:length(ge))),
    MapFromFunc(G, ZZ, x->ZZ(findfirst(isequal(x), ge)),
             y->ge[Int(y)])
end

function natural_gmodule(::Type{FinGenAbGroup}, G::Oscar.GAPGroup, ::ZZRing)
  M = free_abelian_group(order(Int, G))
  ge = collect(G)
  ZG = gmodule(G, [hom(M, M, [M[findfirst(isequal(ge[i]*g), ge)] for i=1:length(ge)]) for g = gens(G)])
  return ZG, C->(x -> sum(x[i]*action(C, ge[i]) for i=1:length(ge))),
    MapFromFunc(G, ZZ, x->ZZ(findfirst(isequal(x), ge)),
             y->ge[Int(y)])
end

Oscar.character_field(C::GModule{<:Any, <:AbstractAlgebra.FPModule{QQFieldElem}}) = QQ

function _character_field(C::GModule{<:Any, <:AbstractAlgebra.FPModule{AbsSimpleNumFieldElem}})
  val = _character(C)
  k, mkK = Hecke.subfield(base_ring(C), [x[2] for x = val])
  return k, mkK
end

function Oscar.character_field(C::GModule{<:Any, <:AbstractAlgebra.FPModule{AbsSimpleNumFieldElem}})
  return _character_field(C)[1]
end

function Oscar.character(C::GModule{<:Any, <:AbstractAlgebra.FPModule{QQAbElem{AbsSimpleNumFieldElem}}})
  return Oscar.class_function(group(C), [x[2] for x = _character(C)])
end

function Oscar.character(C::GModule{<:Any, <:AbstractAlgebra.FPModule{AbsSimpleNumFieldElem}})
  chr = _character(C)
  k, mkK = Hecke.subfield(base_ring(C), [x[2] for x = chr])
  A = maximal_abelian_subfield(ClassField, k)
  c = Hecke.norm(conductor(A)[1])
  QQAb = abelian_closure(QQ)[1]
  K = cyclotomic_field(QQAb, Int(c))[1]
  fl, em = is_subfield(k, K)
  return Oscar.class_function(group(C), [QQAb(em(preimage(mkK, x[2]))) for x = chr])
end

function Oscar.character(C::GModule{<:Any, <:AbstractAlgebra.FPModule{QQFieldElem}})
  QQAb = abelian_closure(QQ)[1]
  return Oscar.class_function(group(C), [QQAb(x[2]) for x = _character(C)])
end

function Oscar.natural_character(C::GModule{<:Any, <:AbstractAlgebra.FPModule{<:FinFieldElem}})
  G = C.G
  tbl = character_table(G)
  k = base_ring(C.M)
  p = characteristic(k)
  modtbl = mod(tbl, p)
  ccl = conjugacy_classes(modtbl)  # p-regular classes
  h = Oscar.iso_oscar_gap(k)

  vals = [GAP.Globals.BrauerCharacterValue(GAP.Obj(map(h, matrix(action(C, representative(x)))))) for x in ccl]

  return Oscar.class_function(modtbl, GAPWrap.ClassFunction(Oscar.GAPTable(modtbl), Oscar.GapObj(vals)))
end

function Oscar.sub(C::GModule{<:Any, <:AbstractAlgebra.FPModule{T}}, m::MatElem{T}) where {T <: FinFieldElem}

  k = base_ring(C)
  h = Oscar.iso_oscar_gap(k)
  s = GAP.Globals.ShallowCopy(GAP.Obj(map(h, m)))
  g = Gap(C)
  x = GAP.Globals.MTX.SubGModule(g, s)
  b = matrix([preimage(h, x[i, j]) for i in 1:GAPWrap.NrRows(x), j in 1:GAPWrap.NrCols(x)])

  y = GAP.Globals.MTX.InducedActionSubmoduleNB(g, x)
  F = free_module(k, nrows(b))
  return gmodule(F, Group(C), [hom(F, F, matrix([preimage(h, x[i, j]) for i in 1:GAPWrap.NrRows(x), j in 1:GAPWrap.NrCols(x)])) for x = y.generators]), hom(F, C.M, b)
  return b
end

# Compute the restriction of the `M.G`-action from `M.M`
# to the submodule given by the embedding `f`.
function Oscar.sub(M::GModule{<:Any, <:AbstractAlgebra.FPModule{T}}, f::AbstractAlgebra.Generic.ModuleHomomorphism{T}) where T
  @assert codomain(f) == M.M
  S = domain(f)
  Sac = [hom(S, S, [preimage(f, h(f(x))) for x in gens(S)]) for h in M.ac]
  return gmodule(S, M.G, Sac)
end

function gmodule(k::Nemo.FinField, C::GModule{<:Any, <:AbstractAlgebra.FPModule{<:FinFieldElem}})
  @assert absolute_degree(k) == 1
  F = free_module(k, dim(C)*absolute_degree(base_ring(C)))
  return GModule(F, group(C), [hom(F, F, hvcat(dim(C), [absolute_representation_matrix(x) for x = transpose(matrix(y))]...)) for y = C.ac])
end

function Hecke.frobenius(K::FinField, i::Int=1)
  MapFromFunc(K, K, x->Hecke.frobenius(x, i), y -> Hecke.frobenius(x, degree(K)-i))
end

function Hecke.absolute_frobenius(K::FinField, i::Int=1)
  MapFromFunc(K, K, x->Hecke.absolute_frobenius(x, i), y -> Hecke.absolute_frobenius(x, absolute_degree(K)-i))
end

@doc raw"""
    gmodule_minimal_field(C::GModule)
    gmodule_minimal_field(k::Field, C::GModule)

Return TODO
"""
function gmodule_minimal_field(C::GModule{<:Any, <:AbstractAlgebra.FPModule{fpFieldElem}})
  return C
end

function gmodule_minimal_field(C::GModule{<:Any, <:AbstractAlgebra.FPModule{<:FinFieldElem}})
  #always over char field
  K =  base_ring(C)
  d = 0
  while d < absolute_degree(K)-1
    d += 1
    absolute_degree(K) % d == 0 || continue
    k = GF(characteristic(K), d)
    D = gmodule_over(k, C, do_error = false)
    D === nothing || return D
  end
  return C
end

function gmodule_minimal_field(C::GModule{<:Any, <:AbstractAlgebra.FPModule{AbsSimpleNumFieldElem}})
  return _minimize(C)
end

"""
    gmodule_over(k::Field, C::GModule)
"""
function gmodule_over(k::FinField, C::GModule{<:Any, <:AbstractAlgebra.FPModule{<:FinFieldElem}}; do_error::Bool = false)
  #mathematically, k needs to contain the character field
  #only works for irreducible modules
  #requires rel cyclic Galois group, not really finite field...
  #
  K = base_ring(C)
  @assert absolute_degree(K) != absolute_degree(k)
  #method: let s = sigma be a generator for Gal(K/k), and rho the representation
  #attached to C, then if there is A s.th.
  #    A^-1 rho(g) A in GL(k)
  # then
  #    A^-s rho(g)^s A^s = A^-1 rho(g) A
  # so
  #    rho(g)^s = A^sA^-1 rho(g) A A^-s
  # Let thus B s.th.
  #    rho(g)^s = B^-1 rho(g) B
  # so
  #    rho(g)^(s^2) = (rho(g)^s)^s =
  #                 = (B^-1 rho(g) B)^s
  #                 = B^-s rho(g)^s B^s
  #                 = B^-s B^-1 rho(g) B B^s
  # inductively:
  #    rho(g)^(s^i) = B^-(s^i-1) B^-(s^(i-1)) ... B^-1 rho(g) B B^s ...
  # From s^n = 1, we obtain N(B) = prod_i=0^n-1 B^(s^i) = lambda I
  # (since rho is irreducible and thus the matrix unique up to scalar)
  # IF B is, as above A^(1-s), then N(B) = 1, so there should be
  # alpha s.th. N(alpha) = lambda
  # thus N(alpha^-1 B) = I
  # Hilbert 90: alpha^-1 B = S^(1-s) and we can use this S to conjugate down

  # ALGO
  s = absolute_frobenius(K, absolute_degree(k))
  mkK = embed(k, K)
  os = divexact(absolute_degree(K), absolute_degree(k))
  hB = hom_base(C, gmodule(C.M, Group(C),
                      [hom(C.M, C.M, map_entries(s, matrix(x))) for x = C.ac]))
  if length(hB) != 1
    if do_error
      length(hB) == 0 && error("Module cannot be written over $k")
      length(hB) > 1  && error("Module not irreducible, hom too large")
    end
    return nothing
  end
  B = hB[1]
  D = norm(B, s, os)
  lambda = D[1,1]
  @hassert :MinField 2 D == lambda*identity_matrix(K, dim(C))
  alpha = norm_equation(K, preimage(mkK, lambda))
  B *= inv(alpha)
  @hassert :MinField 2 isone(norm(B, s, os))
  D = hilbert90_cyclic(B, s, os)
  Di = inv(D)
  F = free_module(k, dim(C))
  return gmodule(F, Group(C), [hom(F, F, map_entries(x -> preimage(mkK, x), Di*matrix(x)*D)) for x = C.ac])
  # return C^-1 x C for x = action_gens(C), coerced into k
end

#...now the same for number fields - and non-cyclic fields.
function gmodule_over(em::Map{AbsSimpleNumField, AbsSimpleNumField}, C::GModule{<:Any, <:AbstractAlgebra.FPModule{AbsSimpleNumFieldElem}}; do_error::Bool = true)
  K = base_ring(C)
  k = domain(em)
  @assert codomain(em) == K
  gk = em(gen(k))

  A, mA = automorphism_group(PermGroup, K)
  s, ms = sub(A, [a for a = A if mA(a)(gk) == gk])
  ac =  _two_cocycle(ms*mA, C, do_error = do_error)
  if ac === nothing
    if do_error
      error("cannot do over this field")
    else
      return nothing
    end
  end
  F = free_module(k, dim(C))
  return gmodule(F, group(C), [hom(F, F, map_entries(t->preimage(em, t), x)) for x = ac])
end

function gmodule_over(::QQField, C::GModule{<:Any, <:AbstractAlgebra.FPModule{QQAbElem{AbsSimpleNumFieldElem}}}; do_error::Bool = true)
  return gmodule_over(QQ, gmodule(CyclotomicField, C); do_error)
end

function gmodule_over(::QQField, C::GModule{<:Any, <:AbstractAlgebra.FPModule{AbsSimpleNumFieldElem}}; do_error::Bool = true)
  K = base_ring(C)
  A, mA = automorphism_group(PermGroup, K)
  ac = _two_cocycle(mA, C, do_error = do_error)
  if ac === nothing
    if do_error
      error("cannot do over this field")
    else
      return nothing
    end
  end
  F = free_module(QQ, dim(C))
  return gmodule(F, group(C), [hom(F, F, map_entries(QQ, x)) for x = ac])
end

function Oscar.hom(M::MultGrp, N::MultGrp, h::Map)
  return MapFromFunc(M, N, x->N(h(x.data)), y->M(preimage(h, y.data)))
end

function Oscar.content_ideal(M::MatElem{AbsSimpleNumFieldElem})
  zk = maximal_order(base_ring(M))
  C = fractional_ideal(zk, 1*zk)
  if nrows(M)*ncols(M) == 0
    return C
  end

  for i=1:nrows(M)
    for j=1:ncols(M)
      if !iszero(M[i,j])
        C += M[i,j]*zk
      end
    end
  end
  return C
end

"""
Compute the factor set or 2-cochain defined by `C` as a Galois
module of the automorphism group over the character field.
If `mA` is given, it needs to map the automorphism group over the
character field into the the automorphisms of the base ring.
"""
function factor_set(C::GModule{<:Any, <:AbstractAlgebra.FPModule{AbsSimpleNumFieldElem}}, mA::Union{Map, Nothing} = nothing)
  K = base_ring(C)
  if mA === nothing
    k, mkK = _character_field(C)
    A, mA = automorphism_group(PermGroup, K)
    if degree(k) > 1
      gk = mkK(gen(k))
      s, ms = sub(A, [g for g = A if g(gk) == gk])
      mA = ms*mA
    end
  end

  c = _two_cocycle(mA, C, do_error = true, two_cycle = true)
  return c
end

function _two_cocycle(mA::Map, C::GModule{<:Any, <:AbstractAlgebra.FPModule{AbsSimpleNumFieldElem}}; do_error::Bool = true, two_cycle::Bool = false)
  G = domain(mA)
  K = base_ring(C)

  homs = []
  @vprint :MinField 1 "Gathering Galois images of the generators...\n"
  for g = gens(G)
    @vprint :MinField 2 "gen: $g\n"
    @vtime :MinField 2 hb = hom_base(C^mA(g), C)
    #C^g * hb == hb * C
    if length(hb) == 0
      do_error && return nothing
      error("field too small")
    end
    if length(hb) > 1
      do_error && return nothing
      error("rep. not abs. irr.")
    end
    #as the matrices are only unique up to scalars, try to
    #"reduce" them via the content(ideal)...
    @vprint :MinField 2 "trying to (size) reduce matrix...\n"
    @vprint :MinField 3 "from\n$(hb[1])\n"
    c = content_ideal(hb[1])
    d = Hecke.short_elem(inv(c))
    @vprint :MinField 3 "via $d to\n"
    push!(homs, d*hb[1])
    @vprint :MinField 3 "$(homs[end])\n"
  end
  I = identity_matrix(K, dim(C))

  @vprint :MinField 1 "computing un-normalized 1-chain (of matrices)\n"
  # pairs: (g, X_g) with operation (g, X_g)(h, X_h) = (gh, X_g^h * X_h)
  @vtime :MinField 2
    c = closure([(gen(G, i), homs[i]) for i=1:ngens(G)],
              (a, b) -> (a[1]*b[1], map_entries(mA(b[1]), a[2])*b[2]),
              (one(G), I),
              eq = (a,b) -> a[1] == b[1])
  X = Dict(x[1] => x[2] for x = c)
  X[one(G)] = I

  #now we need a 2-cycle:
  #X[g] X[h] = sigma(g, h) X[gh] should hold...

  @vprint :MinField 1 "now the 2-cocycle (scalars)\n"
  MK = MultGrp(K)
  sigma = Dict{Tuple{PermGroupElem, PermGroupElem}, AbsSimpleNumFieldElem}()
  for g = G
    for h = G
      if isone(g)
        sigma[(g, h)] = (one(K))
      elseif isone(h)
        sigma[(g, h)] = (one(K))
      else
        lf = findfirst(x->!iszero(x), X[g*h])
        sigma[(g, h)] = (X[g*h][lf]//(map_entries(mA(h), X[g])*X[h])[lf])
#        sigma[(g, h)] = MK(X[g*h][lf]//(X[h]*map_entries(mA(h), X[g]))[lf])
      end
    end
  end
#  istwo_cocycle(sigma, mA)

  @vprint :MinField 1 "test for co-boundary\n"
  D = gmodule(G, [hom(MK, MK, mA(x)) for x = gens(G)])
  Sigma = Oscar.GrpCoh.CoChain{2,PermGroupElem, MultGrpElem{AbsSimpleNumFieldElem}}(D, Dict(k=>MK(v) for (k,v) = sigma))
  if two_cycle
    return Sigma
  end
  @vtime :MinField 2 fl, cb = Oscar.GrpCoh.is_coboundary(Sigma)

  if !fl
    do_error || return nothing
    error("field too small")
  end

  cc = Dict(k => cb(k).data for k = G)
  for g = G
    for h = G
      @assert mA(h)(cc[g])*cc[h] == sigma[(g, h)]*cc[g*h]
    end
  end

  if !fl
    do_error || return nothing
    error("field too small")
  end
  for g = G
    X[g] *= (cb(g).data)
  end
  isone_cochain(X, mA)

  #now X should be in H^1(G, Gl(n, K)) which is trivial
  #hence a Hilbert-90 should find A s.th. A^(1-g) = X[g] for all g

  @vprint :MinField 1 "calling Hilbert-90 on matrices\n"
  @vtime :MinField 2 A, Ai = hilbert90_generic(X, mA)
  c = content_ideal(A)
  d = Hecke.short_elem(inv(c))
  A *= d
  Ai *= inv(d)

  @vprint :MinField 1 "conjugating the generators\n"
  @vtime :MinField 2 r = [A*matrix(x)*Ai for x = C.ac]
  return r
end

function isone_cochain(X::Dict{<:GAPGroupElem, <:MatElem{AbsSimpleNumFieldElem}}, mA)
  G = domain(mA)
  for g = G
    for h = G
      @assert X[g*h] == map_entries(mA(h), X[g])*X[h]
    end
  end
end

function isone_cochain(X::Dict{<:GAPGroupElem, AbsSimpleNumFieldElem}, mA)
  G = domain(mA)
  for g = G
    for h = G
      @assert X[g*h] == mA(h)(X[g])*X[h]
    end
  end
end

function istwo_cocycle(X::Dict, mA, op = *)
  G = domain(mA)
  for g = G
    for h = G
      for k = G
        #= if (g*h)(x) = g(h(x)), then the cocycle should be
             X[(g*h, k)] X[(g, h)] == mA(g)(X[(h, k)]) X[(g, hk)]
           if (g*h)(x) = h(g(x)) then we should get
             X[(g, hk)] X[(h, k)]  == mA(k)(X[(g, h)]) X[(gh, k)]

             (Debeerst, PhD, (1.1) & (1.2))

             However, if we mix the conventions, all bets are off...
        =#
        a = op(X[(g, h*k)], X[(h, k)]) - op(mA(k)(X[(g, h)]), X[(g*h, k)])
        @show a, iszero(a) || valuation(a)
      end
    end
  end
end

"""
  Hilbert-90: H^1(G, Gl(n, K)) = 1
for G = aut(K) and any number field K.

`X` is a 1-chain, X_g = X[g]. This will find a matrix S s.th.
  S^(1-g) = X_g
for all g.

`G` is both the keys of `X` and the domain of `mA`.

Call at your peril. Used in writing a gmodule over a different
number field.
"""
function hilbert90_generic(X::Dict, mA)
  G = domain(mA)
  K = domain(mA(one(G))) #can map parent do this better?
  n = nrows(first(values(X)))
  cnt = 0
  while true
    local Y
    while true #TODO: choose Y more sparse
      #Glasby shows that this approach, over a finite field,
      #has a high success probability.
      Y = matrix(K, n, n, [rand(K, -5:5) for i=1:n*n])
      fl = is_invertible(Y)
      fl && break
      cnt += 1
      if cnt > 10 error("s.th. weird") end
    end
    S = sum(map_entries(mA(g), Y)*v for (g,v) = X)
    fl, Si = is_invertible_with_inverse(S)
    fl && return S, Si
  end
end

"""
Norm of A wrt. s. s acts on the entries of A, this computes
  A * A^s * A^s^2 ... A^s^(os-1)
os is meant to be the order of s
"""
function norm(A::MatElem, s, os::Int)
  B = A
  C = map_entries(s, A)
  for i=1:os-1
    B *= C
    C = map_entries(s, C)
  end
  return B
end

"""
A needs to have norm 1 wrt. to s, so
  A A^s, .. A^(s^n) = I
for ord(s) = os = n+1. Then this will find B s.th.
  A = B^(1-s)
"""
function hilbert90_cyclic(A::MatElem{<:FieldElem}, s, os::Int)
  #apart from rand, this would also work over a number field
  k= base_ring(A)
  if isa(k, FinField)
    rnd = ()->rand(parent(A))
  elseif isa(k, QQField) || isa(k, AbsSimpleNumField)
    rnd = ()->rand(parent(A), -10:10)
  end
  cnt = 1
  while true
    B = rnd()
    Bs = map_entries(s, B)
    As = A
    for i=1:os-1
      B += As*Bs
      As = A*map_entries(s, As)
      Bs = map_entries(s, Bs)
    end
    if !iszero(det(B))
      @hassert :MinField 2 A == B*inv(map_entries(s, B))
      return B
    else
      if cnt > 10 && get_assertion_level(:MinField) > 1
        error("")
      end
      cnt += 1
    end
  end
end

function gmodule(k::Nemo.fpField, C::GModule{<:Any, <:AbstractAlgebra.FPModule{QQFieldElem}})
  F = free_module(k, dim(C))
  return GModule(group(C), [hom(F, F, map_entries(k, matrix(x))) for x=C.ac])
end

function gmodule(mk::Map{AbsSimpleNumField, <:FinField}, C::GModule{<:Any, <:AbstractAlgebra.FPModule{AbsSimpleNumFieldElem}})
  k = codomain(mk)
  @assert domain(mk) == base_ring(C)
  F = free_module(k, dim(C))
  return GModule(group(C), [hom(F, F, map_entries(mk, matrix(x))) for x=C.ac])
end

function Hecke.modular_proj(C::GModule{T, <:AbstractAlgebra.FPModule{AbsSimpleNumFieldElem}}, me::Hecke.modular_env) where T
  R = []
  z = map(x->(Hecke.modular_proj(x.matrix, me)), C.ac)
  for i=1:length(z[1])
    F = free_module(base_ring(z[1][i]), dim(C))
    @assert all(j->base_ring(z[j][i]) == base_ring(z[1][i]), 1:length(z))
    push!(R, GModule(group(C), [hom(F, F, t[i]) for t = z]))
    @assert all(i->base_ring(matrix(R[end].ac[i])) == base_ring(R[end]), 1:length(R[end].ac))
  end
  return R
end

@attr Oscar.GapObj function Gap(C::GModule{<:Any, <:AbstractAlgebra.FPModule{<:FinFieldElem}})
  h = Oscar.iso_oscar_gap(base_ring(C))
  mats = [GAP.Obj(map(h, Matrix(matrix(x)))) for x in C.ac]
  return GAP.Globals.GModuleByMats(GAP.Obj(mats), codomain(h))
end

function Oscar.is_irreducible(C::GModule{<:Any, <:AbstractAlgebra.FPModule{<:FinFieldElem}})
  G = Gap(C)
  return GAP.Globals.MTX.IsIrreducible(G)
end

function Oscar.is_absolutely_irreducible(C::GModule{<:Any, <:AbstractAlgebra.FPModule{<:FinFieldElem}})
  G = Gap(C)
  return GAP.Globals.MTX.IsAbsolutelyIrreducible(G)
end

function Oscar.splitting_field(C::GModule{<:Any, <:AbstractAlgebra.FPModule{<:FinFieldElem}})
  fl = is_irreducible(C)
  fl || error("module must be irreducible")
  is_absolutely_irreducible(C) #for GAP to actually compute the field
  G = Gap(C)
  d = GAP.Globals.MTX.DegreeSplittingField(G)
  return GF(characteristic(base_ring(C)), d)
end


function is_decomposable(C::GModule{<:Any, <:AbstractAlgebra.FPModule{<:FinFieldElem}})
  G = Gap(C)
  return !GAP.Globals.MTX.IsIndecomposable(G)
end

"""
    composition_factors_with_multiplicity(C::GModule{<:Any, <:AbstractAlgebra.FPModule{<:FinFieldElem}})

Return the composition factors of `C` with their frequency.
"""
function Oscar.composition_factors_with_multiplicity(C::GModule{<:Any, <:AbstractAlgebra.FPModule{<:FinFieldElem}})
  G = Gap(C)
  z = GAP.Globals.MTX.CollectedFactors(G)
  g = Group(C)
  k = base_ring(C)

  CF = []
  for c = z
    m = GAP.Globals.MTX.Generators(c[1])
    mm = [matrix(k, x) for x = m]
    F = free_module(k, nrows(mm[1]))
    push!(CF, (gmodule(F, Group(C), [hom(F, F, x) for x = mm]), c[2]))
  end
  return CF
end

"""
    indecomposition(C::GModule{<:Any, <:AbstractAlgebra.FPModule{<:FinFieldElem}})

Return a decomposition of the module `C` into indecomposable summands as a list
of pairs:
 - a direct indecomposable summand
 - a homomorphism (embedding) of the underlying free modules
"""
function indecomposition(C::GModule{<:Any, <:AbstractAlgebra.FPModule{<:FinFieldElem}})
  G = Gap(C)
  z = GAP.Globals.MTX.Indecomposition(G)
  k = base_ring(C)
  CF = []
  for c = z
    m = GAP.Globals.MTX.Generators(c[2])
    mm = [matrix(k, x) for x = m]
    F = free_module(k, nrows(mm[1]))
    push!(CF, (gmodule(F, Group(C), [hom(F, F, x) for x = mm]), hom(F, C.M, matrix(k, c[1]))))
  end
  return CF
end

"""
The group of Z[G]-homomorphisms as a k-module, not a k[G] one. (The G operation
is trivial)
"""
function Oscar.hom(C::T, D::T) where T <: GModule{<:Any, <:AbstractAlgebra.FPModule{<:FieldElem}}
  b = hom_base(C, D)
  H, mH = hom(C.M, D.M)
  s, ms = sub(H, [H(vec(collect(x))) for x = b])
  return GModule(group(C), [hom(s, s, [preimage(ms, H(vec(collect(inv(matrix(C.ac[i]))*g*matrix(D.ac[i]))))) for g = b]) for i=1:ngens(group(C))]), ms * mH
end

#a bad implementation of hom: as the fix module of ghom
#but we don't have he meataxe in general.
function Hecke.hom(C::GModule{T, FinGenAbGroup}, D::GModule{T, FinGenAbGroup}) where {T}

 H, mH = Oscar.GModuleFromGap.ghom(C, D)
 q, mq = H_zero(H)
 mmH = hom(q, H.M, [mq(x)() for x = gens(q)])
 return q, mmH*mH
 return gmodule(C.G, [hom(q, q, gens(q)) for x = gens(C.G)]), mmH*mH
end

"""
The G-module of all Z-module homomorphisms
"""
function ghom(C::GModule, D::GModule)
  @assert C.G === D.G
  H, mH = hom(C.M, D.M)
  return gmodule(H, C.G, [hom(H, H, [preimage(mH, action(C, inv(g))*mH(h)*action(D, g)) for h = gens(H)]) for g = gens(C.G)]), mH
end

function is_G_hom(C::GModule, D::GModule, H::Map)
  return all([H(action(C, g, h)) == action(D, g, H(h)) for g = gens(C.G) for h = gens(C.M)])
end

#=
G = cyclic_group(PermGroup, 3)
A = abelian_group([3, 3])
C = gmodule(G, [hom(A, A, [A[1], A[1]+A[2]])])
is_consistent(C)

zg, ac = Oscar.GModuleFromGap.natural_gmodule(G, ZZ)
zg = gmodule(FinGenAbGroup, zg)
H, mH = Oscar.GModuleFromGap.ghom(zg, C)
inj = hom(C.M, H.M, [preimage(mH, hom(zg.M, C.M, [ac(C)(g)(c) for g = gens(zg.M)])) for c = gens(C.M)])

q, mq = quo(H, image(inj)[2])
=#

function hom_base(C::GModule{S, <:AbstractAlgebra.FPModule{T}}, D::GModule{S, <:AbstractAlgebra.FPModule{T}}) where {S <: Oscar.GAPGroup, T <: FinFieldElem}
  @assert base_ring(C) == base_ring(D)
  h = Oscar.iso_oscar_gap(base_ring(C))
  hb = GAP.Globals.MTX.BasisModuleHomomorphisms(Gap(C), Gap(D))
  n = length(hb)
  b = dense_matrix_type(base_ring(C))[matrix([preimage(h, x[i, j]) for i in 1:GAPWrap.NrRows(x), j in 1:GAPWrap.NrCols(x)]) for x in hb]

#  @show [matrix(C.ac[i])*b[1] == b[1]*matrix(D.ac[i]) for i=1:length(C.ac)]
  return b
end

function hom_base(C::T, D::T) where T <: GModule{<:Any, <:Generic.FreeModule{<:AlgClosureElem{<:FinField}}}

  C1 = gmodule(FinField, C)
  D1 = gmodule(FinField, D)
  Cf = degree(base_ring(C1))
  Df = degree(base_ring(D1))
  l = lcm(Cf, Df)
  K = ext_of_degree(base_ring(C), l)
  if l != Cf
    C1 = gmodule(K, C1)
  end
  if l != Df
    D1 = gmodule(K, D1)
  end
  h = Oscar.GModuleFromGap.hom_base(C1, D1)
  if length(h) == 0
    return h
  end
  return map(x->map_entries(base_ring(C), x), h)
end

#TODO: in ctx of MeatAxe & Gap: we're mostly having a rref,
#      but for a different ordering of entries
function _rref!(V::Vector{<:MatElem{<:FieldElem}})
  #@show :in, V
  @assert all(x->size(x) == size(V[1]), V)
  n = nrows(V[1])
  @assert ncols(V[1]) == n

  o = 1
  for i = CartesianIndices((1:n, 1:n))
    j = findall(x->!iszero(x[i]), V)
    isempty(j) && continue
    if j[1] != o
      V[o], V[j[1]] = V[j[1]], V[o]
      j[1] = o
    end
    if !isone(V[o][i])
      V[o] *= inv(V[o][i])
    end
    for k=o+1:length(V)
      iszero(V[k][i]) && continue
      V[k] -= V[k][i] * V[o]
    end
    o += 1
    if o>length(V)
      return
    end
  end
end

"""
  C*T[i] = T[i]*D
on return.

Currently assumes no bad primes.
"""
function hom_base(C::GModule{<:Any, <:AbstractAlgebra.FPModule{AbsSimpleNumFieldElem}}, D::GModule{<:Any, <:AbstractAlgebra.FPModule{AbsSimpleNumFieldElem}})
  @assert base_ring(C) == base_ring(D)

  p = Hecke.p_start
  p = 2^10
#  p = 127
  m_in = map(matrix, C.ac)
  m_out = map(matrix, D.ac)
  local T
  pp = ZZRingElem(1)
  k = base_ring(C)
  @assert base_ring(m_in[1]) == k
  @assert base_ring(m_in[1]) == k
  while true
    p = next_prime(p)
    me = modular_init(k, p, deg_limit = 1)
    isempty(me) && continue
    z1 = Hecke.modular_proj(C, me)
    if C === D
      z2 = z1
    else
      z2 = Hecke.modular_proj(D, me)
    end
    t = []
    for i=1:length(z1)
      mp = hom_base(z1[i], z2[i])
      if isempty(mp)
        return dense_matrix_type(base_ring(C))[]
      end
      _rref!(mp)
      push!(t, mp)
    end
    #should actually compute an rref of the hom base to make sure
    if any(x->length(x) != length(t[1]), t)
      @show :BP
      #bad prime...
      continue
    end
    if length(t[1]) == 0
      return []
    end

    tt = [Hecke.modular_lift([t[i][j] for i=1:length(z1)], me) for j=1:length(t[1])]
    @assert base_ring(tt[1]) == k
    if isone(pp)
      pp = ZZRingElem(p)
      T = tt
    else
      T = [induce_crt(tt[i], T[i], ZZRingElem(p), pp) for i=1:length(T)]
      @assert base_ring(T[1]) == k
      pp *= p
      S = []
      for t = T
        fl, s = induce_rational_reconstruction(t, pp, ErrorTolerant = true)
        fl || break
        push!(S, s)
      end
      length(S) == 0 && continue
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

function hom_base(C::_T, D::_T) where _T <: GModule{<:Any, <:AbstractAlgebra.FPModule{QQFieldElem}}
  @assert base_ring(C) == base_ring(D)

  p = Hecke.p_start
  p = 2^10
  p = 127
  m_in = map(matrix, C.ac)
  m_out = map(matrix, D.ac)
  local T
  pp = ZZRingElem(1)
  k = base_ring(C)
  @assert base_ring(m_in[1]) == k
  @assert base_ring(m_in[1]) == k
  @assert k == QQ
  #a heuristic when to try to call reconstruct...
  bt = maximum(maximum(nbits, matrix(x)) for x = vcat(C.ac, D.ac)) * dim(C)
  reco = 10
  while true
    p = next_prime(p)
    z1 = gmodule(Native.GF(p), C)
    if C === D
      z2 = z1
    else
      z2 = gmodule(base_ring(z1), D)
    end

    t = hom_base(z1, z2)

    isempty(t) && return QQMatrix[]
    _rref!(t)
    tt = [lift(s)  for s=t]
    @assert base_ring(tt[1]) == ZZ
    if isone(pp)
      pp = ZZRingElem(p)
      T = tt
    else
      T = [induce_crt(tt[i], ZZRingElem(p), T[i], pp)[1] for i=1:length(T)]
      @assert base_ring(T[1]) == ZZ
      pp *= p
      S = []
      if nbits(pp) > min(reco, bt)
        if nbits(pp) > reco
          reco *= 2
        end
        for t = T
          fl, s = induce_rational_reconstruction(t, pp, ErrorTolerant = true)
          fl || break
          push!(S, s)
        end
      end
      if nbits(pp) > 1000 && get_assertion_level(:MinField) > 1
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

function gmodule(K::AbsSimpleNumField, M::GModule{<:Any, <:AbstractAlgebra.FPModule{AbsSimpleNumFieldElem}})
  F = free_module(K, dim(M))
  return gmodule(F, group(M), [hom(F, F, map_entries(K, matrix(x))) for x = M.ac])
end

function hom_base(C::_T, D::_T) where _T <: GModule{<:Any, <:AbstractAlgebra.FPModule{<:QQAbElem}}
  C1 = gmodule(CyclotomicField, C)
  D1 = gmodule(CyclotomicField, D)
  fl, Cf = Hecke.is_cyclotomic_type(base_ring(C1))
  @assert fl
  fl, Df = Hecke.is_cyclotomic_type(base_ring(D1))
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

function gmodule(::QQField, C::GModule{<:Any, <:AbstractAlgebra.FPModule{ZZRingElem}})
  F = free_module(QQ, dim(C))
  return GModule(group(C), [hom(F, F, map_entries(QQ, matrix(x))) for x = C.ac])
end

function hom_base(C::_T, D::_T) where _T <: GModule{<:Any, <:AbstractAlgebra.FPModule{ZZRingElem}}

  h = hom_base(gmodule(QQ, C), gmodule(QQ, D))
  H = reduce(vcat, [integral_split(matrix(QQ, 1, dim(C)^2, vec(collect(x))), ZZ)[1] for x = h])
  H = Hecke.saturate(H)
  return [matrix(ZZ, dim(C), dim(C), vec(collect(H[i, :]))) for i=1:nrows(H)]
end

function gmodule(::ZZRing, C::GModule{<:Any, <:AbstractAlgebra.FPModule{QQFieldElem}})
  ma = map(matrix, C.ac)
  M = identity_matrix(QQ, dim(C))
  if dim(C) == 0
    F = free_module(ZZ, 0)
    return gmodule(F, group(C), [hom(F, F, matrix(ZZ, 0, 0, ZZRingElem[])) for g = gens(group(C))])
  end
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

function Base.transpose(C::GModule{<:Any, <:AbstractAlgebra.FPModule})
  return gmodule(group(C), [transpose(x) for x = action_matrices(C)])
end

function Oscar.dual(C::GModule{<:Any, <:AbstractAlgebra.FPModule})
  D = gmodule(group(C), [inv(transpose(x)) for x = action_matrices(C)])
  return D
end

#if C is abs. irr <=> hom is 1 dim => this always works
#otherwise, this gives a basis of the symmetric matrices, stabilized
#by the group - but possibly not pos. definite
#
#to always get pos. def. one needs a different algorithm
# - sum over group
# - approximate reynolds as in invar thy (for number fields and Q/ Z)
function invariant_forms(C::GModule{<:Any, <:AbstractAlgebra.FPModule})
  D = Oscar.dual(C)
  h = hom_base(C, D)
  k = kernel(transpose(reduce(vcat, [matrix(base_ring(C), 1, dim(C)^2, _vec(x-transpose(x))) for x = h])))
  return [sum(h[i]*k[i, j] for i=1:length(h)) for j=1:ncols(k)]
end

function Oscar.is_isomorphic(A::GModule{T, <:AbstractAlgebra.FPModule{ZZRingElem}}, B::GModule{T, <:AbstractAlgebra.FPModule{ZZRingElem}}) where T

  h = hom_base(gmodule(QQ, A), gmodule(QQ, B))

  if length(h) == 0
    return false
  end

  if length(h) > 1
    error("modules are not abs. irred.")
  end
  S = h[1]
  x = findfirst(!iszero, S[1, :])
  S *= inv(S[1, x])
  return denominator(S) == 1 && abs(det(S)) == 1
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

function Oscar.gmodule(::Type{FinGenAbGroup}, C::GModule{T, <:AbstractAlgebra.FPModule{FqFieldElem}}) where {T <: Oscar.GAPGroup}
  k = base_ring(C.M)
  A = abelian_group([characteristic(k) for i=1:rank(C.M)*absolute_degree(k)])
  return GModule(group(C), [hom(A, A, map_entries(x->lift(ZZ, x), hvcat(dim(C), [absolute_representation_matrix(x) for x = transpose(matrix(y))]...))) for y = C.ac])
end

function Oscar.gmodule(::Type{FinGenAbGroup}, C::GModule{T, <:AbstractAlgebra.FPModule{<:Union{FpFieldElem, fpFieldElem}}}) where {T <: Oscar.GAPGroup}
  A = abelian_group([characteristic(base_ring(C)) for i=1:rank(C.M)])
  return Oscar.gmodule(A, Group(C), [hom(A, A, map_entries(lift, matrix(x))) for x = C.ac])
end

function Oscar.gmodule(::Type{FinGenAbGroup}, C::GModule{T, <:AbstractAlgebra.FPModule{ZZRingElem}}) where {T <: Oscar.GAPGroup}
  A, mA = abelian_group(C.M)
  return Oscar.gmodule(A, Group(C), [hom(A, A, matrix(x)) for x = C.ac])
end

function Oscar.gmodule(chi::Oscar.GAPGroupClassFunction)
  f = GAP.Globals.IrreducibleAffordingRepresentation(chi.values)
  K = abelian_closure(QQ)[1]
  g = GAP.Globals.List(GAP.Globals.GeneratorsOfGroup(GapObj(group(chi))), f)
  z = map(x->matrix(map(y->map(K, y), g[x])), 1:GAP.Globals.Size(g))
  F = free_module(K, degree(Int, chi))
  return gmodule(group(chi), [hom(F, F, x) for x = z])
end

function Oscar.gmodule(T::Union{Type{CyclotomicField}, Type{AbsSimpleNumField}}, chi::Oscar.GAPGroupClassFunction)
  return gmodule(T, gmodule(chi))
end

function Oscar.gmodule(::Type{FinGenAbGroup}, C::GModule{T, AbstractAlgebra.FPModule{fqPolyRepFieldElem}}) where {T <: Oscar.GAPGroup}
  k = base_ring(C)
  A, mA = abelian_group(C.M)

  return Oscar.gmodule(A, Group(C), [hom(A, A, [preimage(mA, x(mA(a))) for a = gens(A)]) for x = C.ac])
end

function (f::Map{<:Oscar.GAPGroup, <:Oscar.GAPGroup})(C::GModule)
  @assert codomain(f) == Group(C)
  return Oscar.gmodule(C.M, domain(f), [action(C, f(x)) for x = gens(domain(f))])
end

#TODO: cover all finite fields
#      make the Modules work

function Oscar.simplify(C::GModule{<:Any, <:AbstractAlgebra.FPModule{QQFieldElem}})
  return gmodule(QQ, Oscar.simplify(gmodule(ZZ, C))[1])
end

function action_matrices(C::GModule{<:Any, <:AbstractAlgebra.FPModule})
  return map(matrix, action(C))
end

function Oscar.simplify(C::GModule{<:Any, <:AbstractAlgebra.FPModule{ZZRingElem}})
 f = invariant_forms(C)[1]
 # global last_f = f
 @assert all(i->det(f[1:i, 1:i])>0, 1:nrows(f))
 m = map(matrix, C.ac)
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

export extension_of_scalars
export factor_set
export ghom
export indecomposition
export irreducible_modules
export is_decomposable
export is_G_hom
export restriction_of_scalars
export trivial_gmodule
export gmodule_minimal_field
export gmodule_over

## Fill in some stubs for Hecke

function _to_gap(h, x::Vector)
  return GAP.Globals.GModuleByMats(GAP.Obj([GAP.Obj(map(h, Matrix(y))) for y in x]), codomain(h))
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
    h = Oscar.iso_oscar_gap(F)
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
    h = Oscar.iso_oscar_gap(F)
    @assert base_ring(x[1]) == base_ring(y[1])
    @assert length(x) == length(y)
    hb = GAP.Globals.MTX.BasisModuleHomomorphisms(_to_gap(h, x), _to_gap(h, y))
    hbb = [_gap_matrix_to_julia(h, g) for g in GAP.gap_to_julia(hb)]
    return hbb
  end
end

end #module GModuleFromGap

using .GModuleFromGap

export extension_of_scalars
export factor_set
export ghom
export indecomposition
export irreducible_modules
export is_decomposable
export is_G_hom
export restriction_of_scalars
export trivial_gmodule
export gmodule_minimal_field
export gmodule_over

include("Brueckner.jl")
