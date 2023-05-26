module GModuleFromGap
using Oscar
using Hecke
import Hecke: data


"""
    restriction_of_scalars(M::GModule, phi::Map)

Return the `S`-module obtained by restricting the scalars of `M`
from `R` to `S`, where `phi` is an embedding of `S` into `R`.

If `R` has `S`-rank `d` and `M` has rank `n` then the returned module
has rank `d*n`.

# Examples

(cases that `R` is a finite field or a number field, and `S` is a subfield)
"""
function restriction_of_scalars(M::GModule, phi::Map)
  error("not yet ...")
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
  error("not yet ...")
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
function descent_to_minimal_degree_field(M::GModule)
  error("not yet ...")
end


"""
    invariant_lattice_classes(M::GModule, phi::Map)

Return representatives of the equivalence classes of `G`-invariant
`S`-lattices in the `R`-module `M`,
where `phi` is an embedding of `S` into `R`.

# Examples

(case that `R` is a number field and `S` is the ring of integers in `R`)
"""
function invariant_lattice_classes(M::GModule, phi::Map)
  error("not yet ...")
end


#XXX: clash of names!
#   gmodule(k, C) vs gmodule_ver(k, C)
#   first does a "restriction of scalars" or blow up with the rep mat
#   second tries to conjugate down to k

import Oscar:gmodule, GAPWrap
import Oscar.GrpCoh: MultGrp, MultGrpElem

import AbstractAlgebra: Group, Module
import Base: parent

function __init__()
  add_verbose_scope(:BruecknerSQ)
  set_verbose_level(:BruecknerSQ, 0)

  add_assert_scope(:BruecknerSQ)
  set_assert_level(:BruecknerSQ, 0)

  add_assert_scope(:MinField)
  set_assert_level(:MinField, 0)

  add_verbose_scope(:MinField)
  set_verbose_level(:MinField, 0)
end

function Hecke.number_field(::QQField, chi::Oscar.GAPGroupClassFunction; cached::Bool = false)
  return number_field(QQ, map(x->GAP.gap_to_julia(QQAbElem, x), chi.values), cached = cached)
end

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

function Oscar.gmodule(::Type{AnticNumberField}, M::GModule{<:Oscar.GAPGroup, Generic.FreeModule{nf_elem}})
  k, mk = Hecke.subfield(base_ring(M), vec(collect(vcat(map(mat, M.ac)...))))
  if k != base_ring(M)
    F = free_module(k, dim(M))
    return gmodule(group(M), [hom(F, F, map_entries(pseudo_inv(mk), mat(x))) for x = M.ac])
  end
  return M
end

function Oscar.gmodule(::Type{AnticNumberField}, M::GModule{<:Oscar.GAPGroup, Generic.FreeModule{QQAbElem{nf_elem}}})
  return gmodule(AnticNumberField, gmodule(CyclotomicField, M))
end

function irreducible_modules(::Type{AnticNumberField}, G::Oscar.GAPGroup; minimal_degree::Bool = false)
  z = irreducible_modules(G)
  Z = GModule[]
  for m in z
    a = gmodule(CyclotomicField, m)
    b = gmodule(AnticNumberField, a)
    push!(Z, b)
  end

  if !minimal_degree
    return Z
  else
    res = GModule[]
    return map(_minimize, Z)
  end
end

function _minimize(V::GModule{<:Oscar.GAPGroup, Generic.FreeModule{nf_elem}})
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
    s = [x for x = s if degree(x[1]) >= d*degree(k)]
    sort!(s, lt = (a,b) -> degree(a[1]) < degree(b[1]))
    for (m, mm) = s
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
  z = irreducible_modules(CyclotomicField, G)
  return [gmodule(QQ, m) for m in z]
end

function irreducible_modules(::ZZRing, G::Oscar.GAPGroup)
  z = irreducible_modules(QQ, G)
  return [gmodule(ZZ, m) for m in z]
end

function gmodule(::typeof(CyclotomicField), C::GModule)
  @assert isa(base_ring(C), QQAbField)
  d = dim(C)
  l = 1
  for g = C.ac
    l = lcm(l, lcm(collect(map_entries(x->Hecke.is_cyclotomic_type(parent(x.data))[2], mat(g)))))
  end
  K = cyclotomic_field(base_ring(C), l)[1]
  F = free_module(K, dim(C))
  if d == 0 
    h = hom(F, F, elem_type(F)[])
    return gmodule(F, group(C), typeof(h)[hom(F, F, map_entries(x->K(x.data), mat(x))) for x = C.ac])
  end
  return gmodule(F, group(C), [hom(F, F, map_entries(x->K(x.data), mat(x))) for x = C.ac])
end

function gmodule(k::Nemo.fpField, C::GModule{PermGroup, GrpAbFinGen})
  q, mq = quo(C.M, characteristic(k))
  s, ms = snf(q)

  r = ngens(s)
  F = free_module(k, r)
  mp = [GrpAbFinGenMap(ms*pseudo_inv(mq)*x*mq*pseudo_inv(ms)) for x= C.ac]
  return gmodule(F, group(C), [hom(F, F, map_entries(k, x.map)) for x = mp])
end

function gmodule(k::Nemo.fpField, mC::Hecke.MapClassGrp)
  return gmodule(k, gmodule(ray_class_field(mC)))
end


import Base: ^
function ^(C::GModule{<:Any, Generic.FreeModule{nf_elem}}, phi::Map{AnticNumberField, AnticNumberField})
  F = free_module(codomain(phi), dim(C))
  return GModule(group(C), [hom(F, F, map_entries(phi, mat(x))) for x = C.ac])
end

function ^(C::GModule{<:Any, T}, h::Map{S, S}) where T <: S where S
  return GModule(group(C), [inv(h)*x*h for x = C.ac])
end

function ^(C::GModule{<:Any, Generic.FreeModule{QQAbElem}}, phi::Map{QQAbField, QQAbField})
  F = free_module(codomain(phi), dim(C))
  return GModule(F, group(C), [hom(F, F, map_entries(phi, mat(x))) for x = C.ac])
end

function gmodule(::QQField, C::GModule{<:Any, Generic.FreeModule{nf_elem}})
  F = free_module(QQ, dim(C)*degree(base_ring(C)))
  return GModule(F, group(C), [hom(F, F, hvcat(dim(C), [representation_matrix(x) for x = transpose(mat(y))]...)) for y = C.ac])
end

gmodule(k::Nemo.fpField, C::GModule{<:Any, Generic.FreeModule{fpFieldElem}}) = C

function _character(C::GModule{<:Any, <:Generic.FreeModule{<:AbstractAlgebra.FieldElem}})
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
    p = preimage(phi, r)
    T = map_word(p, ac; genimgs_inv = iac)
    push!(chr, (c, trace(mat(T))))
  end
  return chr
end

Oscar.character_field(C::GModule{<:Any, <:Generic.FreeModule{QQFieldElem}}) = QQ

function _character_field(C::GModule{<:Any, <:Generic.FreeModule{nf_elem}})
  val = _character(C)
  k, mkK = Hecke.subfield(base_ring(C), [x[2] for x = val])
  return k, mkK
end

function Oscar.character_field(C::GModule{<:Any, <:Generic.FreeModule{nf_elem}})
  return _character_field(C)[1]
end

function Oscar.character(C::GModule{<:Any, <:Generic.FreeModule{nf_elem}})
  chr = _character(C)
  k, mkK = Hecke.subfield(base_ring(C), [x[2] for x = chr])
  A = maximal_abelian_subfield(ClassField, k)
  c = Hecke.norm(conductor(A)[1])
  QQAb = abelian_closure(QQ)[1]
  K = cyclotomic_field(QQAb, Int(c))[1]
  fl, em = is_subfield(k, K)
  return Oscar.group_class_function(group(C), [QQAb(em(preimage(mkK, x[2]))) for x = chr])
end

function Oscar.character(C::GModule{<:Any, <:Generic.FreeModule{QQFieldElem}})
  QQAb = abelian_closure(QQ)[1]
  return Oscar.group_class_function(group(C), [QQAb(x[2]) for x = _character(C)])
end

function Oscar.character(C::GModule{<:Any, <:Generic.FreeModule{<:AbstractAlgebra.FieldElem}})
  return Oscar.group_class_function(group(C), [base_ring(C)(x[2]) for x = _character(C)])
end


function gmodule(k::Nemo.fpField, C::GModule{<:Any, <:Generic.FreeModule{<:FinFieldElem}})
  F = free_module(k, dim(C)*degree(base_ring(C)))
  return GModule(F, group(C), [hom(F, F, hvcat(dim(C), [representation_matrix(x) for x = transpose(mat(y))]...)) for y = C.ac])
end

function Hecke.frobenius(K::FinField, i::Int=1)
  MapFromFunc(x->Hecke.frobenius(x, i), y -> Hecke.frobenius(x, degree(K)-i), K, K)
end

function gmodule_minimal_field(C::GModule{<:Any, <:Generic.FreeModule{fpFieldElem}})
  return C
end

function gmodule_minimal_field(C::GModule{<:Any, <:Generic.FreeModule{<:FinFieldElem}})
  #always over char field
  K =  base_ring(C)
  d = 0
  while d < degree(K)-1
    d += 1
    degree(K) % d == 0 || continue
    k = GF(Int(characteristic(K)), d)
    D = gmodule_over(k, C, do_error = false)
    D === nothing || return D
  end
  return C
end

function gmodule_minimal_field(C::GModule{<:Any, <:Generic.FreeModule{nf_elem}})
  return _minimize(C) 
end

function gmodule_over(k::FinField, C::GModule{<:Any, <:Generic.FreeModule{<:FinFieldElem}}; do_error::Bool = false)
  #mathematically, k needs to contain the character field
  #only works for irreducible modules
  #requires rel cyclic Galois group, not really finite field...
  #
  K = base_ring(C)
  @assert degree(K) != degree(k)
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
  s = frobenius(K, degree(k))
  mkK = embed(k, K)
  os = divexact(degree(K), degree(k))
  hB = hom_base(C, gmodule(C.M, Group(C), 
                      [hom(C.M, C.M, map_entries(s, mat(x))) for x = C.ac]))
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
  return gmodule(F, Group(C), [hom(F, F, map_entries(x -> preimage(mkK, x), Di*mat(x)*D)) for x = C.ac])
  # return C^-1 x C for x = action_gens(C), coerced into k
end

#...now the same for number fields - and non-cyclic fields.
function gmodule_over(em::Map{AnticNumberField, AnticNumberField}, C::GModule{<:Any, <:Generic.FreeModule{nf_elem}}; do_error::Bool = true)
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

function gmodule_over(::QQField, C::GModule{<:Any, <:Generic.FreeModule{nf_elem}}; do_error::Bool = true)
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
  return MapFromFunc(x->N(h(x.data)), y->M(preimage(h, y.data)), M, N)
end

function Oscar.content_ideal(M::MatElem{nf_elem})
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
module of the autmorphism group over the character field.
If `mA` is given, it needs to map the automorphism group over the
character field into the the automorphisms of the base ring.
"""
function factor_set(C::GModule{<:Any, <:Generic.FreeModule{nf_elem}}, mA::Union{Map, Nothing} = nothing)
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

function _two_cocycle(mA::Map, C::GModule{<:Any, <:Generic.FreeModule{nf_elem}}; do_error::Bool = true, two_cycle::Bool = false)
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

  @vprint :MinField 1 "computing un-normalised 1-chain (of matrices)\n"
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
  sigma = Dict{Tuple{PermGroupElem, PermGroupElem}, nf_elem}()
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
  Sigma = Oscar.GrpCoh.CoChain{2,PermGroupElem, MultGrpElem{nf_elem}}(D, Dict(k=>MK(v) for (k,v) = sigma))
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
  @vtime :MinField 2 r = [A*mat(x)*Ai for x = C.ac]
  return r
end

function isone_cochain(X::Dict{<:GAPGroupElem, <:MatElem{nf_elem}}, mA)
  G = domain(mA)
  for g = G
    for h = G
      @assert X[g*h] == map_entries(mA(h), X[g])*X[h]
    end
  end
end

function isone_cochain(X::Dict{<:GAPGroupElem, nf_elem}, mA)
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
        #= if (g*h)(x) = h(g(x)), then the cocycle should be
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
      if cnt > 10 error("s.th. wiered") end
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
function hilbert90_cyclic(A::MatElem{<:FinFieldElem}, s, os::Int)
  #apart from rand, this would also work over a number field
  cnt = 1
  while true
    B = rand(parent(A))
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
      if cnt > 10 && get_assert_level(:MinField) > 1
        error("")
      end
      cnt += 1
    end
  end
end

function gmodule(k::Nemo.fpField, C::GModule{<:Any, Generic.FreeModule{QQFieldElem}})
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

function Gap(C::GModule{<:Any, <:Generic.FreeModule{<:FinFieldElem}}, h=Oscar.iso_oscar_gap(base_ring(C)))
  z = get_attribute(C, :Gap)
  if z !== nothing
    return z
  end
  z = GAP.Globals.GModuleByMats(GAP.Obj([GAP.Obj(map(h, Matrix(mat(x)))) for x = C.ac]), codomain(h))
  set_attribute!(C, :Gap=>z)
  return z
end

function Oscar.is_irreducible(C::GModule{<:Any, <:Generic.FreeModule{<:FinFieldElem}})
  G = Gap(C)
  return GAP.Globals.MTX.IsIrreducible(G)
end

function Oscar.is_absolutely_irreducible(C::GModule{<:Any, <:Generic.FreeModule{<:FinFieldElem}})
  G = Gap(C)
  return GAP.Globals.MTX.IsAbsolutelyIrreducible(G)
end

function is_decomposable(C::GModule{<:Any, <:Generic.FreeModule{<:FinFieldElem}})
  G = Gap(C)
  return !GAP.Globals.MTX.IsIndecomposable(G)
end

"""
    composition_factors_with_multiplicity(C::GModule{<:Any, <:Generic.FreeModule{<:FinFieldElem}})

Return the composition factors of `C` with their frequency.
"""
function Oscar.composition_factors_with_multiplicity(C::GModule{<:Any, <:Generic.FreeModule{<:FinFieldElem}})
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
    indecomposition(C::GModule{<:Any, <:Generic.FreeModule{<:FinFieldElem}})

Return a decomposition of the module `C` into indecomposable summands as a list
of pairs: 
 - a direct indecomposable summand
 - a homomorphism (embedding) of the underlying free modules
"""
function indecomposition(C::GModule{<:Any, <:Generic.FreeModule{<:FinFieldElem}})
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
  h = Oscar.iso_oscar_gap(base_ring(C))
  hb = GAP.Globals.MTX.BasisModuleHomomorphisms(Gap(C, h), Gap(D, h))
  n = length(hb)
  b = [matrix([preimage(h, x[i, j]) for i in 1:GAPWrap.NrRows(x), j in 1:GAPWrap.NrCols(x)]) for x in hb]
#  @show [mat(C.ac[i])*b[1] == b[1]*mat(D.ac[i]) for i=1:length(C.ac)]
  return b
end

"""
  C*T[i] = T[i]*D
on return.

Currently assumes no bad primes.
"""
function hom_base(C::_T, D::_T) where _T <: GModule{<:Any, <:Generic.FreeModule{nf_elem}}
  @assert base_ring(C) == base_ring(D)

  p = Hecke.p_start
  p = 2^10
#  p = 127
  m_in = map(mat, C.ac)
  m_out = map(mat, D.ac)
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
      push!(t, hom_base(z1[i], z2[i]))
    end
    tt = [Hecke.modular_lift([t[i][j] for i=1:length(z1)], me) for j=1:length(t[1])]
    if length(tt) == 0
      return []
    end
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
        fl, s = induce_rational_reconstruction(t, pp)
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

Oscar.nbits(a::QQFieldElem) = nbits(numerator(a)) + nbits(denominator(a))

function hom_base(C::_T, D::_T) where _T <: GModule{<:Any, <:Generic.FreeModule{QQFieldElem}}
  @assert base_ring(C) == base_ring(D)

  p = Hecke.p_start
  p = 2^10
  p = 127
  m_in = map(mat, C.ac)
  m_out = map(mat, D.ac)
  local T
  pp = ZZRingElem(1)
  k = base_ring(C)
  @assert base_ring(m_in[1]) == k
  @assert base_ring(m_in[1]) == k
  @assert k == QQ
  #a heuristic when to try to call reconstruct...
  bt = maximum(maximum(nbits, mat(x)) for x = vcat(C.ac, D.ac)) * dim(C)
  reco = 10
  while true
    p = next_prime(p)
    z1 = gmodule(GF(p), C)
    if C === D
      z2 = z1
    else
      z2 = gmodule(GF(p), D)
    end
    
    t = hom_base(z1, z2)  #TODO: here and elsewhere: get a rref of the hom
                          #base to combine!!!!
                          #(replaced by using fault tolerant lifting)
                          #should use vector reconstruction - once we have it
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
          fl, s = induce_rational_reconstruction(t, pp)
          fl || break
          push!(S, s)
        end
      end
      if nbits(pp) > 1000 && get_assert_level(:MinField) > 1
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

function hom_base(C::_T, D::_T) where _T <: GModule{<:Any, <:Generic.FreeModule{<:QQAbElem}}
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

function gmodule(::QQField, C::GModule{<:Any, <:Generic.FreeModule{ZZRingElem}})
  F = free_module(QQ, dim(C))
  return GModule(group(C), [hom(F, F, map_entries(QQ, mat(x))) for x = C.ac])
end

function hom_base(C::_T, D::_T) where _T <: GModule{<:Any, <:Generic.FreeModule{ZZRingElem}}

  h = hom_base(gmodule(QQ, C), gmodule(QQ, D))
  H = vcat([integral_split(matrix(QQ, 1, dim(C)^2, vec(collect(x))), ZZ)[1] for x = h]...)
  H = Hecke.saturate(H)
  return [matrix(ZZ, dim(C), dim(C), vec(collect(H[i, :]))) for i=1:nrows(H)]
end

function gmodule(::ZZRing, C::GModule{<:Any, <:Generic.FreeModule{QQFieldElem}})
  ma = map(mat, C.ac)
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
#to always get pos. def. one needs a different algorithm
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


function Oscar.gmodule(::Type{GrpAbFinGen}, C::GModule{T, Generic.FreeModule{ZZRingElem}}) where {T <: Oscar.GAPGroup}
  A = free_abelian_group(rank(C.M))
  return Oscar.gmodule(Group(C), [hom(A, A, mat(x)) for x = C.ac])
end

function Oscar.gmodule(::Type{GrpAbFinGen}, C::GModule{T, Generic.FreeModule{FpFieldElem}}) where {T <: Oscar.GAPGroup}
  A = abelian_group([characteristic(base_ring(C)) for i=1:rank(C.M)])
  return Oscar.gmodule(A, Group(C), [hom(A, A, map_entries(lift, mat(x))) for x = C.ac])
end

function Oscar.gmodule(::Type{GrpAbFinGen}, C::GModule{T, Generic.FreeModule{fpFieldElem}}) where {T <: Oscar.GAPGroup}
  A = abelian_group([characteristic(base_ring(C)) for i=1:rank(C.M)])
  return Oscar.gmodule(A, Group(C), [hom(A, A, map_entries(lift, mat(x))) for x = C.ac])
end

function Oscar.abelian_group(M::Generic.FreeModule{fqPolyRepFieldElem})
  k = base_ring(M)
  A = abelian_group([characteristic(k) for i = 1:dim(M)*degree(k)])
  n = degree(k)
  function to_A(m::Generic.FreeModuleElem{fqPolyRepFieldElem})
    a = ZZRingElem[]
    for i=1:dim(M)
      c = m[i]
      for j=0:n-1
        push!(a, coeff(c, j))
      end
    end
    return A(a)
  end
  function to_M(a::GrpAbFinGenElem)
    m = fqPolyRepFieldElem[]
    for i=1:dim(M)
      push!(m, k([a[j] for j=(i-1)*n+1:i*n]))
    end
    return M(m)
  end
  return A, MapFromFunc(to_M, to_A, A, M)
end

function Oscar.group(chi::Oscar.GAPGroupClassFunction)
  return chi.table.GAPGroup
end

function Oscar.gmodule(chi::Oscar.GAPGroupClassFunction)
  f = GAP.Globals.IrreducibleAffordingRepresentation(chi.values)
  K = abelian_closure(QQ)[1]
  g = GAP.Globals.List(GAP.Globals.GeneratorsOfGroup(group(chi).X), f)
  z = map(x->matrix(map(y->map(K, y), g[x])), 1:GAP.Globals.Size(g))
  F = free_module(K, degree(Int, chi))
  return gmodule(group(chi), [hom(F, F, x) for x = z])
end

function Oscar.gmodule(::Type{GrpAbFinGen}, C::GModule{T, Generic.FreeModule{fqPolyRepFieldElem}}) where {T <: Oscar.GAPGroup}
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

function Oscar.simplify(C::GModule{<:Any, <:Generic.FreeModule{QQFieldElem}})
  return gmodule(QQ, Oscar.simplify(gmodule(ZZ, C))[1])
end

function action_matrices(C::GModule{<:Any, <:Generic.FreeModule})
  return map(mat, action(C))
end

function Oscar.simplify(C::GModule{<:Any, <:Generic.FreeModule{ZZRingElem}})
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

function Hecke.induce_crt(a::Generic.MatSpaceElem{nf_elem}, b::Generic.MatSpaceElem{nf_elem}, p::ZZRingElem, q::ZZRingElem)
  c = parent(a)()
  pi = invmod(p, q)
  mul!(pi, pi, p)
  pq = p*q
  z = ZZRingElem(0)

  for i=1:nrows(a)
    for j=1:ncols(a)
      c[i,j] = Hecke.induce_inner_crt(a[i,j], b[i,j], pi, pq, z)
    end
  end
  return c
end

function Hecke.induce_rational_reconstruction(a::Generic.MatSpaceElem{nf_elem}, pg::ZZRingElem)
  c = parent(a)()
  for i=1:nrows(a)
    for j=1:ncols(a)
      fl, c[i,j] = rational_reconstruction(a[i,j], pg)
      fl || return fl, c
    end
  end
  return true, c
end

function Hecke.induce_rational_reconstruction(a::ZZMatrix, pg::ZZRingElem)
  c = zero_matrix(QQ, nrows(a), ncols(a))
  for i=1:nrows(a)
    for j=1:ncols(a)
      fl, n, d = rational_reconstruction(a[i,j], pg, ErrorTolerant = true)
      fl || return fl, c
      c[i,j] = n//d
    end
  end
  return true, c
end

export factor_set
export indecomposition
export irreducible_modules
export is_decomposable

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

export factor_set
export indecomposition
export irreducible_modules
export is_decomposable

