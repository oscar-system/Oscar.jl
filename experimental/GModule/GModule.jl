module GModuleFromGap
using Oscar
using Hecke

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

gmodule(k::Nemo.GaloisField, C::GModule{<:Any, Generic.FreeModule{gfp_elem}}) = C

function Oscar.representation_matrix(a::fq_nmod)
  K = parent(a)
  k = GF(Int(characteristic(K)))
  m = zero_matrix(k, degree(K), degree(K))
  b = basis(K)
  for i=1:degree(K)
    c = a*b[i]
    for j=1:degree(K)
      m[i,j] = coeff(c, j-1)
    end
  end
  return m
end

function _character(C::GModule{<:Any, <:Generic.FreeModule{<:Union{nf_elem, fmpq}}})
  G = group(C)
  phi = GAP.Globals.EpimorphismFromFreeGroup(G.X)
  ac = Oscar.GrpCoh.action(C)
  iac = Oscar.GrpCoh.inv_action(C)

  n = dim(C)
  K = base_ring(C)

  chr = []
  for c = conjugacy_classes(G)
    r = representative(c)
    if isone(r)
      push!(chr, (c, K(n)))
      continue
    end
    p = GAP.Globals.PreImagesRepresentative(phi, r.X)
    w = map(Int, GAP.Globals.LetterRepAssocWord(p))
    T = w[1]<0 ? iac[-w[1]] : ac[w[1]]
    for i=2:length(w)
      T = T*(w[i]<0 ? iac[-w[i]] : ac[w[i]])
    end
    push!(chr, (c, trace(mat(T))))
  end
  return chr
end

character_field(C::GModule{<:Any, <:Generic.FreeModule{fmpq}}) = QQ

function character_field(C::GModule{<:Any, <:Generic.FreeModule{nf_elem}})
  val = _character(C)
  k, mkK = Hecke.subfield(base_ring(C), [x[2] for x = val])
  return k
end

function character(C::GModule{<:Any, <:Generic.FreeModule{nf_elem}})
  chr = _character(C)
  k, mkK = Hecke.subfield(base_ring(C), [x[2] for x = chr])
  A = maximal_abelian_subfield(ClassField, k)
  c = Hecke.norm(conductor(A)[1])
  K = cyclotomic_field(Int(c))[1]
  fl, em = issubfield(k, K)
  return [(x[1], em(preimage(mkK, x[2]))) for x = chr]
end

function character(C::GModule{<:Any, <:Generic.FreeModule{fmpq}})
  return _character(C)
end

function gmodule(k::Nemo.GaloisField, C::GModule{<:Any, <:Generic.FreeModule{<:FinFieldElem}})
  F = free_module(k, dim(C)*degree(base_ring(C)))
  return GModule(F, group(C), [hom(F, F, hvcat(dim(C), [representation_matrix(x) for x = transpose(mat(y))]...)) for y = C.ac])
end

function Hecke.frobenius(K::FinField, i::Int=1)
  MapFromFunc(x->Hecke.frobenius(x, i), y -> Hecke.frobenius(x, degree(K)-i), K, K)
end

function gmodule_minimal_field(C::GModule{<:Any, <:Generic.FreeModule{gfp_elem}})
  return C
end

function gmodule_minimal_field(C::GModule{<:Any, <:Generic.FreeModule{<:FinFieldElem}})
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
  # (since rho is irreducible and thus the matrix unque up to scalar)
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
  alpha = _norm_equation(K, preimage(mkK, lambda))
  B *= inv(alpha)
  @hassert :MinField 2 isone(norm(B, s, os))
  D = hilbert90_cyclic(B, s, os)
  Di = inv(D)
  F = free_module(k, dim(C))
  return gmodule(F, Group(C), [hom(F, F, map_entries(x -> preimage(mkK, x), Di*mat(x)*D)) for x = C.ac])
  # return C^-1 x C for x = action_gens(C), coerced into k
end

#...now the same for number fields - and non-cyclic fields.
function gmodule_over(em::Map{AnticNumberField, AnticNumberField}, C::GModule{<:Any, <:Generic.FreeModule{nf_elem}}; do_error::Bool = false)
  K = base_ring(C)
  k = domain(em)
  @assert codomain(em) == K
  gk = em(gen(k))

  A, mA = automorphism_group(PermGroup, K)
  s, ms = sub(A, [a for a = A if mA(a)(gk) == gk])
  ac =  _two_cocycle(ms*mA, C, do_error = do_error)  
  F = free_module(k, dim(C))
  return gmodule(F, group(C), [hom(F, F, map_entries(t->preimage(em, t), x)) for x = ac])
end

function gmodule_over(::FlintRationalField, C::GModule{<:Any, <:Generic.FreeModule{nf_elem}}; do_error::Bool = false)
  K = base_ring(C)
  A, mA = automorphism_group(PermGroup, K)
  ac = _two_cocycle(mA, C, do_error = do_error)  
  F = free_module(QQ, dim(C))
  return gmodule(F, group(C), [hom(F, F, map_entries(QQ, x)) for x = ac])
end

function Oscar.hom(M::MultGrp, N::MultGrp, h::Map)
  return MapFromFunc(x->N(h(x.data)), y->M(preimage(h, N.data)), M, N)
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

function _two_cocycle(mA::Map, C::GModule{<:Any, <:Generic.FreeModule{nf_elem}}; do_error::Bool = false)
  G = domain(mA)
  K = base_ring(C)

  homs = []
  @vprint :MinField 1 "Gathering Galois images of the generators...\n"
  for g = gens(G)
    @vprint :MinField 2 "gen: $g\n"
    @vtime :MinField 2 hb = hom_base(C, C^mA(g))
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
  # pairs: (g, X_g) with operation (g, X_g)(h, X_h) = (gh, X_h*X_g^h)
  @vtime :MinField 2 
    c = closure([(gen(G, i), homs[i]) for i=1:ngens(G)], 
              (a, b) -> (a[1]*b[1], b[2]*map_entries(mA(b[1]), a[2])),
              (one(G), I),
              eq = (a,b) -> a[1] == b[1])
  X = Dict(x[1] => x[2] for x = c)
  X[one(G)] = I

  #now we need a 2-cycle:
  #X[g] X[h] = sigma(g, h) X[gh] should hold...

  @vprint :MinField 1 "now the 2-cocycle (scalars)\n"
  MK = MultGrp(K)
  sigma = Dict{Tuple{PermGroupElem, PermGroupElem}, MultGrpElem{nf_elem}}()
  for g = G
    for h = G
      if isone(g)
        sigma[(g, h)] = MK(one(K))
      elseif isone(h)
        sigma[(g, h)] = MK(one(K))
      else
        lf = findfirst(x->!iszero(x), X[g*h])
        sigma[(g, h)] = MK(X[g*h][lf]//(X[h]*map_entries(mA(h), X[g]))[lf])
      end
    end
  end

  @vprint :MinField 1 "test for co-boundary\n"
  D = gmodule(G, [hom(MK, MK, mA(x)) for x = gens(G)])
  Sigma = Oscar.GrpCoh.CoChain{2,PermGroupElem, MultGrpElem{nf_elem}}(D, sigma)
  @vtime :MinField 2 fl, cb = Oscar.GrpCoh.iscoboundary(Sigma)

  if !fl
    do_error || return nothing
    error("field too small")
  end
  for g = G
    X[g] *= inv(cb(g).data)
  end

  #now X should be in H^1(G, Gl(n, K)) which is trivial
  #hence a Hilbert-90 should find A s.th. A^(1-g) = X[g] for all g

  @vprint :MinField 1 "calling Hilbert-90 on matrices\n"
  @vtime :MinField 2 A, Ai = hilbert90_generic(X, mA)
  c = content_ideal(A)
  d = Hecke.short_elem(inv(c))
  A *= d
  Ai *= inv(d)

  @vprint :MinField 1 "conjugating the generators\n"
  @vtime :MinField 2 r = [Ai*mat(x)*A for x = C.ac]
  return r
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
      fl = isinvertible(Y)
      fl && break
      cnt += 1
      if cnt > 10 error("s.th. wiered") end
    end
    S = sum(v*map_entries(mA(g), Y) for (g,v) = X)
    fl, Si = isinvertible_with_inverse(S)
    fl && return S, Si
  end
end

"""
Hunt for b s.th. N(b) == a
"""
#TODO conflicts with the "same" algo in Hecke - where it does not
#     work due to a missing function or so.
function _norm_equation(K::FinField, a::FinFieldElem)
  # a in K, k = fix(K, <s>)
  # we want irr. poly over k of degree(K:k) with constant term \pm a
  # then a root in K has the norm...
  k = parent(a)
  fkK = embed(k, K)
  os = divexact(degree(K), degree(k))
  if isodd(os)
    a = -a
  end
  kt, t = PolynomialRing(k, cached = false)
  while true
    f = t^os + a + sum(t^rand(1:os-1)*rand(k) for i=1:rand(1:os-1))
    isirreducible(f) || continue
    r = roots(map_coefficients(fkK, f))[1]
    return r
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
  #apart form rand, this would also work over a number field
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

function Gap(C::GModule{<:Any, <:Generic.FreeModule{<:FinFieldElem}}, h=Oscar.iso_oscar_gap(base_ring(C)))
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

function Oscar.gmodule(::Type{GrpAbFinGen}, C::GModule{T, Generic.FreeModule{gfp_elem}}) where {T <: Oscar.GAPGroup}
  A = abelian_group([characteristic(base_ring(C)) for i=1:rank(C.M)])
  return Oscar.gmodule(A, Group(C), [hom(A, A, map_entries(lift, mat(x))) for x = C.ac])
end

function Oscar.abelian_group(M::Generic.FreeModule{fq_nmod})
  k = base_ring(M)
  A = abelian_group([characteristic(k) for i = 1:dim(M)*degree(k)])
  n = degree(k)
  function to_A(m::Generic.FreeModuleElem{fq_nmod})
    a = fmpz[]
    for i=1:dim(M)
      c = m[i]
      for j=0:n-1
        push!(a, coeff(c, j))
      end
    end
    return A(a)
  end
  function to_M(a::GrpAbFinGenElem)
    m = fq_nmod[]
    for i=1:dim(M)
      push!(m, k([a[j] for j=(i-1)*n+1:i*n]))
    end
    return M(m)
  end
  return A, MapFromFunc(to_M, to_A, A, M)
end

function Oscar.gmodule(::Type{GrpAbFinGen}, C::GModule{T, Generic.FreeModule{fq_nmod}}) where {T <: Oscar.GAPGroup}
  k = base_ring(C)
  A, mA = abelian_group(C.M)

  return Oscar.gmodule(A, Group(C), [hom(A, A, [preimage(mA, x(mA(a))) for a = gens(A)]) for x = C.ac])
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

export irreducible_modules, isabsolutely_irreducible, isdecomposable

module RepPc
using Oscar

export coimage

Base.pairs(M::MatElem) = Base.pairs(IndexCartesian(), M)
Base.pairs(::IndexCartesian, M::MatElem) = Base.Iterators.Pairs(M, CartesianIndices(axes(M)))

function Hecke.roots(a::FinFieldElem, i::Int)
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
"""
  For K a finite field, Q, a number field or Qab, find all
abs. irred. representations of G.

Note: the reps are NOT neccessarily over the smallest field.

Note: the field is NOT extended - but it throws an error if it was too small.

Implements: Brueckner, Chap 1.2.3
"""
function reps(K, G::Oscar.GAPGroup)
  isfinite(G) || error("the group is not finite")
  if order(G) == 1
    F = free_module(K, 1)
    h = hom(F, F, [F[1]])
    return [gmodule(F, G, typeof(h)[])]
  end

  pcgs = GAP.Globals.Pcgs(G.X)
  pcgs == GAP.Globals.fail && error("the group is not polycyclic")

  gG = [Oscar.group_element(G, x) for x = pcgs]
  s, ms = sub(G, [gG[end]])
  o = Int(order(s))
  @assert isprime(o)
  z = roots(K(1), o)
  @assert characteristic(K) == o || length(z) == o
  F = free_module(K, 1)
  R = [gmodule(F, s, [hom(F, F, [r*F[1]])]) for r = z]
  @hassert :BruecknerSQ 2 Oscar.GrpCoh.isconsistent(R[1])

  for i=length(gG)-1:-1:1
    h = gG[i]
    ns, mns = sub(G, gG[i:end])
    @assert mns(ns[1]) == h
    p = Int(divexact(order(ns), order(s)))
    @assert isprime(p)
    new_R = []
    todo = trues(length(R)) # which entries in `R` have to be handled
    #TODO: use extend below
    for pos in 1:length(R)
      if todo[pos]
        r = R[pos]
        F = r.M
        @assert group(r) == s
        rh = gmodule(group(r), [action(r, preimage(ms, x^h)) for x = gens(s)])
        @hassert :BruecknerSQ 2 Oscar.GrpCoh.isconsistent(rh)
        l = Oscar.GModuleFromGap.hom_base(r, rh)
        @assert length(l) <= 1
        Y = mat(action(r, preimage(ms, h^p)))
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
          rt = roots(C, p)
          @assert characteristic(K) == p || length(rt) == p
          Y = r.ac
          for x = rt
            nw = gmodule(F, ns,  vcat([hom(F, F, x*X)], Y))
            @hassert :BruecknerSQ 2 Oscar.GrpCoh.isconsistent(nw)
            push!(new_R, nw)
          end
        else #need to extend dim
          n = dim(r)
          F = free_module(K, dim(r)*p)

          # a block permutation matrix for the element `h`
          z = zero_matrix(K, dim(F), dim(F))
          z[1:n,(p-1)*n+1:end] = inv(Y)
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
              m = mat(Y)
              z[(j-1)*n+1:j*n, (j-1)*n+1:j*n] = m
              push!(conjreps[j], hom(M, M, m))
              g = preimage(ms, ms(g)^h)
            end
            push!(md, hom(F, F, z))
          end

          # Find the positions of the equiv. classes of the `h`-conjugate
          # representations, we need not deal with them later on
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
          @hassert :BruecknerSQ 2 Oscar.GrpCoh.isconsistent(new_R[end])
        end
      end
    end
    s, ms = ns, mns
    R = new_R
  end
  return R
end

#TODO: do this properly, eventually...
#      hnf is a rref. maybe we could use a ref s.w.?
#      in this file, all(?) hnf are actually over GF(p), so 
#      maybe use this as well?
function Nemo._hnf(x::fmpz_mat)
  if nrows(x) * ncols(x) > 100
    s = sparse_matrix(x)
    if sparsity(s) > 0.7
      return matrix(Hecke.hnf(s))
    end
  end
  return Nemo.__hnf(x) # ist die original Nemo flint hnf
end

function Nemo._hnf_with_transform(x::fmpz_mat)
  if nrows(x) * ncols(x) > 100
    s = sparse_matrix(x)
    if sparsity(s) > 0.7
      s = hcat(s, identity_matrix(SMat, ZZ, nrows(x)))
      m = matrix(Hecke.hnf(s))
      return m[:, 1:ncols(x)], m[:, ncols(x)+1:end]
    end
  end
  return Nemo.__hnf_with_transform(x) # ist die original Nemo flint hnf
end


function extend(C::GModule, m::Map)
  #N acts and is normal in <h, N> 
  #m injects N into <h, N>  (at least h, maybe more)
  #h = gen(domain(m), 1)
  #h has order p in <h, N>/N
  #Satz 10 in Brueckner, Chap 1.2.3
  #TODO: should be used above in reps!

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
    #TODO: see above in reps, this is wrong.
    #TODO: fuse
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

"""
Brueckner Chap 1.3.1

Given 
  mp: G ->> Q

Find a set of primes suth that are any irreducible F_p module M
s.th. there is an epimorphism of G onto the extension of Q by M,
the p is in the set.
"""
function find_primes(mp::Map{FPGroup, PcGroup})
  G = domain(mp)
  Q = codomain(mp)
  I = irreducible_modules(ZZ, Q) 
  lp = Set(collect(keys(factor(order(Q)).fac)))
  for i = I
    ib = gmodule(i.M, G, [action(i, mp(g)) for g = gens(G)])
    ia = gmodule(GrpAbFinGen, ib)
    a, b = Oscar.GrpCoh.H_one_maps(ia)
    da = Oscar.dual(a)
    db = Oscar.dual(b)
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
    q = quo(kernel(da)[1], image(db)[1])[1]
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

If neccessary, the prime(s) p that can be used are computed as well.
"""
function brueckner(mQ::Map{FPGroup, PcGroup}; primes::Vector=[])
  Q = codomain(mQ)
  G = domain(mQ)
  @vprint :BruecknerSQ 1 "lifting $mQ using SQ\n"
  if length(primes) == 0
    @vprint :BruecknerSQ 1 "primes not provided, searching...\n"
    lp = find_primes(mQ)
  else
    lp = map(fmpz, primes)
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
      append!(allR, [x for x in l])# if issurjective(x)])
    end
  end
  return allR
end

Base.getindex(M::AbstractAlgebra.FPModule, i::Int) = i==0 ? zero(M) : gens(M)[i]
Oscar.gen(M::AbstractAlgebra.FPModule, i::Int) = M[i]

Oscar.isfree(M::Generic.FreeModule) = true
Oscar.isfree(M::Generic.DirectSumModule) = all(isfree, M.m)

function Oscar.dual(h::Map{GrpAbFinGen, GrpAbFinGen})
  A = domain(h)
  B = codomain(h)
  @assert isfree(A) && isfree(B)
  return hom(B, A, transpose(h.map))
end

function Oscar.dual(h::Map{<:AbstractAlgebra.FPModule{fmpz}, <:AbstractAlgebra.FPModule{fmpz}})
  A = domain(h)
  B = codomain(h)
  @assert isfree(A) && isfree(B)
  return hom(B, A, transpose(mat(h)))
end

function coimage(h::Map)
  return quo(domain(h), kernel(h)[1])
end

function Base.iterate(M::Generic.Submodule{<:FinFieldElem})
  k = base_ring(M)
  if dim(M) == 0
    return zero(M), iterate([1])
  end
  p = Base.Iterators.ProductIterator(Tuple([k for i=1:dim(M)]))
  f = iterate(p)
  return M(elem_type(k)[f[1][i] for i=1:dim(M)]), (f[2], p)
end

function Base.iterate(::AbstractAlgebra.Generic.Submodule{fq_nmod}, ::Tuple{Int64, Int64})
  return nothing
end

function Base.iterate(M::Generic.Submodule{<:FinFieldElem}, st::Tuple{<:Tuple, <:Base.Iterators.ProductIterator})
  n = iterate(st[2], st[1])
  if n === nothing
    return n
  end
  return M(elem_type(base_ring(M))[n[1][i] for i=1:dim(M)]), (n[2], st[2])
end

"""
  mp: G ->> Q
  C a F_p[Q]-module
  Find all extensions of Q my C s.th. mp can be lifted to an epi.
"""
function lift(C::GModule, mp::Map)
  #m: G->group(C)
  #compute all(?) of H^2 that will descibe groups s.th. m can be lifted to

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
  z = get_attribute(C, :H_two)[1]

  seen = Set{Tuple{elem_type(D), elem_type(codomain(mH2))}}()
  #TODO: the projection maps seem to be rather slow - in particular
  #      as they SHOULD be trivial...
  for x = k
    epi = pDE[1](mk(x)) #the map
    chn = pDE[2](mk(x)) #the tail data
    if (epi,mH2(chn)) in seen
      continue
    else
      push!(seen, (epi, mH2(chn)))
    end
    #TODO: not all "chn" yield distinct groups - the factoring by the 
    #      co-boundaries is missing
    #      not all "epi" are epi, ie. surjective. The part of the thm
    #      is missing...
    # (Thm 15, part b & c) (and the weird lemma)
    @hassert :BruecknerSQ 2 all(x->all(y->sc(x, y)(chn) == last_c(x, y), gens(N)), gens(N))
    @hassert :BruecknerSQ 2 preimage(z, z(chn)) == chn
    GG, GGinj, GGpro, GMtoGG = Oscar.GrpCoh.extension(PcGroup, z(chn))
    if get_assert_level(:BruecknerSQ) > 1
      _GG, _ = Oscar.GrpCoh.extension(z(chn))
      @assert isisomorphic(GG, _GG)[1]
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

    h = hom(G, GG, gens(G), [GMtoGG(reduce(gen(G, i)), pro[i](epi)) for i=1:ngens(G)])
    if !issurjective(h)
      @show :darn
      continue
    else
      @show :bingo
    end
    push!(allG, h)
  end
  return allG
end

end #module RepPc

using .RepPc

export coimage
