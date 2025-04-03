module Misc
using Oscar
import Base: ==, parent

export relative_field


function primitive_element(a::Vector{QQBarFieldElem})
  pe = a[1]
  f = minpoly(pe)
  Qx = parent(f)
  for i = 2:length(a)
    g = minpoly(a[i])
    f = minpoly(pe)
    k, _ = number_field(f, check = false, cached = false)
    lf = collect(keys(factor(k, g).fac))
    for j = 1:length(lf)
      h = map_coefficients(x->Qx(x)(pe), lf[j])
      if is_zero(h(a[i]))
        d = degree(f) * degree(h)
        mu = 0
        while degree(minpoly(pe+mu*a[i])) != d
          mu += 1
          if mu > 10
            error("too bad")
          end
        end
        pe += mu*a[i]
      end
    end
  end
  return pe
end

function Hecke.number_field(::QQField, a::Vector{QQBarFieldElem}; cached::Bool = false)
  return number_field(QQ, primitive_element(a))
end

function Hecke.number_field(::QQField, a::QQBarFieldElem; cached::Bool = false)
  f = minpoly(a)
  k, b = number_field(f, check = false, cached = cached)
  Qx = parent(k.pol)
  function to_k(x::QQBarFieldElem)
    if x == a
      return b
    end
    f = minpoly(x)
    r = roots(k, f)
    pr = 10
    while true
      C = AcbField(pr)
      CalciumFieldElem = C(a)
      lp = findall(i->contains_zero(Qx(i)(CalciumFieldElem) - C(x)), r)
      if length(lp) == 1
        return r[lp[1]]
      end
      if length(lp) == 0
        error("not in the image")
      end
      pr *= 2
      @assert pr < 2^16
    end
  end
  function to_qqbar(x::AbsSimpleNumFieldElem)
    return Qx(x)(a)
  end
  #TODO: make map canonical?
  # ... and return gen(k) instead?
  return k, MapFromFunc(k, parent(a), to_qqbar, to_k)
end

Base.getindex(::QQField, a::QQBarFieldElem) = number_field(QQ, a)
Base.getindex(::QQField, a::Vector{QQBarFieldElem}) = number_field(QQ, a)

function cyclo_fixed_group_gens(a::AbsSimpleNumFieldElem)
  C = parent(a)
  fl, f = Hecke.is_cyclotomic_type(C)
  @assert fl
  if isone(f)
    return [(1,1)]
  end
  p = first(PrimesSet(1000, -1, f, 1))
  k = GF(p)
  o = k(1)
  lf = factor(f)
  test = [div(f, x) for x = keys(lf.fac)]
  while any(i->isone(o^i), test)
    o = rand(k)^divexact(p-1, f)
  end
  poly_a = numerator(parent(C.pol)(a))
  #the conjugates will be o^j for all j coprime to f
  #the plan is to check conjugates being equal mod p and then
  #verify exactly. Alternatively one could make sure that p is
  #large enough...
  bs = Dict{typeof(o), Set{Int}}()
  co = Set([j for j=1:f-1 if gcd(j, f) == 1])
  all_aut = Set([1])
  gn = Tuple{Int, Int}[]
  while length(co) > 0
    j = pop!(co)
    c = poly_a(o^j)
    if haskey(bs, c)
      d = poly_a(gen(C)^j)
      @assert d == poly_a(gen(C)^first(bs[c]))
      #now we have that conj[i] == conj[j], so 
      #I'd like all autmorphisms mapping zeta^i -> zeta^j
      #given aut: zeta -> zeta^k, so
      #zeta^i -> zeta^(ki) and ki = j mod f
      #so k = j * modinv(i, f)
      #but if aut fixes, then all powers will also fix.
      #that should deal with many conjugates...
      m = Set{Int}()
      for i = bs[c]
        k = (invmod(i, f)*j) % f
        if k in all_aut
#          @show "already have $k"
          continue
        end
        push!(gn, (k, f))
        while true
          n = [(k*x) % f for x = all_aut]
          sz = length(all_aut)
          push!(all_aut, n...)
          push!(m, n...)
          if length(all_aut) == sz
            break
          end
        end
      end
      for i = copy(bs[c])
        for k = m
          ki = (k*i) % f
          if ki in co
            pop!(co, ki)
          end
          push!(bs[c], ki)
        end
      end
    else
      bs[c] = Set([j])
    end
  end
  return gn
end

function cyclo_fixed_group_gens(A::AbstractArray{AbsSimpleNumFieldElem})
  if length(A) == 0
    return [(1,1)]
  end
  G = map(cyclo_fixed_group_gens, A)
  F = lcm([x[1][2] for x = G])
  R, mR = unit_group(quo(ZZ, F)[1])
  Q, mQ = sub(R, gens(R))
  for s = G
    zf = quo(ZZ, s[1][2])[1]
    S, mS = unit_group(zf)
    u, mu = sub(S, [preimage(mS, zf(x[1])) for x = s])
    h = hom(R, S, [preimage(mS, mR(R[i])) for i=1:ngens(R)])
    Q = intersect(Q, preimage(h, u)[1])
  end
  s, ms = snf(Q)
  _, sR = is_subgroup(Q, R)
  Qgen = map(sR, map(ms, gens(s)))
  lf = factor(F)
  for p = keys(lf.fac)
    while true
      G = div(F, p)
      zg = quo(ZZ, G)[1]
      S, mS = unit_group(zg)
      hRS = hom(R, S, [preimage(mS, zg(mR(R[i]))) for i = 1:ngens(R)])
      if length(quo(R, Qgen)[1]) == length(quo(S, map(hRS, Qgen))[1])
        R, mR = S, mS
        Qgen = map(hRS, Qgen)
        F = G
      else
        break
      end
    end
  end
  return [(mR(sR(ms(x))), F) for x = gens(s)]
end


#############################################################################
##
## functions that will eventually get defined in Hecke.jl,
## and then should get removed here

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

function Oscar.cokernel(h::Map)
  return quo(codomain(h), image(h)[1])
end

is_sub_with_data(M::FinGenAbGroup, N::FinGenAbGroup) = is_subgroup(M, N)

function Oscar.direct_product(M::AbstractAlgebra.Module...; task::Symbol = :none)
  D, inj, pro = direct_sum(M...)
  if task == :none
    return D
  elseif task == :both
    return D, pro, inj
  elseif task == :sum
    return D, inj
  elseif task == :prod
    return D, pro
  end
  error("illegal task")
end

function Oscar.id_hom(A::AbstractAlgebra.FPModule)
  return Generic.ModuleHomomorphism(A, A, identity_matrix(base_ring(A), ngens(A)))
end

Oscar.elem_type(::Type{Hecke.NfMorSet{T}}) where {T <: Hecke.LocalField} = Hecke.LocalFieldMor{T, T}
parent(f::Hecke.LocalFieldMor) = Hecke.NfMorSet(domain(f))

function (G::FinGenAbGroup)(x::FinGenAbGroupElem)
  fl, m = is_subgroup(parent(x), G)
  @assert fl
  return m(x)
end

#trivia for QQ
Base.minimum(::Map{QQField, AbsSimpleNumField}, I::Union{Hecke.AbsNumFieldOrderIdeal, Hecke.AbsNumFieldOrderFractionalIdeal}) = minimum(I)*ZZ

Hecke.extend(::Hecke.QQEmb, mp::MapFromFunc{QQField, AbsSimpleNumField}) = complex_embeddings(codomain(mp))

Hecke.restrict(::Hecke.NumFieldEmb, ::Map{QQField, AbsSimpleNumField}) = complex_embeddings(QQ)[1]

#XXX: have a type for an implicit field - in Hecke?
#     add all(?) the other functions to it
function relative_field(m::Map{<:AbstractAlgebra.Field, <:AbstractAlgebra.Field})
  k = domain(m)
  K = codomain(m)
  if k == base_field(K)
    return defining_polynomial(K), Hecke.coordinates, representation_matrix
  end
  @assert base_field(k) == base_field(K)
  kt, t = polynomial_ring(k, cached = false)
  f = defining_polynomial(K)
  Qt = parent(f)
  #the Trager construction, works for extensions of the same field given
  #via primitive elements
  h = gcd(gen(k) - map_coefficients(k, Qt(m(gen(k))), parent = kt), map_coefficients(k, f, parent = kt))
  coordinates = function(x::FieldElem)
    @assert parent(x) == K
    c = collect(Hecke.coefficients(map_coefficients(k, Qt(x), parent = kt) % h))
    c = vcat(c, zeros(k, degree(h)-length(c)))
    return c
  end
  rep_mat = function(x::FieldElem)
    @assert parent(x) == K
    c = map_coefficients(k, Qt(x), parent = kt) % h
    m = collect(Hecke.coefficients(c))
    m = vcat(m, zeros(k, degree(h) - length(m)))
    r = m
    for i in 2:degree(h)
      c = shift_left(c, 1) % h
      m = collect(Hecke.coefficients(c))
      m = vcat(m, zeros(k, degree(h) - length(m)))
      r = hcat(r, m)
    end
    return transpose(matrix(r))
  end
  return h, coordinates, rep_mat
end

Oscar.parent(H::AbstractAlgebra.Generic.ModuleHomomorphism{<:FieldElem}) = Hecke.MapParent(domain(H), codomain(H), "homomorphisms")

function Oscar.hom(F::AbstractAlgebra.FPModule{T}, G::AbstractAlgebra.FPModule{T}) where T
  k = base_ring(F)
  @assert base_ring(G) == k
  H = free_module(k, dim(F)*dim(G))
  return H, MapFromFunc(H, Hecke.MapParent(F, G, "homomorphisms"), x->hom(F, G, matrix(k, dim(F), dim(G), vec(collect(x.v)))), y->H(vec(collect(transpose(matrix(y))))))
end

function Oscar.abelian_group(M::Generic.FreeModule{ZZRingElem})
  A = free_abelian_group(rank(M))
  return A, MapFromFunc(A, M, x->M(x.coeff), y->A(y.v))
end

#TODO: for modern fin. fields as well
function Oscar.abelian_group(M::AbstractAlgebra.FPModule{fqPolyRepFieldElem})
  k = base_ring(M)
  A = abelian_group([characteristic(k) for i = 1:dim(M)*degree(k)])
  n = degree(k)
  function to_A(m::AbstractAlgebra.FPModuleElem{fqPolyRepFieldElem})
    a = ZZRingElem[]
    for i=1:dim(M)
      c = m[i]
      for j=0:n-1
        push!(a, coeff(c, j))
      end
    end
    return A(a)
  end
  function to_M(a::FinGenAbGroupElem)
    m = fqPolyRepFieldElem[]
    for i=1:dim(M)
      push!(m, k([a[j] for j=(i-1)*n+1:i*n]))
    end
    return M(m)
  end
  return A, MapFromFunc(A, M, to_M, to_A)
end

function Hecke.induce_crt(a::Generic.MatSpaceElem{AbsSimpleNumFieldElem}, b::Generic.MatSpaceElem{AbsSimpleNumFieldElem}, p::ZZRingElem, q::ZZRingElem)
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

function Hecke.induce_rational_reconstruction(a::Generic.MatSpaceElem{AbsSimpleNumFieldElem}, pg::ZZRingElem; error_tolerant::Bool = false)
  c = parent(a)()
  for i=1:nrows(a)
    for j=1:ncols(a)
      fl, c[i,j] = rational_reconstruction(a[i,j], pg)#, error_tolerant = error_tolerant)
      fl || return fl, c
    end
  end
  return true, c
end

function Hecke.induce_rational_reconstruction(a::ZZMatrix, pg::ZZRingElem; error_tolerant::Bool = false, unbalanced::Bool = false)
  return Nemo._induce_rational_reconstruction_nosplit(a, pg; error_tolerant, unbalanced)
end

#############################################################################
##
## functions that will eventually get defined in AbstractAlgebra.jl,
## and then should get removed here

Oscar.is_free(M::Generic.DirectSumModule) = all(is_free, M.m)

function Oscar.pseudo_inv(h::Generic.ModuleHomomorphism)
  return MapFromFunc(codomain(h), domain(h), x->preimage(h, x))
end

end # module
using .Misc
export relative_field
