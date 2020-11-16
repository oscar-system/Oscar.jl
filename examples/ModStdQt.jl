module MPolyInterpolate

using Oscar
import Oscar.Nemo
import Oscar.Hecke

function Oscar.evaluate(f::FracElem{<:MPolyElem}, v::Vector; ErrorTolerant::Bool = false)
  n = evaluate(numerator(f), v)
  d = evaluate(denominator(f), v)
  if iszero(d) 
    if ErrorTolerant
      return n//1
    else
      throw(DivideError())
    end
  end
  return n//d
end

Oscar.normalise(f::fmpq_poly, ::Int64) = error("no normlise") #length(f)
Oscar.set_length!(f::fmpq_poly, ::Int64) = error("no sert_length") #f

function points(p::Vector{T}, j::Int, z::Vector{T}, s::Vector{T}) where {T}

  pp = map(fmpz, p)
  zz = map(fmpz, z)
  ss = map(fmpz, s)

  return [[[pp[i]^l*zz[k]+ss[i] for i=1:length(pp)] for k=1:length(z)] for l=0:j]
end

function eval(f::FracElem{fmpq_mpoly}, pt; ErrorTolerant::Bool = false)
  return [[evaluate(f, x, ErrorTolerant= ErrorTolerant) for x = y] for y = pt]
end

mutable struct InterpolateCtx{T}
  C::Hecke.crt_env{T}
  function InterpolateCtx(a::Vector{S}, parent = PolynomialRing(parent(a[1]), cached = false)[1]) where {S}
    bt = parent
    t = gen(bt)
    r = new{elem_type(bt)}()
    r.C = Hecke.crt_env([t - i for i = a])
    return r
  end
end

function Oscar.base_ring(I::InterpolateCtx)
  return base_ring(base_ring(I.C))
end
function Oscar.base_ring(C::Hecke.crt_env)
  return parent(C.pr[1])
end

function Oscar.interpolate(v::Vector, I::InterpolateCtx)
  bt = base_ring(I.C)
  return crt([bt(i) for i = v], I.C)
end

function Oscar.prod(C::Hecke.crt_env)
  return C.pr[end]
end

#does not exist in general... why
Hecke.rem!(a::T, b::T, c::T) where {T <: RingElem} = rem(b, c)

function Oscar.rational_reconstruction(f::PolyElem, I::InterpolateCtx; ErrorTolerant::Bool = false)
  return rational_reconstruction(f, prod(I.C), ErrorTolerant = ErrorTolerant)
end

mutable struct MPolyPt{T}
  pt_z::Array{T, 1}
  pt_p::Array{T, 1}
  pt_s::Array{T, 1}
  j::Int

  function MPolyPt(n::Int)
    p = Array{Int, 1}()
    for q = PrimesSet(11, -1)
      push!(p, q)
      if length(p) == n
        break
      end
    end
    return new{Int}(Int[], p, Int[rand(-100:100) for i=1:n], 0)
  end
  function MPolyPt(lp::Vector{Int}, z::AbstractVector, s::Vector{Int}, j::Int)
    r = new{Int}()
    r.pt_z = collect(z)
    r.pt_p = lp
    r.pt_s = s
    r.j = j
    return r
  end
end

function eval(f::FracElem{fmpq_mpoly}, P::MPolyPt; ErrorTolerant::Bool = false)
  return eval(f, points(P.pt_p, P.j, P.pt_z, P.pt_s), ErrorTolerant = ErrorTolerant)
end

function Base.push!(pt::MPolyPt, p::Int)
  if p in pt.pt_z
    error("point already there")
  end
  push!(pt.pt_z, p)
end

mutable struct MPolyInterpolateCtx{T}
  R#::parent_type(T)
  parent#::FracField{T}

  pt::MPolyPt
  I::InterpolateCtx#{T}
  status::Symbol

  function MPolyInterpolateCtx(F::FracField{<:MPolyElem{S}}, pt::MPolyPt) where {S}
    r = new{S}()
    r.parent = F
    r.R = base_ring(F)
    r.pt = pt
    r.I = InterpolateCtx(map(base_ring(r.R), pt.pt_z))
    return r
  end
end

function set_status!(M::MPolyInterpolateCtx, s::Symbol)
  @show s
  M.status = s
end

#Ben-Or Tiwari...
function Oscar.interpolate(val::Array{Array{T, 1}, 1}, M::MPolyInterpolateCtx) where {T}

#  nd = Array{Tuple{Bool, fmpq_poly, fmpq_poly}, 1}()
  nd = []
  for x = val
    f = interpolate(x, M.I)
    mu = rational_reconstruction(f, M.I, ErrorTolerant = true)
    if !mu[1]
      set_status!(M, :univariate_failed)
      return false, zero(M.R) # more z
    end
    t = inv(trail(mu[3]))
    push!(nd, (true, mu[2]*t, mu[3]*t))
  end
  Qx = parent(nd[1][2])
  @assert parent(nd[1][3]) == Qx

  R = [M.R(), M.R()]
  for ii = 1:2
    cor = elem_type(Qx)[zero(Qx) for i=1:length(nd)]
    for i=maximum(degree(x[ii+1]) for x = nd):-1:0
      bm = berlekamp_massey([coeff(nd[x][ii+1]-cor[x], i) for x = 1:length(nd)])
      if !bm[1]
        set_status!(M, :BM_failed)
        return false, zero(M.R)
      end
      if isconstant(bm[2])
        continue
      end

      r = roots(bm[2])
      if length(r) != degree(bm[2]) 
        set_status!(M, :no_roots)
        return false, zero(M.R)
      end
      if !all(isinteger, r) || any(iszero, r)
        set_status!(M, :wrong_roots)
        return false, zero(M.R)
      end

      V = Vandermonde(r)
      s = solve_left(V, matrix(QQ, 1, length(r), [coeff(nd[x][ii+1]-cor[x], i) for x=1:length(r)]))
      G = MPolyBuildCtx(M.R)
      for i=1:length(r)
        push_term!(G, s[1,i], [valuation(r[i], p) for p = M.pt.pt_p])
      end
      g = finish(G)
      R[ii] += g
      for j=0:length(nd)-1
        t = evaluate(g, [M.pt.pt_p[l]^j*gen(Qx) + M.pt.pt_s[l] for l=1:length(M.pt.pt_p)])
  #      @show t, lead(t), lead(nd[j+1][2] - cor[j+1])
        cor[j+1] += t
      end
  #    @show cor
    end
  end
  if iszero(R[2])
    set_status!(M, :bad)
    return false, zero(M.R)
  end
  set_status!(M, :OK)
  return true, R[1]//R[2]
end

struct Vals{T}
  v::Array{Array{T, 1}, 1}
end

function Base.setindex!(v::Vals{T}, c::T, i::Int, j::Int) where {T}
  if i > length(v.v)
    push!(v.v, zeros(parent(v.v[1][1]), length(v.v[1])))
  end
  if j > length(v.v[i])
    while length(v.v[i]) < j-1
      push!(v.v[i], parent(v.v[1][1])(0))
    end
    @assert j == length(v.v[i])+1
    push!(v.v[i], c)
    return
  end
  v.v[i][j] = c
end

function Groebner_basis(I::Oscar.MPolyIdeal{<:Generic.MPoly{<:Generic.Frac{<:MPolyElem{T}}}}) where {T <: RingElem}
  S = base_ring(I)
  R = base_ring(S)
  n = ngens(base_ring(R))
  P = MPolyPt(n)
  Q = base_ring(base_ring(R))
  Qy, _ = PolynomialRing(Q, ngens(S))
  frst = true
  lst = []
  do_z = true
  local lst
  while true
    if do_z
      push!(P, length(P.pt_z)+1)
    else
      P.j += 1
    end

    idx = 1
    while true
      gens_J = elem_type(Qy)[]
      if do_z
        if idx > P.j+1
          break
        end
        pt = T[Q(P.pt_p[i])^(idx-1) * Q(P.pt_z[end]) + Q(P.pt_s[i]) for i=1:n]
#        println( "z: j = ", idx -1, ", z = ", P.pt_z[end])
      else
        if idx > length(P.pt_z)
          break
        end
        pt = T[Q(P.pt_p[i])^P.j * Q(P.pt_z[idx]) + Q(P.pt_s[i]) for i=1:n]
#        println("P: j = ", P.j, ", z = ", P.pt_z[idx])
      end
      idx += 1
      for g = gens(I)
        G = MPolyBuildCtx(Qy)
        for (c, e) = Base.Iterators.zip(Generic.MPolyCoeffs(g), Generic.MPolyExponentVectors(g))
          push_term!(G, evaluate(c, pt, ErrorTolerant = true), e)
        end
        push!(gens_J, finish(G))
      end
      J = ideal(Qy, gens_J)
      gJ = groebner_basis(J, :degrevlex, complete_reduction = true)
      if frst
        lst = []
        for _g = gJ
          g = inv(lead(_g))*_g
          f = []
          for (c, e) = Base.Iterators.zip(Generic.MPolyCoeffs(g), Generic.MPolyExponentVectors(g))
            push!(f, (e, Vals(Array{T, 1}[[c]])))
          end
          push!(lst, f)
        end
        frst = false
      else
        for ig = 1:length(gJ)
          g = gJ[ig]
          g *= inv(lead(g))
          jg = 1
          for (c, e) = Base.Iterators.zip(Generic.MPolyCoeffs(g), Generic.MPolyExponentVectors(g))
            @assert lst[ig][jg][1] == e #TODO: sort and match
            if do_z
              lst[ig][jg][2][idx-1, length(P.pt_z)] = c
            else
              lst[ig][jg][2][P.j+1, idx-1] = c
            end
            jg += 1
          end
        end
      end
    end
    
    res = []
    M = MPolyInterpolateCtx(R, P)
    local fl = true
    for F = lst
      f = MPolyBuildCtx(S)
      for C = F
        fl, c = interpolate(C[2].v, M)
        if fl
          push_term!(f, c, C[1])
        else
          break
        end
      end
      fl || break
      push!(res, finish(f))
    end
    fl && return res

    do_z = !do_z
    if length(P.pt_z) > 5 || P.j > 5
      break
    end
  end
  return lst, P
end

###############################################################
#
# disc_log for finite fields
# TODO
# extend to non-prime
# use the Conway property if applicable
# store the dlog data on the field to re-use
# do more intelligent (sub-exponential) algorithms
# (maybe not all, but at least some?)
# see Hecke/src/Misc/UnitsModM: fro BS-GS, cache the BS Array
# sort the types...
###############################################################
#TODO: consolidate with Hecke: isprimitive_root?
function Oscar.isprimitive(a::Nemo.gfp_elem, f::Fac{fmpz} = factor(length(parent(a))-1))
  n = length(parent(a))-1
  for p = keys(f.fac)
    if a^divexact(n, p) == 1
      return false
    end
  end
  return true
end

function generator(K::Nemo.GaloisField, f::Fac{fmpz} = factor(length(K)-1))
  a = rand(K)
  while !isprimitive(a, f)
    a = rand(K)
  end
  return a
end

function disc_log(a::Nemo.gfp_elem, b::Nemo.gfp_elem)
  return disc_log_bs_gs(a, b, fmpz(length(parent(a))-1))
end

mutable struct DiscLogCtx
  K::Nemo.GaloisField
  f::Fac{fmpz}
  g::Nemo.gfp_elem
  function DiscLogCtx(g::Nemo.gfp_elem, f::Fac{fmpz} = factor(length(parent(g))-1))
    r = new()
    r.K = parent(g)
    r.g = g
    r.f = f
    return r
  end
end

function disc_log(A::DiscLogCtx, a::Nemo.gfp_elem)
  res = Array{Tuple{fmpz, fmpz}, 1}()
#  return Hecke.disc_log_bs_gs(A.g, a, size(parent(a))-1)
  Pm1 = size(parent(a))-1
  for (p, k) = A.f.fac
    r = divexact(Pm1, p^k)
    push!(res, (p^k, Hecke.disc_log_ph(A.g^r, a^r, fmpz(p), k)))
  end
  return crt([x[2] for x = res], [x[1] for x= res])
end
#####################################################################
function ben_or(lp::Array{Int, 1}, a::Array{fmpz, 1}, P::fmpz)
  #tries to write a as a product of powers of lp

  K = GF(Int(P))
  R = quo(ZZ, length(K)-1)[1]
  f = factor(P-1)
  g = generator(K, f)
  A = DiscLogCtx(g, f)

  dl = [disc_log(A, K(x)) for x = lp]
  @show ddl = [disc_log(A, K(aa)) for aa = a]
  z = ([identity_matrix(ZZ, length(lp)) matrix(ZZ, length(lp), 1, dl); zero_matrix(ZZ, 1, length(lp)) matrix(ZZ, 1, 1, [P-1])])
  return z, ddl
end
#= external memory

B = {p_1, ..., p_k}, m = prod p_i^n_i
task: find n_i

for any other prime q we have p_i = g^x_i mod q, m = g^y mod q, so

  sum n_i x_i = y mod q-1
      equivalently
  prod p_i^n_i = m mod q

furthermore, n_i are "small" since fixed.

suppose we have 
  prod p_i^r_i = 1 mod q_j, ie. for several j, then
    (CRT) prod q_j | (prod p_i^r_i) - 1
  so, in particular prod p_i^r_i > prod q_j
  thus
  sum r_i log_2 p_i > sum log q_j
  or
  sum log q_j < sum r_i log p_i <= |r_i, i|_2 |log p_i, i|_2
  via Cauchy-Schwarz, so this gives a LOWER BOUND on |r_i, i|_2

  fix q_i, then { r_i | q_j | (prod p_i^r_i) -1 } forms a lattice,
  we have lower bounds on the minimum (above)

  Thus LLL can find the n_i from above - if the minimum is > 2*n_i I'd guess.

  So lets add more primes!
=#

#############################################################################
# some special lin. alg
#############################################################################
struct Vandermonde{T} <: MatElem{T}
  val::Vector{T}
end

Base.size(V::Vandermonde) = (length(V.val), V.val)
Base.getindex(V::Vandermonde, i::Int, j::Int) = V.val[i]^(j-1)
Oscar.nrows(V::Vandermonde) = size(V)[1]
Oscar.ncols(V::Vandermonde) = size(V)[1]
Oscar.base_ring(V::Vandermonde) = parent(V.val[1])

struct VandermondeT{T} <: MatElem{T}
  V::Vandermonde{T}
end

Base.size(V::VandermondeT) = size(V.V)
Base.getindex(V::VandermondeT, i::Int, j::Int) = V.V[j,i]
Oscar.nrows(V::VandermondeT) = size(V)[1]
Oscar.ncols(V::VandermondeT) = size(V)[1]
Oscar.base_ring(V::VandermondeT) = base_ring(V.V)

Base.transpose(V::Vandermonde) = VandermondeT(V)
Base.transpose(V::VandermondeT) = V.V

struct Hankel{T} <: MatElem{T}
  v::Vector{T}
end
Base.size(V::Hankel) = (div(length(V.v)+1, 2), div(length(V.v)+1, 2))
Oscar.nrows(V::Hankel) = size(V)[1]
Oscar.ncols(V::Hankel) = size(V)[1]
Oscar.base_ring(V::Hankel) = parent(V.v[1])
Base.getindex(V::Hankel, i::Int, j::Int)  = V.v[j+i-1]

function solve_left(V::Vandermonde{T}, b::MatElem{T}) where {T <: Oscar.FieldElem}
  Qx = Oscar.Hecke.Globals.Qx
  x = gen(Qx)
  f = reduce(*, [x-i for i=V.val], init = one(Qx)) #stupid: this is the poly from BM
  z = zero(base_ring(V))
  H = Hankel(vcat([coeff(f, i) for i=1:nrows(V)], [z for i=nrows(V):2*nrows(V)-1]))
  Vt = V'
  #TODO: for fun, use mult-point evaluation to find D
  # and use a special type of matrix again
  D = diagonal_matrix([inv(derivative(f)(x)) for x = V.val])
  return ((b*H)*Vt)*D
end


end

#= from Dereje:

 Qt, t = PolynomialRing(QQ, :t=>1:2)
 Qtx, x = PolynomialRing(FractionField(Qt), :x => 1:3)
 f1 = x[1]^2*x[2]^3*x[3] + 2*t[1]*x[1]*x[2]*x[3]^2 + 7*x[2]^3
 f2 = x[1]^2*x[2]^4*x[3] + (t[1]-7*t[2])*x[1]^2*x[2]*x[3]^2 - x[1]*x[2]^2*x[3]^2+2*x[1]^2*x[2]*x[3] - 12*x[1] + t[2]*x[2]
 f3 = (t[1]^2+t[2]-2)*x[2]^5*x[3] + (t[1]+5*t[2])*x[1]^2*x[2]^2*x[3] - t[2]*x[1]*x[2]^3*x[3] - t[2]*x[1]*x[2]^3*x[3] - x[1]*x[2]^3+x[2]^4+2*t[1]^2*x[2]^2*x[3]
 f4 = t[1]*x[2]^2*x[2]^2*x[3] - x[1]*x[2]^3 *x[3] + (-t[1]+4)*x[2]^3*x[3]^2 + 3*t[1]*x[1]*x[2]*x[3]^3 + 4*x[3]^2 - t[2]*x[1]

 z = MPolyInterpolate.Groebner_basis(ideal(Qtx, [f1, f2, f3, f4]))


=#



