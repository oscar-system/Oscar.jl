module ModStdQt

using Oscar
import Oscar.Nemo
import Oscar.Hecke

function __init__()
  Hecke.add_verbose_scope(:ModStdQt)
end

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

Oscar.normalise(f::fmpq_poly, ::Int64) = error("no normalise") #length(f)
#Oscar.set_length!(f::fmpq_poly, ::Int64) = error("no set_length") #f

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
    return new{Int}(Int[], p, Int[rand(-10000:10000) for i=1:n], 0)
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
  @vprint :ModStdQt 2 "setting status to $s\n"
  @show s
  M.status = s
end

mutable struct Vals{T}
  v::Array{Array{T, 1}, 1}
  nd::Array{Tuple{<:PolyElem{T}, <:PolyElem{T}}, 1}
  G::Generic.Frac{<:MPolyElem{T}}
  function Vals(v::Array{Array{S, 1}, 1}) where {S}
    r = new{S}()
    r.v = v
    return r
  end
end

#Ben-Or Tiwari...
function Oscar.interpolate(Val::Vals{T}, M::MPolyInterpolateCtx) where {T}

  if isdefined(Val, :G)
    set_status!(M, :OK)
    return true, Val.G
  end

  if !isdefined(Val, :nd) || length(Val.nd) < length(Val.v)
    nd = []
    val = Val.v
    i = 1
    for x = val
      f = interpolate(x, M.I)
      i += 1
      mu = rational_reconstruction(f, M.I, ErrorTolerant = true)
      if !mu[1]
        set_status!(M, :univariate_failed)
        return false, zero(M.R) # more z
      end
      t = inv(trail(mu[3]))
      push!(nd, (mu[2]*t, mu[3]*t))
    end
    Val.nd = nd
  else
    nd = Val.nd
  end
  Qx = parent(nd[1][1])
  @assert parent(nd[1][2]) == Qx

  R = [M.R(), M.R()]
  for ii = 1:2
    cor = elem_type(Qx)[zero(Qx) for i=1:length(nd)]
    for i=maximum(degree(x[ii]) for x = nd):-1:0
      bm = berlekamp_massey([coeff(nd[x][ii]-cor[x], i) for x = 1:length(nd)], ErrorTolerant = true)
      if !bm[1]
        set_status!(M, :BM_failed)
        return false, zero(M.R)
      end
      if isconstant(bm[2])
        continue #why? but is neccessary
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
      s = solve_left(V, matrix(QQ, 1, length(r), [coeff(nd[x][ii]-cor[x], i) for x=1:length(r)]))
      G = MPolyBuildCtx(M.R)
      for i=1:length(r)
        push_term!(G, s[1,i], [valuation(r[i], p) for p = M.pt.pt_p])
      end
      g = finish(G)
      R[ii] += g
      for j=0:length(nd)-1
        t = evaluate(g, [T(M.pt.pt_p[l])^j*gen(Qx) + M.pt.pt_s[l] for l=1:length(M.pt.pt_p)])
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
  Val.G = R[1]//R[2]

  return true, Val.G
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

function Oscar.groebner_basis(I::Oscar.MPolyIdeal{<:Generic.MPoly{<:Generic.Frac{fmpq_mpoly}}}; ord::Symbol = :degrevlex, complete_reduction::Bool = true)
  Oscar.groebner_assure(I, ord, complete_reduction = complete_reduction)
end

function Oscar.groebner_assure(I::Oscar.MPolyIdeal{<:Generic.MPoly{<:Generic.Frac{fmpq_mpoly}}}, ord::Symbol = :degrevlex; complete_reduction::Bool = true)
  T = fmpq
  if ord == :degrevlex && isdefined(I, :gb)
    return I.gb
  end
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
      @vprint :ModStdQt 1 "adding z-evaluation $(length(P.pt_z)+1)\n"
      push!(P, length(P.pt_z)+1)
    else
      @vprint :ModStdQt 1 "increasing power-of-prime to $(P.j+1)\n"
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
      @vprint :ModStdQt 2 "using evaluation point $pt\n"
      @vtime :ModStdQt 2 for g = gens(I)
        G = MPolyBuildCtx(Qy)
        for (c, e) = zip(Generic.MPolyCoeffs(g), Generic.MPolyExponentVectors(g))
          push_term!(G, evaluate(c, pt, ErrorTolerant = true), e)
        end
        push!(gens_J, finish(G))
      end
      J = ideal(Qy, gens_J)
      @vtime :ModStdQt 2 gJ = groebner_basis(J, ord = ord, complete_reduction = true, Proof = false)
      if frst
        lst = []
        for _g = gJ
          g = inv(lead(_g))*_g
          f = []
          for (c, e) = zip(Generic.MPolyCoeffs(g), Generic.MPolyExponentVectors(g))
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
          for (c, e) = zip(Generic.MPolyCoeffs(g), Generic.MPolyExponentVectors(g))
            if lst[ig][jg][1] != e #TODO: sort and match
              @assert lst[ig][jg][1] == e
              @show ig, jg
              for y = 1:length(lst)
                @show y, [x[1] for x = lst[y]]
              end
              @show gJ[ig]
              @show gJ
              if do_z
                lst[ig][jg][2][idx-1, length(P.pt_z)] = T(0)
              else
                lst[ig][jg][2][P.j+1, idx-1] = T(0)
              end
              jg += 1
              @assert lst[ig][jg][1] == e
            end
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
    
    res = elem_type(S)[]
    @vprint :ModStdQt 1 "starting interpolation...\n"
    M = MPolyInterpolateCtx(R, P)
    local fl = true
    @vtime :ModStdQt 2 for F = lst
      f = MPolyBuildCtx(S)
      for C = F
        fl, c = interpolate(C[2], M)
        if fl
          push_term!(f, c, C[1])
        else
          break
        end
      end
      fl || break
      push!(res, finish(f))
    end
    if fl
      @vprint :ModStdQt 1 "success\n"
      if ord == :degrevlex
        I.gb = Oscar.BiPolyArray(res)
      end
      return res
    end

    @vprint :ModStdQt 1 "failed, more points\n"
    
    if M.status == :univariate_failed
      @show "needs more z"
      do_z = true
      continue
    end

    if M.status == :BM_failed && P.j > 10 
      @show "needs more j"
      do_z = !true
      continue
    end


    do_z = !do_z
#    if length(P.pt_z) > 15 || P.j > 15
#      error("did not finish")
#      break
#    end
  end
  return lst, P
end


#=
 From
 Sparse Polynomial Interpolation and
  Berlekamp/Massey Algorithms That Correct Outlier
  Errors in Input Values

  https://dl.acm.org/doi/10.1145/2442829.2442852

=#

function cleanup(Lambda::PolyElem{T}, E::Int, c::Array{T, 1}, k::Int) where {T}
  c = copy(c)
  t = degree(Lambda)
  i = k+2*t
  e = 0
  n = length(c)
  while i <= n-1
    d = sum(coeff(Lambda, j)*c[i+j-t+1] for j=0:t-1)
    if d != -c[i+1] 
      c[i+1] = -d
      e += 1
    end
    i += 1
  end
  i = k-1
  while i>= 0 && e <= E
    d = sum(coeff(Lambda, j)*c[i+j+1] for j=1:t)
    if d+coeff(Lambda, 0)*c[i+1] != 0
      c[i+1] = -d//coeff(Lambda, 0)
      e += 1
    end
    i -= 1
  end
  return c, e
end

function FTBM(c::Array{Q, 1}, T::Int, E::Int) where {Q <: FieldElem}
  R = parent(c[1])
  Rt, t = PolynomialRing(R, cached = false)
  L = []
  n = length(c)
  for i=0:div(n, 2*T)-1
    Lambda = berlekamp_massey(c[2*T*i+1:2*T*(i+1)], parent = Rt)[2]
    if constant_coefficient(Lambda) != 0
      a, e = cleanup(Lambda, E, c, 2*T*i)
      if e <= E
        push!(L, (a, e))
      end
    end
  end
  return L
end

function MRBM(c::Array{Q, 1}, T::Int, E::Int) where {Q <: FieldElem}
  R = parent(c[1])
  Rt, t = PolynomialRing(R, cached = false)
  L = Dict{elem_type(Rt), Set{Int}}()
  n = length(c)
  for i=0:div(n, 2*T)-1
    Lambda = berlekamp_massey(c[2*T*i+1:2*T*(i+1)], parent = Rt)[2]
    if constant_coefficient(Lambda) == 0 || 2*degree(Lambda)+1 >length(c)
      continue
    end
    @assert parent(Lambda) == Rt
    if haskey(L, Lambda)
      push!(L[Lambda], i)
    else
      L[Lambda] = Set(i)
    end
  end
  if length(L) == 0
    return false, t
  end
  m = 0
  local Lambda
  for (k, v) = L
    if length(v) > m
      m = length(v)
      Lambda = k
    end
  end
  for i = L[Lambda]
    a, e = cleanup(Lambda, E, c, 2*T*i)
    if e <= E
      return true, Lambda
    end
  end
  return false, t
end
function MRBM(c::Array{Q, 1}) where {Q <: FieldElem}
  n = length(c)
  #we need, by paper: (2T)(2E+1) <= n
  #assuming E < T/2 => (2T)(T+1) = 2T^2+2T <= n
  #so T <= root(n)/2
  T = div(root(n, 2), 2)
  E = floor(Int, (n/2/T-1)/2)
  local L
  for e = 0:E
    T = floor(Int, n/(2*e+1)/2)
    fl, L = MRBM(c, T, e)
    if fl
      return fl, L
    end
  end
  return false, L
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

 z = groebner_basis(ideal(Qtx, [f1, f2, f3, f4]))
 z = groebner_basis(ideal(Qtx, [f1, f2, f3, f4]), ord = :lex)


=#



