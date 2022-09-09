module ModStdNF

using Oscar
import Hecke
import Oscar: MPolyIdeal, IdealGens, Hecke, AbstractAlgebra
import Hecke: modular_lift, modular_proj, modular_env, RecoCtx, 
              induce_rational_reconstruction

function __init__()
  Hecke.add_verbose_scope(:ModStdNF)
end


function stdhilb(I::Singular.sideal, h::Vector{Int32}; complete_reduction::Bool=false)
  R = base_ring(I)
#  @show "stdhilb", I
  @time ptr = Singular.libSingular.id_StdHilb(I.ptr, R.ptr, h, complete_reduction)
  Singular.libSingular.idSkipZeroes(ptr)
  z = Singular.Ideal(R, ptr)
  z.isGB = true
  return z
end

function Oscar.binomial(a::RingElem, k::Int)
  p = parent(a)
  return prod([a-i for i=0:k-1])*inv(p(factorial(k)))
end

function exp_groebner_basis(B::IdealGens{nmod_mpoly}, h::HilbertData; ord::Symbol = :degrevlex, complete_reduction::Bool = false)
  if ord != :degrevlex
    R = Oscar.singular_poly_ring(B.Ox, ord)
    i = stdhilb(Singular.Ideal(R, [convert(R, x) for x = B]), h.data, complete_reduction = complete_reduction)
    return IdealGens(B.Ox, i)
  end
  if !isdefined(B, :S)
    B.S = Singular.Ideal(B.Sx, [convert(B.Sx, x) for x = B.O])
  end 
  return IdealGens(B.Ox, stdhilb(B.S, h.data, complete_reduction = complete_reduction), keep_ordering = false, isGB = true)
end

#TODO (to dream)
#  the groeber bases in Singular in parallel
#  the rat-reco for different polys in parallel
#  crt in parallel
#  use walk, tracing, ...
function exp_groebner_assure(I::MPolyIdeal{Generic.MPoly{nf_elem}}, ord::Symbol = :degrevlex; use_hilbert::Bool = false)
  if isdefined(I, :gb) && ord == :degrevlex
    return I.gb
  end
  ps = Hecke.PrimesSet(Hecke.p_start, -1)
  ps = Hecke.PrimesSet(2^28+2^20, -1)

  p = iterate(ps)[1]
  Kt = base_ring(I)
  K = base_ring(Kt)
  max_stable = 2
  stable = max_stable

  local gc::Vector{Generic.MPoly{nf_elem}}
  local gd::Vector{Generic.MPoly{nf_elem}}
  Zx = Hecke.Globals.Zx
  R = RecoCtx(K)

  fl = true
  d = fmpz(1)
  very_first = true
  local H::HilbertData
  while true
    p = iterate(ps, p)[1]
    @vprint :ModStdNF 2 "Main loop: using $p\n"
    @vtime :ModStdNF 3 me = Hecke.modular_init(K, p, deg_limit = 1)
#    nbits(d) > 1700 && error("too long")
    if isempty(me)
      continue
    end

    @vtime :ModStdNF 3 Ip = Hecke.modular_proj(I.gens, me)
    Jp = typeof(Ip[1])[]
    @vtime :ModStdNF 2 for fp = Ip
      if use_hilbert
        if very_first
          @show H = HilbertData(fp)
          very_first = false
        end
        push!(Jp, exp_groebner_basis(fp, H, ord = ord, complete_reduction = true))
      else
        push!(Jp, exp_groebner_basis(fp, ord = ord, complete_reduction = true))
      end
    end
    @vtime :ModStdNF 2 IP = Hecke.modular_lift(Jp, me)
    if d == 1
      d = fmpz(p)
      gc = IP
      push!(R, d, lift(Zx, me.ce.pr[end]))
      fl = true
      gd = []
      @vtime :ModStdNF 2 for f = gc
        fl, fQ = Hecke.induce_rational_reconstruction(f, R, integral = false)
        fl || break
        push!(gd, fQ)
      end
      for i = length(gd)+1:length(gc)
        push!(gd, Kt(0))
      end
#      @show gd
    else
      new_idx = [any(x -> any(x->!iszero(x), Hecke.modular_proj(x, me)), coefficients(gd[i] - IP[i])) for i=1:length(gc)]
      @vprint :ModStdNF 1 "new information in $new_idx\n"
      push!(R, fmpz(p), lift(Zx, me.ce.pr[end]))
      fl = !any(new_idx)
      if !fl
        @vtime :ModStdNF 2 for i = 1:length(gc)
          if new_idx[i]
            gc[i], _ = induce_crt(gc[i], d, IP[i], fmpz(p), true)
          end
        end
        d *= fmpz(p)
        fl = true
        @vtime :ModStdNF 2 for i = 1:length(gc)
          if new_idx[i]
            fl, gd[i] = Hecke.induce_rational_reconstruction(gc[i], R, integral = false)
            fl || break
          end
        end
        stable = max_stable
#        @show gd
      else
        d *= fmpz(p)
        stable -= 1
        if stable <= 0
          if ord == :degrevlex
            I.gb = IdealGens(gd, keep_ordering = false, isGB = true)
          end
          return gd
        end
      end
    end
    @vprint :ModStdNF 1 "Information now at $(nbits(d)) bits\n"
  end
end

function exp_groebner_basis(I::MPolyIdeal{Generic.MPoly{nf_elem}}; ord::Symbol = :degrevlex, complete_reduction::Bool = true)
  return Oscar.exp_groebner_assure(I, ord)
end


#TODO? directly project down to Singular???
#      operate on IdealGenss or arrays rather than ideals?
# definitely: IdealGenss as list is what we do
# for induced stuff and majority voting and such think of data structures
#   that allow to match monomials effectively.
function Hecke.modular_proj(B::IdealGens{Generic.MPoly{nf_elem}}, me::Hecke.modular_env)
  g = [Vector{nmod_mpoly}() for i = me.fld]
  for i=B
    h = Hecke.modular_proj(i, me)
    for j = 1:length(h)
      push!(g[j], h[j])
    end
  end
  return [IdealGens(x, keep_ordering = false) for x = g] 
end

function Hecke.modular_lift(f::Vector{IdealGens{nmod_mpoly}}, me::Hecke.modular_env)
  g = []
  @assert all(x -> length(x) == length(f[1]), f)
  for i=1:length(f[1])
    lp = nmod_mpoly[ f[j][Val(:O), i] for j=1:length(f)]
    push!(g, Hecke.modular_lift(lp, me))
  end
  return g
end

function homogenize(f::MPolyElem, S::MPolyRing)
  d = total_degree(f)
  g = MPolyBuildCtx(S)
  for (c, e) = Base.Iterators.zip(Generic.MPolyCoeffs(f), Generic.MPolyExponentVectors(f))
    push_term!(g, c, push!(e, d-sum(e)))
  end
  return finish(g)
end

function homogenize(i::MPolyIdeal)
  R = base_ring(i)
  S, _ = PolynomialRing(base_ring(R), push!([string(x) for x = symbols(R)], "H"))
  return homogenize(i, S)
end

function homogenize(i::MPolyIdeal, S::MPolyRing)
  return ideal([homogenize(x, S) for x = gens(i)])
end

#Hecke.induce_crt and such from MPolyGcd...

end

module DerejeGB
using Oscar

function example_1()
  k, a = quadratic_field(-1)
  kt, (x,y,z) = PolynomialRing(k, ["x", "y", "z"])
  f1 = (a+8)*x^2*y^2+5*x*y^3+(-a+3)*x^3*z+x^2*y*z
  f2 = x^5+2*y^3*z^2+13*y^2*z^3+5*y*z^4
  f3 = 8*x^3+(a+12)*y^3+x*z^2+3
  f4 = (-a+7)*x^2*y^4+y^3*z^3+18*y^3*z^2

  return ideal([f1, f2, f3, f4])
end

function example_2()
  t = gen(Hecke.Globals.Qx)
  k, a = number_field(t^5+t^2+2)
  kt, (x, y, z) = PolynomialRing(k, ["x", "y", "z"])

  f1 = 2*x*y^4*z^2+(a-1)*x^2*y^3*z+(2*a)*x*y*z^2+7*y^3+(7*a+1)
  f2 = 2*x^2*y^4*z+(a)*x^2*y*z^2-x*y^2*z^2+(2*a+3)*x^2*y*z-12*x+(12*a)*y
  f3 = (2*a)*y^5*z+x^2*y^2*z-x*y^3*z+(-a)*x*y^3+y^4+2*y^2*z
  f4 = (3*a)*x*y^4*z^3+(a+1)*x^2*y^2*z-x*y^3*z+4*y^3*z^2+(3*a)*x*y*z^3+4*z^2-x+(a)*y

  return ideal([f1, f2, f3, f4])
end

function example_3()
  t = gen(Hecke.Globals.Qx)
  k, a = number_field(t^7-7*t+3)

  kt, (v, w, x, y, z) = PolynomialRing(k, ["v", "w", "x", "y", "z"])
  f1 = (a)*v+(a-1)*w+x+(a+2)*y+z
  f2 = v*w+(a-1)*w*x+(a+2)*v*y+x*y+(a)*y*z
  f3 = (a)*v*w*x+(a+5)*w*x*y+(a)*v*w*z+(a+2)*v*y*z+(a)*x*y*z
  f4 = (a-11)*v*w*x*y+(a+5)*v*w*x*z+(a)*v*w*y*z+(a)*v*x*y*z+(a)*w*x*y*z
  f5 = (a+3)*v*w*x*y*z+(a+23)

  return ideal([f1, f2, f3, f4, f5])
end

function example_4()
  t = gen(Hecke.Globals.Qx)
  k, a = number_field(t^7-7*t+3)

  kt, (u, v, w, x, y, z) = PolynomialRing(k, ["u", "v", "w", "x", "y", "z"])

  f1 = (a)*u+(a+2)*v+w+x+y+z
  f2 = u*v+v*w+w*x+x*y+(a+3)*u*z+y*z
  f3 = u*v*w+v*w*x+(a+1)*w*x*y+u*v*z+u*y*z+x*y*z
  f4 = (a-1)*u*v*w*x+v*w*x*y+u*v*w*z+u*v*y*z+u*x*y*z+w*x*y*z
  f5 = u*v*w*x*y+(a+1)*u*v*w*x*z+u*v*w*y*z+u*v*x*y*z+u*w*x*y*z+v*w*x*y*z
  f6 = u*v*w*x*y*z+(-a+2)

  return ideal([f1, f2, f3, f4, f5, f6])
end

function example_5()
  k, a = cyclotomic_field(7)
  kt, (w, x, y, z) = PolynomialRing(k, ["w", "x", "y", "z"])

  f1 = (a+5)*w^3*x^2*y+(a-3)*w^2*x^3*y+(a+7)*w*x^2*y^2
  f2 = (a)*w^5+(a+3)*w*x^2*y^2+(a^2+11)*x^2*y^2*z
  f3 = (a+7)*w^3+12*x^3+4*w*x*y+(a)*z^3
  f4 = 3*w^3+(a-4)*x^3+x*y^2

  return ideal([f1, f2, f3, f4])
end

function example_6()
  t = gen(Hecke.Globals.Qx)
  k, a = number_field(t^12-5*t^11+24*t^10-115*t^9+551*t^8-2640*t^7+12649*t^6-2640*t^5+551*t^4-115*t^3+24*t^2-5*t+1)

  kt, (w, x, y, z) = PolynomialRing(k, ["w", "x", "y", "z"])

  f1 = (2*a+3)*w*x^4*y^2+(a+1)*w^2*x^3*y*z+2*w*x*y^2*z^3+(7*a-1)*x^3*z^4
  f2 = 2*w^2*x^4*y+w^2*x*y^2*z^2+(-a)*w*x^2*y^2*z^2+(a+11)*w^2*x*y*z^3-12*w*z^6+12*x*z^6
  f3 = 2*x^5*y+w^2*x^2*y*z-w*x^3*y*z-w*x^3*z^2+(a)*x^4*z^2+2*x^2*y*z^3
  f4 = 3*w*x^4*y^3+w^2*x^2*y*z^3-w*x^3*y*z^3+(a+4)*x^3*y^2*z^3+3*w*x*y^3*z^3+(4*a)*y^2*z^6-w*z^7+x*z^7

  return ideal([f1, f2, f3, f4])
end

function example_7()
  t = gen(Hecke.Globals.Qx)
  k, a = number_field(t^2+5*t+1)

  kt, (u, v, w, x, y, z) = PolynomialRing(k, ["u", "v", "w", "x", "y", "z"])

  f1 = u+v+w+x+y+z+(a)
  f2 = u*v+v*w+w*x+x*y+y*z+(a)*u+(a)*z
  f3 = u*v*w+v*w*x+w*x*y+x*y*z+(a)*u*v+(a)*u*z+(a)*y*z
  f4 = u*v*w*x+v*w*x*y+w*x*y*z+(a)*u*v*w+(a)*u*v*z+(a)*u*y*z+(a)*x*y*z
  f5 = u*v*w*x*y+v*w*x*y*z+(a)*u*v*w*x+(a)*u*v*w*z+(a)*u*v*y*z+(a)*u*x*y*z+(a)*w*x*y*z
  f6 = u*v*w*x*y*z+(a)*u*v*w*x*y+(a)*u*v*w*x*z+(a)*u*v*w*y*z+(a)*u*v*x*y*z+(a)*u*w*x*y*z+(a)*v*w*x*y*z
  f7 = (a)*u*v*w*x*y*z-1

  return ideal([f1, f2, f3, f4, f5, f6, f7])
end

function example_8()
  t = gen(Hecke.Globals.Qx)
  k, a = number_field(t^8-16*t^7+19*t^6-t^5-5*t^4+13*t^3-9*t^2+13*t+17)

  kt, (w, x, y, z) = PolynomialRing(k, ["w", "x", "y", "z"])

  f1 = (-a^2-1)*x^2*y+2*w*x*z-2*w+(a^2+1)*y
  f2 = (a^3-a-3)*w^3*y+4*w*x^2*y+4*w^2*x*z+2*x^3*z+(a)*w^2-10*x^2+4*w*y-10*x*z+(2*a^2+a)
  f3 = (a^2+a+11)*x*y*z+w*z^2-w-2*y
  f4 = -w*y^3+4*x*y^2*z+4*w*y*z^2+2*x*z^3+(2*a^3+a^2)*w*y+4*y^2-10*x*z-10*z^2+(3*a^2+5);

  return ideal([f1, f2, f3, f4])
end

function example_9()
  t = gen(Hecke.Globals.Qx)
  k, a = number_field(t^7+10*t^5+5*t^3+10*t+1)

  kt, (t, u, v, w, x, y, z) = PolynomialRing(k, ["t", "u", "v", "w", "x", "y", "z"])

  f1 = v*x+w*y-x*z-w-y
  f2 = v*w-u*x+x*y-w*z+v+x+z
  f3 = t*w-w^2+x^2-t
  f4 = (-a)*v^2-u*y+y^2-v*z-z^2+u
  f5 = t*v+v*w+(-a^2-a-5)*x*y-t*z+w*z+v+x+z+(a+1)
  f6 = t*u+u*w+(-a-11)*v*x-t*y+w*y-x*z-t-u+w+y
  f7 = w^2*y^3-w*x*y^3+x^2*y^3+w^2*y^2*z-w*x*y^2*z+x^2*y^2*z+w^2*y*z^2-w*x*y*z^2+x^2*y*z^2+w^2*z^3-w*x*z^3+x^2*z^3
  f8 = t^2*u^3+t^2*u^2*v+t^2*u*v^2+t^2*v^3-t*u^3*x-t*u^2*v*x-t*u*v^2*x-t*v^3*x+u^3*x^2+u^2*v*x^2+u*v^2*x^2+v^3*x^2;

  return ideal([f1, f2, f3, f4, f5, f6, f7, f8])
end

end #DerejeGB

module Decker

using Oscar

function example_dim(k, n::Int, d::Int, nc::Int, range)
  kt, t = PolynomialRing(k, 2*n+1)

  id = typeof(t[1])[]
  for i=1:n
    z = kt(0)
    for j=1:nc
      m = kt(1)
      for k=1:d
        m *= rand(t)
      end
      z += rand(k, range)*m
    end
    push!(id, z)
  end

  # the benchmark is for the dimension of the ideal
  return ideal(id)
end

function book_page_80(k)
  Qx, (x, y, z) = PolynomialRing(k, ["x", "y", "z"])
  i = ideal([3x^3*y+x^3+x*y^3+y^2*z^2,2x^3*z-x*y-x*z^3-y^4-z^2,2x^2*y*z-2x*y^2+x*z^2-y^4])
  return i
end

end #Decker


module Fieker
using Oscar

function example_bad_lex()
  k, i = quadratic_field(-1)
  kx, (x1, x2, x3, x4, x5, x6, x7) = PolynomialRing(k, 7)
  i = ideal(
      [(5*i + 7)*x1^2*x5 + (4*i + 9)*x1*x2*x4 + (10*i + 2)*x2*x3*x6 + (10*i + 9)*x2*x4*x6,
       (6*i + 8)*x1*x2*x3 + (7*i + 3)*x2^2*x5 + (i + 3)*x2*x4*x6 + (3*i + 10)*x2*x6^2,
       (i + 7)*x1*x4*x7 + (8*i + 10)*x2*x6*x7 + (6*i + 4)*x3*x4*x6 + (8*i + 7)*x3*x6^2])

  return i
end

end #Fieker
