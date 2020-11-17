module ModStdQ

using Oscar
import Hecke
import Oscar: MPolyIdeal, BiPolyArray, Hecke, AbstractAlgebra
import Hecke.MPolyGcd: RecoCtx, rational_reconstruct


function (R::AbstractAlgebra.Generic.MPolyRing{Nemo.gfp_elem})(f::fmpq_mpoly)
  g  = MPolyBuildCtx(R)
  S = base_ring(R)
  for (c, v) in zip(coeffs(f), exponent_vectors(f))
    push_term!(g, S(c), v)
  end
  return finish(g)
end

function (S::Nemo.NmodRing)(a::fmpq)
  return S(numerator(a))//S(denominator(a))
end

#TODO (to dream)
#  the groeber bases in Singular in parallel
#  the rat-reco for different polys in parallel
#  crt in parallel
#  use walk, tracing, ...
function Oscar.groebner_assure(I::MPolyIdeal{fmpq_mpoly}, ord::Symbol = :degrevlex; use_hilbert::Bool = false)
  if isdefined(I, :gb) && ord == :degrevlex
    return I.gb
  end
  ps = Hecke.PrimesSet(Hecke.p_start, -1)
  ps = Hecke.PrimesSet(2^28+2^20, -1)

  p = iterate(ps)[1]
  Qt = base_ring(I)
  Q = base_ring(Qt)
  Zt = PolynomialRing(ZZ, [string(s) for s = symbols(Qt)], cached = false)[1]
  max_stable = 2
  stable = max_stable

  local gc::Array{fmpz_mpoly, 1}
  local gd::Array{fmpq_mpoly, 1}

  fl = true
  d = fmpz(1)
  very_first = true

  while true
    p = iterate(ps, p)[1]
    @vprint :MPolyGcd 2 "Main loop: using $p\n"
#    nbits(d) > 1700 && error("too long")
    R = ResidueRing(ZZ, Int(p))
    Rt, t = PolynomialRing(R, [string(s) for s = symbols(Qt)], cached = false)
    @vtime :MPolyGcd 3 Ip = Oscar.BiPolyArray([Rt(x) for x = gens(I)], keep_ordering = false)
    Jp = map(x->lift(Zt, x), (Oscar.groebner_basis(Ip, ord = ord, complete_reduction = true)))
    if d == 1
      d = fmpz(p)
      gc = Jp
      fl = true
      gd = []
      @vtime :MPolyGcd 2 for f = gc
        fl, fQ = rational_reconstruction(f, d, false, parent = Qt)
        fl || break
        push!(gd, fQ)
      end
      for i = length(gd)+1:length(gc)
        push!(gd, Qt(0))
      end
    else
      @assert length(Jp) == length(gc)
      new_idx = [any(x -> !iszero(R(x)), coefficients(map_coeffs(QQ, Jp[i], parent = Qt) - gd[i])) for i=1:length(gc)]
      @vprint :MPolyGcd 1 "new information in $new_idx\n"
      fl = !any(new_idx)
      if !fl
        @vtime :MPolyGcd 2 for i = 1:length(gc)
          if new_idx[i]
            gc[i], _ = induce_crt(gc[i], d, Jp[i], fmpz(p), true)
          end
        end
        d *= fmpz(p)
        fl = true
        @vtime :MPolyGcd 2 for i = 1:length(gc)
          if new_idx[i]
            fl, gd[i] = rational_reconstruction(gc[i], d, false, parent = Qt)
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
            I.gb = BiPolyArray(gd)
          end
          return gd
        end
      end
    end
    @vprint :MPolyGcd 1 "Information now at $(nbits(d)) bits\n"
  end
end

function Oscar.lift(R::Nemo.Ring, f::nmod_mpoly)
  g = MPolyBuildCtx(R)
  for (c, v) in zip(coeffs(f), exponent_vectors(f))
    push_term!(g, lift(c), v)
  end
  return finish(g)
end

function Hecke.rational_reconstruction(f::fmpz_mpoly, d::fmpz, b::Bool; parent=1)
  g = MPolyBuildCtx(parent)
  for (c, v) in zip(coeffs(f), exponent_vectors(f))
    fl, r, s = Hecke.rational_reconstruction(c, d)
    if !fl
      return false, finish(g)
    end
    push_term!(g, r//s, v)
  end
  return true, finish(g)
end

function Hecke.induce_crt(f::fmpz_mpoly, d::fmpz, g::fmpz_mpoly, p::fmpz, b::Bool)
  mu = MPolyBuildCtx(parent(f))
  for i=1:length(f)
    e = exponent_vector(f, i)
    @assert e == exponent_vector(g, i)
    push_term!(mu, crt(coeff(f, i), d, coeff(g, i), p, b), e)
  end
  return finish(mu), d*p
end

function Oscar.groebner_basis(I::MPolyIdeal{fmpq_mpoly}; ord::Symbol = :degrevlex, complete_reduction::Bool = false)
  H = Oscar.groebner_assure(I, ord)
  return H
end


#TODO? directly project down to Singular???
#      operate on BiPolyArrays or arrays rather than ideals?
# definitely: BiPolyArrays as list is what we do
# for induced stuff and majority voting and such think of data structures
#   that allow to match monomials effectively.

end 
