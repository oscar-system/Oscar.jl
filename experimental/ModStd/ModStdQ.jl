module ModStdQ

using Oscar
import Hecke
import Oscar: MPolyIdeal, BiPolyArray, IdealGens, Hecke, AbstractAlgebra
import Hecke: induce_rational_reconstruction, induce_crt

function __init__()
  Hecke.add_verbose_scope(:ModStdQ)
end

function (R::GFPMPolyRing)(f::fmpq_mpoly)
  g  = MPolyBuildCtx(R)
  S = base_ring(R)
  for (c, v) in zip(coefficients(f), exponent_vectors(f))
    push_term!(g, S(c), v)
  end
  return finish(g)
end

function (S::Union{Nemo.NmodRing, Nemo.GaloisField})(a::fmpq)
  return S(numerator(a))//S(denominator(a))
end

#TODO (to dream)
#  the groeber bases in Singular in parallel
#  the rat-reco for different polys in parallel
#  crt in parallel
#  use walk, tracing, ...

#= TODO: Currently we had to "disable" modular GB stuff due to introducing dictionaries of GBs for ideals.
 =     Next step is to reenable modular Singular.std and modular f4 again. =#
function exp_groebner_assure(I::MPolyIdeal{fmpq_mpoly}, ord::Symbol = :degrevlex; use_hilbert::Bool = false, Proof::Bool = true)
  if isdefined(I, :gb) && ord == :degrevlex
    return collect(I.gb)
  end
  if Proof
    return Oscar.groebner_basis_with_transform(I, ordering = ord)[1]
  end

  ps = Hecke.PrimesSet(Hecke.p_start, -1)
  ps = Hecke.PrimesSet(2^28+2^20, -1)

  p = iterate(ps)[1]
  Qt = base_ring(I)
  Q = base_ring(Qt)
  Zt = PolynomialRing(ZZ, [string(s) for s = symbols(Qt)], cached = false)[1]
  max_stable = 2
  stable = max_stable

  local gc::Vector{fmpz_mpoly}
  local gd::Vector{fmpq_mpoly}

  fl = true
  d = fmpz(1)
  very_first = true
  gI = gens(I)

  while true
    p = iterate(ps, p)[1]
    @vprint :ModStdQ 2 "Main loop: using $p\n"
#    nbits(d) > 1700 && error("too long")
    R = ResidueRing(ZZ, Int(p)) #gfp_mpoly missing...
    Rt, t = PolynomialRing(R, [string(s) for s = symbols(Qt)], cached = false)
    @vtime :ModStdQ 3 Ip = Oscar.IdealGens([Rt(x) for x = gI], keep_ordering = false)
    Gp = Oscar.exp_groebner_basis(Ip, ord = ord, complete_reduction = true)
    Jp = map(x->lift(Zt, x), Gp)
    if d == 1
      d = fmpz(p)
      gc = Jp
      fl = true
      gd = []
      @vtime :ModStdQ 2 for f = gc
        fl, fQ = induce_rational_reconstruction(f, d, false, parent = Qt)
        fl || break
        push!(gd, fQ)
      end
      for i = length(gd)+1:length(gc)
        push!(gd, Qt(0))
      end
    else
      @assert length(Jp) == length(gc)
      new_idx = [any(x -> !iszero(R(x)), coefficients(map_coefficients(QQ, Jp[i], parent = Qt) - gd[i])) for i=1:length(gc)]
      @vprint :ModStdQ 1 "new information in $new_idx\n"
      fl = !any(new_idx)
      if !fl
        @vtime :ModStdQ 2 for i = 1:length(gc)
          if new_idx[i]
            gc[i], _ = induce_crt(gc[i], d, Jp[i], fmpz(p), true)
          end
        end
        d *= fmpz(p)
        fl = true
        @vtime :ModStdQ 2 for i = 1:length(gc)
          if new_idx[i]
            fl, gd[i] = induce_rational_reconstruction(gc[i], d, false, parent = Qt)
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
            I.gb[degrevlex(gens(Qt))] = IdealGens(gd, keep_ordering = false, isGB = true)
          end
          return gd
        end
      end
    end
    @vprint :ModStdQ 1 "Information now at $(nbits(d)) bits\n"
  end
end

function groebner_basis_with_transform_inner(I::MPolyIdeal{fmpq_mpoly}, ord::MonomialOrdering; complete_reduction::Bool = true, use_hilbert::Bool = false)
  if iszero(I)
    I.gb[ord] = IdealGens(base_ring(I), fmpq_mpoly[], isGB = true, keep_ordering = false)
    singular_assure(I.gb[ord])
    return fmpq_mpoly[], matrix(base_ring(I), ngens(I), 0, fmpq_mpoly[])
  end
    
  ps = Hecke.PrimesSet(Hecke.p_start, -1)
  ps = Hecke.PrimesSet(2^28+2^20, -1)

  p = iterate(ps)[1]
  Qt = base_ring(I)
  Q = base_ring(Qt)
  Zt = PolynomialRing(ZZ, [string(s) for s = symbols(Qt)], cached = false)[1]
  max_stable = 2
  stable = max_stable

  local gc::Vector{fmpz_mpoly}
  local gd::Vector{fmpq_mpoly}
  local length_gc::Int

  fl = true
  d = fmpz(1)
  very_first = true

  gI = gens(I)

  while true
    p = iterate(ps, p)[1]
    @vprint :ModStdQ 2 "Main loop: using $p\n"
#    nbits(d) > 1700 && error("too long")
    R = GF(p)
    Rt, t = PolynomialRing(R, [string(s) for s = symbols(Qt)], cached = false)
    @vtime :ModStdQ 3 Ip = Oscar.IdealGens([Rt(x) for x = gI], keep_ordering = false)
    Gp, Tp = Oscar.groebner_basis_with_transform(Ip, ord, complete_reduction)
    length_gc = length(Gp)
    Jp = vcat(map(x->lift(Zt, x), Gp), map(x->lift(Zt, x), reshape(collect(Tp), :)))

    if d == 1
      d = fmpz(p)
      gc = Jp
      fl = true
      gd = []
      @vtime :ModStdQ 2 for f = gc
        fl, fQ = Hecke.induce_rational_reconstruction(f, d, false, parent = Qt)
        fl || break
        push!(gd, fQ)
      end
      for i = length(gd)+1:length(gc)
        push!(gd, Qt(0))
      end
    else
      @assert length(Jp) == length(gc)
      new_idx = [any(x -> !iszero(R(x)), coefficients(map_coefficients(QQ, Jp[i], parent = Qt) - gd[i])) for i=1:length(gc)]
      @vprint :ModStdQ 1 "new information in $new_idx\n"
      fl = !any(new_idx)
      if !fl
        @vtime :ModStdQ 2 for i = 1:length(gc)
          if new_idx[i]
            gc[i], _ = induce_crt(gc[i], d, Jp[i], fmpz(p), true)
          end
        end
        d *= fmpz(p)
        fl = true
        @vtime :ModStdQ 2 for i = 1:length(gc)
          if new_idx[i]
            fl, gd[i] = Hecke.induce_rational_reconstruction(gc[i], d, false, parent = Qt)
            fl || break
          end
        end
        stable = max_stable
#        @show gd
      else
        d *= fmpz(p)
        stable -= 1
        if stable <= 0
          G = gd[1:length_gc]
          T = matrix(Qt, length_gc, length(gI), gd[length_gc+1:end])
          #at this point we SHOULD have T*gens(I) == G...
          if T*matrix(Qt, length(gI), 1, gI) == matrix(Qt, length_gc, 1, G)
            if !isdefined(I.gens, :ord)
               I.gens.ord = ord.o
            end
            if ord.o == I.gens.ord && !isdefined(I, :gb)
              I.gb[ord] = IdealGens(gd[1:length_gc], keep_ordering = false, isGB = true)
              singular_assure(I.gb[ord])
            end
            return G, T
          else
            @vprint :ModStdQ 1 "all lifts are stable, result is wrong (via trafo)"
          end
        end
      end
    end
    @vprint :ModStdQ 1 "Information now at $(nbits(d)) bits\n"
  end
end

#= function Oscar.groebner_basis_with_transform(I::MPolyIdeal{fmpq_mpoly}; ordering::Symbol = :degrevlex, complete_reduction::Bool = true, use_hilbert::Bool = false)
 =    ord = Oscar.Orderings.MonomialOrdering(base_ring(I), Oscar.Orderings.ordering(gens(base_ring(I)), ordering))
 =    return groebner_basis_with_transform_inner(I, ord; complete_reduction=complete_reduction, use_hilbert=use_hilbert)
 = end
 =  =#
 function Oscar.groebner_basis_with_transform(I::MPolyIdeal{fmpq_mpoly}, ord::MonomialOrdering=default_ordering(base_ring(I)); complete_reduction::Bool = true, use_hilbert::Bool = false)
   return groebner_basis_with_transform_inner(I, ord; complete_reduction=complete_reduction, use_hilbert=use_hilbert)
end

function Oscar.lift(R::Nemo.Ring, f::Union{gfp_mpoly, nmod_mpoly})
  g = MPolyBuildCtx(R)
  for (c, v) in zip(coefficients(f), exponent_vectors(f))
    push_term!(g, lift(c), v)
  end
  return finish(g)
end

function induce_rational_reconstruction(f::fmpz_mpoly, d::fmpz, b::Bool; parent=1)
  g = MPolyBuildCtx(parent)
  for (c, v) in zip(coefficients(f), exponent_vectors(f))
    fl, r, s = Hecke.rational_reconstruction(c, d)
    if !fl
      return false, finish(g)
    end
    push_term!(g, r//s, v)
  end
  return true, finish(g)
end

function induce_crt(f::fmpz_mpoly, d::fmpz, g::fmpz_mpoly, p::fmpz, b::Bool)
  mu = MPolyBuildCtx(parent(f))
  for i=1:length(f)
    e = exponent_vector(f, i)
    @assert e == exponent_vector(g, i)
    push_term!(mu, crt(coeff(f, i), d, coeff(g, i), p, b), e)
  end
  return finish(mu), d*p
end

function exp_groebner_basis(I::MPolyIdeal{fmpq_mpoly}; ord::Symbol = :degrevlex, complete_reduction::Bool = false)
  H = exp_groebner_assure(I, ord)
  return H
end


#TODO? directly project down to Singular???
#      operate on IdealGens or arrays rather than ideals?
# definitely: IdealGens as list is what we do
# for induced stuff and majority voting and such think of data structures
#   that allow to match monomials effectively.

end 
