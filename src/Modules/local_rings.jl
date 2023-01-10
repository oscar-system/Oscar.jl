export length
export composition_series

function length(M::ModuleFP{RingElemType}) where {RingElemType<:AbsLocalizedRingElem{<:Any, <:Any, <:MPolyComplementOfPrimeIdeal}}
  return length(composition_series(M))
end

@attr function composition_series(
    M::ModuleFP{RingElemType}
  ) where {RingElemType<:AbsLocalizedRingElem{<:Any, <:Any, <:MPolyComplementOfPrimeIdeal}}
  if iszero(M) 
    return Vector{elem_type(M)}()
  end
  R = base_ring(M)
  P = prime_ideal(inverted_set(R))
  F = FreeMod(R, 1)
  N, _ = quo(F, (R(P)*F)[1])
  k = 0
  MM = M
  p = identity_map(M)
  NtoM, interp = hom(N, MM)
  comp = Vector{elem_type(M)}()
  while !iszero(NtoM)
    g = gens(NtoM)
    i = findfirst(x->!(iszero(x)), g)
    Msub, inc = sub(MM, [interp(g[i])(N[1])])
    push!(comp, preimage(p, inc(Msub[1])))
    MM, pp = quo(MM, Msub)
    p = compose(p, pp)
    NtoM, interp = hom(N, MM)
    k = k + 1
  end
  iszero(MM) || error("module does not have finite length")
  return comp
end

# TODO: 
# For localizations at 𝕜-point ideals 𝔪 = ⟨x₁-a₁, …, xₙ-aₙ⟩ over polynomial 
# rings we should have a special dispatch for the above method. 
# This should then translate the point to the origin and proceed via 
# standard-basis computations, cf. _vdim_hack(I::MPolyIdeal).

function composition_series(
        M::ModuleFP{T}
    ) where {
             T<:MPolyLocalizedRingElem{<:Field, <:FieldElem, 
                                       <:MPolyRing, <:MPolyElem,
                                       <:MPolyComplementOfKPointIdeal
                                      }
            }
    error("composition_series not implemented; see the source code for details")
end

function composition_series(
        M::ModuleFP{T}
    ) where {
             T<:MPolyQuoLocalizedRingElem{<:Field, <:FieldElem, 
                                          <:MPolyRing, <:MPolyElem,
                                          <:MPolyComplementOfKPointIdeal
                                         }
            }
    error("composition_series not implemented; see the source code for details")
end

function length(
        M::ModuleFP{RingElemType}
    ) where {
             RingElemType<:AbsLocalizedRingElem{<:Any, <:Any, 
                                                <:MPolyComplementOfKPointIdeal
                                               }
            }
  return length(composition_series(M))
end
