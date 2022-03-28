export default_ordering, singular_assure 

function singular_assure(M::ModuleGens{T}) where {T<:MPolyLocalizedRingElem{<:Any, <:Any, <:Any, <:Any, <:MPolyComplementOfKPointIdeal}}
  L = base_ring(M)
  R = base_ring(L)
  F = ambient_free_module(M)
  SR = singular_ring(R, singular(default_ordering(F)))
  SF = Singular.FreeModule(SR, ngens(F))
  sgens_and_denoms = [SF(a) for a in oscar_generators(M)]
  M.SF = SF
  M.S = Singular.Module(SR, [a[1] for a in sgens_and_denoms]...)
end

function (SF::Singular.FreeMod)(v::FreeModElem{T}) where {T<:MPolyLocalizedRingElem{<:Any, <:Any, <:Any, <:Any, <:MPolyComplementOfKPointIdeal}}
  F = parent(v)
  L = base_ring(F)
  R = base_ring(L)
  SR = base_ring(SF)
  # find a common denominator for the components of this generator
  d = lcm(denominator.(Vector(v)))
  new_numerators = [numerator(v[i])*divexact(d, denominator(v[i])) for i in 1:ngens(F)]
  return sum([SR(new_numerators[i])*gens(SF)[i] for i in 1:length(new_numerators)]), SR(d)
end

function default_ordering(F::FreeModuleType) where {FreeModuleType<:FreeMod{<:MPolyLocalizedRingElem{<:Any, <:Any, <:Any, <:Any, <:MPolyComplementOfKPointIdeal}}}
  # We need to set up a free module over the polynomial ring 
  # so that the monomial ordering can be given.
  L = base_ring(F)
  R = base_ring(L)
  helperF = FreeMod(R, ngens(F))
  return negdegrevlex(gens(base_ring(base_ring(F))))*lex(gens(helperF))
end

#function ModuleGens{T}(O::Vector{FreeModElemType}, 
#                       F::FreeMod{T}
#                      ) where {
#                               T<:MPolyLocalizedRingElem{<:Any, <:Any, <:Any, <:Any, 
#                                                         <:MPolyComplementOfKPointIdeal},
#                               FreeModElemType <: FreeModElem
#                              }
#  return ModuleGens(O, F)
#end
