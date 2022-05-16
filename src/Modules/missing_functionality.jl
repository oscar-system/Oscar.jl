export subquo_type
########################################################################
#
# This file contains hacks for functionality that was found missing. 
#
# These methods should be improved and/or moved to their appropriate 
# places, eventually.
#

# missing functionality for maps of modules
compose(f::FreeModuleHom, g::FreeModuleHom) = hom(domain(f), codomain(g), representing_matrix(f)*representing_matrix(g))
compose(f::SubQuoHom, g::FreeModuleHom) = hom(domain(f), codomain(g), representing_matrix(f)*representing_matrix(g))
compose(f::SubQuoHom, g::SubQuoHom) = hom(domain(f), codomain(g), representing_matrix(f)*representing_matrix(g))
compose(f::FreeModuleHom, g::SubQuoHom) = hom(domain(f), codomain(g), representing_matrix(f)*representing_matrix(g))

# missing (?) constructors for SubQuos
function SubQuo(F::FreeMod{T}, g::Vector{FreeModElem{T}}, q::Vector{FreeModElem{T}}) where {T<:RingElem} 
  return SubQuo(Oscar.SubModuleOfFreeModule(F, g), Oscar.SubModuleOfFreeModule(F, q))
end

function sub(F::FreeMod{T}, A::MatElem{T}) where {T} 
  M = SubQuo(F, A, zero(MatrixSpace(base_ring(F), 1, rank(F))))
  inc = hom(M, F, ambient_representatives_generators(M))
  inc.matrix = A
  return M, inc
end

function quo(F::FreeMod{T}, A::MatElem{T}) where {T}
  E = one(MatrixSpace(base_ring(F), rank(F), rank(F)))
  M = SubQuo(F, E, A)
  proj = hom(F, M, gens(M))
  proj.matrix = E
  return M, proj
end

#promotion for scalar multiplication
AbstractAlgebra.promote_rule(::Type{RET}, ::Type{MET}) where {RET<:RingElem, MET<:ModuleElem} = MET

# iterators over singular modules
Base.iterate(L::Singular.smodule) = iterate(L, 1)
Base.eltype(::Type{Singular.smodule}) = Singular.svector
Base.length(L::Singular.smodule) = ngens(L)

# the default module ordering assumes that we're computing in a global ring
function default_ordering(F::FreeMod{T}) where {T<:MPolyLocalizedRingElem}
  return default_ordering(base_ring_module(F))
end

# missing functionality to write an element f ∈ I of an ideal as 
# a linear combination of the generators of I
function coordinates(f::MPolyElem, I::MPolyIdeal)
  R = parent(f)
  R == base_ring(I) || error("polynomial does not belong to the base ring of the ideal")
  f in I || error("polynomial does not belong to the ideal")
  singular_assure(I)
  Rsing = I.gens.Sx
  fsing = Singular.Ideal(Rsing, [Rsing(f)])
  a_s, u_s = Singular.lift(I.gens.S, fsing)
  A_s = Matrix(a_s)
  U_s = Matrix(u_s)
  (ncols(U_s) == nrows(U_s) == 1 && iszero(U_s[1,1])) || error("no suitable ordering was used")
  A = zero(MatrixSpace(R, 1, ngens(I)))
  for i in 1:ngens(I)
    A[1, i] = R(A_s[i, 1])
  end
  return A
end

### TODO: The following should not be necessary in the first place! 
# If the module code is supposed to run over arbitrary rings, it also 
# has to be possible to do without orderings. Up to now, defining 
# a FreeMod over an MPolyQuoLocalizedRing requires me to implement this! Why????
function default_ordering(F::FreeMod{T}) where {T<:MPolyQuoLocalizedRingElem}
  return default_ordering(base_ring_module(F))
end
    
subquo_type(::Type{RingType}) where {RingType<:Ring} = SubQuo{elem_type(RingType)}
subquo_type(R::RingType) where {RingType<:Ring} = subquo_type(typeof(R))
