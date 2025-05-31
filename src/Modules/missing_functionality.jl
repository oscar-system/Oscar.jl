########################################################################
#
# This file contains hacks for functionality that was found missing. 
#
# These methods should be improved and/or moved to their appropriate 
# places, eventually.
#

# the default module ordering assumes that we're computing in a global ring
function default_ordering(F::FreeMod{T}) where {T<:MPolyLocRingElem}
  return default_ordering(base_ring_module(F))
end

# missing functionality to write an element f ∈ I of an ideal as 
# a linear combination of the generators of I
function coordinates(f::MPolyRingElem, I::MPolyIdeal)
  iszero(f) && return zero_matrix(base_ring(I), 1, ngens(I))
  R = parent(f)
  R == base_ring(I) || error("polynomial does not belong to the base ring of the ideal")
  f in I || error("polynomial does not belong to the ideal")
  singular_assure(I)
  Rsing = I.gens.gens.Sx
  fsing = Singular.Ideal(Rsing, [Rsing(f)])
  a_s, u_s = Singular.lift(singular_generators(I), fsing)
  A_s = Matrix(a_s)
  U_s = Matrix(u_s)
  (ncols(U_s) == nrows(U_s) == 1 && iszero(U_s[1,1])) || error("no suitable ordering was used")
  A = zero_matrix(R, 1, ngens(I))
  for i in 1:ngens(I)
    A[1, i] = R(A_s[i, 1])
  end
  return A
end

### This is a dirty hack to bring the `coordinates` command to 
# MPolyQuoRing ideals. Should be replaced by something better!
# I tried to at least cache the ideal K via attributes, but 
# even that is not possible.
function coordinates(f::MPolyQuoRingElem, I::MPolyQuoIdeal)
  Q = base_ring(I)
  f in I || error("element is not in the ideal")
  J = modulus(Q)
  R = base_ring(Q)
  K = ideal(R, vcat(lift.(gens(I)), gens(J)))
  return Q.(coordinates(lift(f), K)[1:1, 1:ngens(I)])
end

function lift(f::MPolyRingElem, I::MPolyIdeal, o::MonomialOrdering)
  iszero(f) && return zero_matrix(base_ring(I), 1, ngens(I))
  R = parent(f)
  R == base_ring(I) || error("polynomial does not belong to the base ring of the ideal")
  Rsing = singular_poly_ring(R, o)
  fsing = Singular.Ideal(Rsing, [Rsing(f)])
  gsing = Singular.Ideal(Rsing, Rsing.(gens(I)))
  a_s, rem_s, u_s = lift(gsing, fsing, false, false, false)
  A_s = Matrix(a_s)
  u = R(u_s[1,1])
  A = zero_matrix(R, 1, ngens(I))
  for i in 1:ngens(I)
    A[1, i] = R(A_s[i, 1])
  end
  return A, u
end

### TODO: The following should not be necessary in the first place! 
# If the module code is supposed to run over arbitrary rings, it also 
# has to be possible to do without orderings. Up to now, defining 
# a FreeMod over an MPolyQuoLocRing requires me to implement this! Why????
#=function default_ordering(F::FreeMod{T}) where {T<:MPolyQuoLocRingElem}
  return default_ordering(base_ring_module(F))
end=#
    
subquo_type(::Type{RingType}) where {RingType<:Ring} = SubquoModule{elem_type(RingType)}
subquo_type(R::RingType) where {RingType<:Ring} = subquo_type(typeof(R))

function (==)(M1::SubquoModule{T}, M2::SubquoModule{T}) where {T<:AbsLocalizedRingElem}
  F = ambient_free_module(M1)
  F === ambient_free_module(M2) || error("ambient free modules are incompatible")
  all(x->iszero(M1(x)), relations(M2)) || return false
  all(x->iszero(M2(x)), relations(M1)) || return false
  all(x->x in M1, ambient_representatives_generators(M2)) || return false
  all(x->x in M2, ambient_representatives_generators(M1)) || return false
  return true
end


#to bypass the vec(collect(M)) which copies twice
function _vec(M::Generic.Mat)
  return vec(M.entries)
end

function _vec(M::MatElem)
  r = elem_type(base_ring(M))[]
  sizehint!(r, nrows(M) * ncols(M))
  for j=1:ncols(M)
    for i=1:nrows(M)
      push!(r, M[i, j])
    end
  end
  return r
end
