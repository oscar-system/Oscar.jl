########################################################################
#
# This file contains hacks for functionality that was found missing. 
#
# These methods should be improved and/or moved to their appropriate 
# places, eventually.
#

# iterators over singular modules
Base.iterate(L::Singular.smodule) = iterate(L, 1)
Base.eltype(::Type{Singular.smodule}) = Singular.svector
Base.length(L::Singular.smodule) = ngens(L)

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
  Rsing = I.gens.Sx
  fsing = Singular.Ideal(Rsing, [Rsing(f)])
  a_s, u_s = Singular.lift(I.gens.S, fsing)
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
  return Q.(coordinates(lift(f), K)[1, 1:ngens(I)])
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

function _kbase(I::SubModuleOfFreeModule{<:MPolyRingElem{<:FieldElem}};
    ordering::ModuleOrdering=default_ordering(ambient_free_module(I))
  )
  lead_I = leading_module(I, ordering)
  @show ordering
  @show lead_I
  F = ambient_free_module(I)
  R = base_ring(F)
  # Iterate through all components until we have everything
  result = elem_type(F)[]
  for i in 1:ngens(F)
    m = ideal(R, gens(R))
    pom = ideal(R, one(R))
    while true
      new_mons = [x for x in gens(pom) if !(x*F[i] in lead_I)]
      iszero(length(new_mons)) && break
      result = vcat(result, [x*F[i] for x in new_mons])
      pom = pom * m
    end
  end
  return result
end

function vector_space(
    K::AbstractAlgebra.Field, M::SubquoModule{<:MPolyRingElem{<:FieldElem}};
    ordering::ModuleOrdering=default_ordering(ambient_free_module(M))
  )
  R = base_ring(M)
  @assert K === coefficient_ring(R)

  C = presentation(M)
  F0 = C[0]
  F1 = C[1]
  # TODO: Is there a way to access this directly? If yes, do so.
  I = SubModuleOfFreeModule(F0, map(C, 1).(gens(F1)))
  gb_I = standard_basis(I, ordering=ordering)
  mons = _kbase(I, ordering=ordering)
  l = length(mons)
  V = free_module(K, l)
  function im(a::Generic.FreeModuleElem)
    @assert parent(a) == V
    b = zero(F0)
    for k=1:l
      c = a[k]
      if !iszero(c)
        b += c*mons[k]
      end
    end
    return map(C, 0)(b)
  end

  # The inverse function. We use the fact that for a chosen monomial ordering 
  # the monomials which are not in the leading ideal, form a basis for the 
  # quotient; see Greuel/Pfister "A singular introduction to Commutative Algebra".
  function prim(a::SubquoModuleElem)
    @assert parent(a) === M
    b = preimage(map(C, 0), a)
    b = normal_form(b, gb_I) #TODO: Is this normal_form reduced?
    result = zero(V)
    while !iszero(b)
      m = leading_monomial(b, ordering=ordering)
      c = leading_coefficient(b, ordering=ordering)
      j = findfirst(n->n==m, mons)
      result = result + c * V[j]
      b = b - c * m
    end
    return result
  end
  # Adaptation of the hack in src/Rings/mpolyquo-localizations.jl
  function prim_loc(a::SubquoModuleElem)
    b = preimage(map(C, 0), a)
    result = zero(V)
    while !iszero(b)
      m = leading_monomial(b, ordering=ordering)
      c = leading_coefficient(b, ordering=ordering)
      error("normal form is buggy for local orderings at the moment; see #2155")
      t = normal_form(c*m, gb_I)
      if m in leading_module(I, ordering) 
        b = b - c*m + t
      else
        j = findfirst(n->n==m, mons)
        result = result + c * V[j]
        b = b - c * m
      end
    end
    return map(C, 0)(result)
  end
  if is_global(induced_ring_ordering(ordering))
    return V, MapFromFunc(im, prim, V, M)
  else
    return V, MapFromFunc(im, prim_loc, V, M)
  end
end

function vector_space(
    K::AbstractAlgebra.Field, M::SubquoModule{T}
  ) where {T<:MPolyLocRingElem{<:Any, <:Any, <:Any, <:Any, <:MPolyComplementOfKPointIdeal}}
  L = base_ring(M)::MPolyLocRing
  R = base_ring(L)::MPolyRing
  kk = coefficient_ring(R)
  kk === K || error("change of fields not implemented")
  shift, back_shift = base_ring_shifts(L)
  Mb = pre_saturated_module(M)
  Mb, m_shift, m_backshift = shifted_module(M)
  Fb = ambient_free_module(Mb)
  V, iso = vector_space(kk, Mb, ordering=lex(gens(Fb))*negdegrevlex(gens(R)))
  error("implementation is still work in progress")
  return V, compose(iso, m_backshift)
end

