export AbsAlgebraicCycle
export AlgebraicCycle
export components
export scheme



########################################################################
#
# AbsAlgebraicCycle
#
# Abstract type for algebraic cycles on a separated scheme X.
########################################################################
abstract type AbsAlgebraicCycle{
                                CoveredSchemeType<:AbsCoveredScheme, 
                                CoefficientRingType<:AbstractAlgebra.Ring
                               }
end

### type getters 
scheme_type(D::AbsAlgebraicCycle{S, U}) where {S, U} = S
scheme_type(::Type{AbsAlgebraicCycle{S, U}}) where {S, U} = S
coefficient_ring_type(D::AbsAlgebraicCycle{S, U}) where {S, U} = U
coefficient_ring_type(::Type{AbsAlgebraicCycle{S, U}}) where {S, U} = U
coefficient_ring_elem_type(D::AbsAlgebraicCycle{S, U}) where {S, U} = elem_type(U)
coefficient_ring_elem_type(::Type{AbsAlgebraicCycle{S, U}}) where {S, U} = elem_type(U)

### essential getters and functionality

@doc raw"""
    scheme(D::AbsAlgebraicCycle)

Return the `CoveredScheme` ``X`` on which `D` is defined.
"""
scheme(D::AbsAlgebraicCycle) = scheme(underlying_cycle(D))

# For an element `I` of `components(D)`, this returns the coefficient 
# of `I` in the formal sum for `D`.
getindex(D::AbsAlgebraicCycle, I::IdealSheaf) = getindex(underlying_cycle(D), I)

@doc raw"""
    components(D::AbsAlgebraicCycle)

Return the irreducible components ``Eⱼ`` of the divisor 
``D = Σⱼ aⱼ ⋅ Eⱼ``.
"""
components(D::AbsAlgebraicCycle) = components(underlying_cycle(D))

# Return the coefficient ring over which the cycle is defined
coefficient_ring(D::AbsAlgebraicCycle) = coefficient_ring(underlying_cycle(D))

# All `components` of a cycle `D` must be prime. This returns the supremum 
# of their dimensions.
dim(D::AbsAlgebraicCycle) = dim(underlying_cycle(D))

set_name!(X::AbsAlgebraicCycle, name::String) = set_attribute!(X, :name, name)
name(X::AbsAlgebraicCycle) = get_attribute(X, :name)::String
has_name(X::AbsAlgebraicCycle) = has_attribute(X, :name)

function setindex!(D::AbsAlgebraicCycle, c::RingElem, I::IdealSheaf)
  parent(c) === coefficient_ring(D) || error("coefficient does not belong to the correct ring")
  return setindex!(underlying_cycle(D), c, I)
end

# Non user-facing getters
# We assume every cycle D = ∑ᵢ aᵢ ⋅ 𝒥 ᵢto store the data of its formal 
# sum in an `IdDict` with the 𝒥 ᵢ as keys and their coefficients aᵢ as 
# values. This is a getter to that dictionary to allow for a generic 
# implementation of the arithmetic. 
coefficient_dict(D::AbsAlgebraicCycle) = coefficient_dict(underlying_cycle(D))

### forwarding of the essential functionality

function underlying_cycle(D::AbsAlgebraicCycle)
  error("method `underlying_cycle` not implemented for arguments of type $(typeof(D))")
end


########################################################################
#
# AlgebraicCycle
#
# A minimal implementation of AbsAlgebraicCycle.
########################################################################
 
@attributes mutable struct AlgebraicCycle{
                                          CoveredSchemeType<:AbsCoveredScheme, 
                                          CoefficientRingType<:AbstractAlgebra.Ring, 
                                          CoefficientRingElemType<:AbstractAlgebra.RingElem
                                         } <: AbsAlgebraicCycle{CoveredSchemeType, 
                                                                CoefficientRingType}

  X::CoveredSchemeType # the parent
  R::CoefficientRingType # the ring of coefficients
  coefficients::IdDict{IdealSheaf, CoefficientRingElemType} # the formal linear combination

  function AlgebraicCycle(
      X::AbsCoveredScheme,
      R::CoefficientRingType, 
      coefficients::IdDict{<:IdealSheaf, CoefficientRingElemType};
      check::Bool=true
    ) where {CoefficientRingType, CoefficientRingElemType}
    for D in keys(coefficients)
      space(D) === X || error("component of divisor does not lie in the given scheme")
      parent(coefficients[D]) === R || error("coefficients do not lie in the given ring")
    end
    if check
      is_integral(X) || error("scheme must be integral") 
      #is_separated(X) || error("scheme must be separated") # We need to test this somehow, but how?
      for D in keys(coefficients)
        isprime(D) || error("components of a divisor must be sheaves of prime ideals")
      end
    end
    return new{typeof(X), CoefficientRingType, CoefficientRingElemType}(X, R, coefficients)
  end
end

### implementation of the essential functionality
scheme(D::AlgebraicCycle) = D.X
getindex(D::AlgebraicCycle, I::IdealSheaf) = (D.coefficients)[I]

components(D::AlgebraicCycle) = collect(keys(D.coefficients))
coefficient_dict(D::AlgebraicCycle) = D.coefficients
coefficient_ring(D::AlgebraicCycle) = D.R

set_name!(X::AlgebraicCycle, name::String) = set_attribute!(X, :name, name)
name(X::AlgebraicCycle) = get_attribute(X, :name)::String
has_name(X::AlgebraicCycle) = has_attribute(X, :name)

function setindex!(D::AlgebraicCycle, c::RingElem, I::IdealSheaf)
  parent(c) === coefficient_ring(D) || error("coefficient does not belong to the correct ring")
  coefficient_dict(D)[I] = c
end

@doc raw"""
    AlgebraicCycle(X::CoveredScheme, R::Ring)

Return the zero `AlgebraicCycle` over `X` with coefficients 
in `R`.
"""
function AlgebraicCycle(X::AbsCoveredScheme, R::Ring)
  D = IdDict{IdealSheaf, elem_type(R)}()
  return AlgebraicCycle(X, R, D)
end

# provide non-camelcase methods
@doc raw"""
    algebraic_cycle(X::AbsCoveredScheme, R::Ring)

See the documentation for `AlgebraicCycle`.
"""
algebraic_cycle(X::AbsCoveredScheme, R::Ring) = AlgebraicCycle(X, R)

@doc raw"""
    AlgebraicCycle(I::IdealSheaf, R::Ring)

Return the `AlgebraicCycle` ``D = 1 ⋅ V(I)`` with coefficients 
in ``R`` for a sheaf of prime ideals ``I``.
"""
function AlgebraicCycle(I::IdealSheaf, R::Ring)
  D = AlgebraicCycle(space(I), R)
  D[I] = one(R)
  return D
end

algebraic_cycle(I::IdealSheaf, R::Ring) = AlgebraicCycle(I, R)

@doc raw"""
    AlgebraicCycle(I::IdealSheaf)

Return the `AlgebraicCycle` ``D = 1 ⋅ V(I)`` with coefficients
in ``ℤ`` for a sheaf of prime ideals ``I``.
"""
function AlgebraicCycle(I::IdealSheaf)
  D = AlgebraicCycle(space(I), ZZ)
  D[I] = one(ZZ)
  return D
end

algebraic_cycle(I::IdealSheaf) = AlgebraicCycle(I)

### copy constructor
function copy(D::AlgebraicCycle) 
  new_dict = IdDict{IdealSheaf, elem_type(coefficient_ring_type(D))}()
  for I in keys(coefficient_dict(D))
    new_dict[I] = D[I]
  end
  return AlgebraicCycle(scheme(D), coefficient_ring(D), new_dict)
end

function Base.show(io::IO, D::AlgebraicCycle)
  if has_name(D)
    print(io, name(D))
    return
  end
  if length(components(D)) == 0
    println(io, "the zero Weil divisor on $(scheme(D))")
    return
  end
  println(io, "Weil divisor on $(scheme(D)) given as the formal sum")
  comp = ["$(D[I]) ⋅ $(I)" for I in components(D)]
  out_str = comp[1]
  for c in comp[2:end]
    out_str = out_str* " + " * c
  end
  println(io, out_str)
end

function dim(D::AlgebraicCycle)
  result = -1
  for I in components(D)
    d = dim(I)
    if d > result
      result = d
    end
  end
  return result
end



### half-generic implementation of the arithmetic
# Note that we need one minimal concrete type for the return values, 
# so the implementation can not be truly generic. 

function +(D::T, E::T) where {T<:AbsAlgebraicCycle}
  X = scheme(D)
  X === scheme(E) || error("divisors do not live on the same scheme")
  R = coefficient_ring(D)
  R === coefficient_ring(E) || error("coefficient rings do not coincide")
  dict = IdDict{IdealSheaf, elem_type(R)}()
  for I in keys(coefficient_dict(D))
    dict[I] = D[I]
  end
  for I in keys(coefficient_dict(E))
    if haskey(dict, I)
      c = D[I] + E[I]
      if iszero(c) 
        delete!(dict, I)
      else 
        dict[I] = c
      end
    else
      dict[I] = E[I]
    end
  end
  return AlgebraicCycle(X, R, dict, check=false)
end

function -(D::T) where {T<:AbsAlgebraicCycle}
  dict = IdDict{IdealSheaf, elem_type(coefficient_ring(D))}()
  for I in keys(coefficient_dict(D))
    dict[I] = -D[I]
  end
  return AlgebraicCycle(scheme(D), coefficient_ring(D), dict, check=false)
end

-(D::T, E::T) where {T<:AbsAlgebraicCycle} = D + (-E)

function *(a::RingElem, E::AbsAlgebraicCycle)
  dict = IdDict{IdealSheaf, elem_type(coefficient_ring(E))}()
  for I in keys(coefficient_dict(E))
    c = a*E[I]
    if iszero(c)
      delete!(dict, I)
    else
      dict[I] = c
    end
  end
  return AlgebraicCycle(scheme(E), coefficient_ring(E), dict, check=false)
end

*(a::Int, E::AbsAlgebraicCycle) = coefficient_ring(E)(a)*E
*(a::Integer, E::AbsAlgebraicCycle) = coefficient_ring(E)(a)*E

+(D::AbsAlgebraicCycle, I::IdealSheaf) = D + AbsAlgebraicCycle(I)

function ==(D::AbsAlgebraicCycle, E::AbsAlgebraicCycle) 
  if all(k->k in keys(coefficient_dict(D)), keys(coefficient_dict(E))) && all(k->k in keys(coefficient_dict(E)), keys(coefficient_dict(D))) 
    for I in keys(coefficient_dict(D))
      if haskey(coefficient_dict(E), I)
        D[I] == E[I] || return false
      else
        iszero(D[I]) || return false
      end
    end
    for I in keys(coefficient_dict(E))
      !(I in keys(coefficient_dict(D))) && !(iszero(E[I])) && return false
    end
  else
    keys_D = collect(keys(coefficient_dict(D)))
    keys_E = collect(keys(coefficient_dict(E)))
    for I in keys(coefficient_dict(D))
      I_cand = findall(x->(x==I), keys_D)
      J_cand = findall(x->(x==I), keys_E)
      sum([D[keys_D[i]] for i in I_cand]) == sum([E[keys_E[j]] for j in J_cand]) || return false
    end
    for J in keys(coefficient_dict(E))
      I_cand = findall(x->(x==J), keys_D)
      J_cand = findall(x->(x==J), keys_E)
      sum([D[keys_D[i]] for i in I_cand]) == sum([E[keys_E[j]] for j in J_cand]) || return false
    end
  end
  return true
end

@doc raw"""
    integral(W::AbsAlgebraicCycle)

Assume ``W`` is an algebraic cycle on ``X``. This returns the sum of 
the lengths of all the components of dimension `0` of ``W``.
"""
function integral(W::AbsAlgebraicCycle; check::Bool=true)
  result = zero(coefficient_ring(W))
  X = scheme(W)
  for I in components(W)
    # All components must be prime, so the information is exclusively stored in the coefficients
    if check
      dim(I) == 0 || continue
    end
    result = result + W[I]
  end
  return result
end


