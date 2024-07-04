export AbsAlgebraicCycle
export AlgebraicCycle
export algebraic_cycle
export components
export scheme


########################################################################
#
# AbsAlgebraicCycle
#
# Abstract type for algebraic cycles on a separated scheme X.
########################################################################

### The following declaration has been moved to src/forward_declarations.jl
#abstract type AbsAlgebraicCycle{
#                                CoveredSchemeType<:AbsCoveredScheme, 
#                                CoefficientRingType<:AbstractAlgebra.Ring
#                               }
#end

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

function coeff(D::AbsAlgebraicCycle, I::IdealSheaf)
  d = coefficient_dict(D)
  if I in keys(d)
    return d[I]
  else
    return zero(coefficient_ring(D))
  end
end

function is_effective(A::AbsAlgebraicCycle)
  return all(coeff(A, I)>=0 for I in components(A))
end

function Base.:<=(A::AbsAlgebraicCycle,B::AbsAlgebraicCycle)
  for I in components(A)
    coeff(A, I) <= coeff(B, I) || return false
  end
  for I in components(B)
    coeff(A, I) <= coeff(B, I) || return false
  end
  return true
end
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
        is_equidimensional(D) || error("components of a divisor must be sheaves of equidimensional ideals")
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

function zero(D::AbsAlgebraicCycle) 
  return AlgebraicCycle(scheme(D), coefficient_ring(D))
end

# provide non-camelcase methods
@doc raw"""
    algebraic_cycle(X::AbsCoveredScheme, R::Ring) -> AlgebraicCycle

Return the zero `AlgebraicCycle` over `X` with coefficients
in `R`.

# Examples
```jldoctest
julia> P, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> I = ideal([x^3-y^2*z]);

julia> Y = proj(P, I);

julia> Ycov = covered_scheme(Y);

julia> R = ZZ;

julia> algebraic_cycle(Ycov, R)
Zero algebraic cycle
  on scheme over QQ covered with 3 patches
with coefficients in integer ring
```
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

@doc raw"""
    algebraic_cycle(I::IdealSheaf, R::Ring) -> AlgebraicCycle

Return the `AlgebraicCycle` ``D = 1 ⋅ V(I)`` with coefficients
in ``R`` for a sheaf of prime ideals ``I``.

# Examples
```jldoctest
julia> P, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> I = ideal([x^3-y^2*z]);

julia> Y = proj(P);

julia> II = IdealSheaf(Y, I);

julia> R = ZZ;

julia> algebraic_cycle(II, R)
Effective algebraic cycle
  on scheme over QQ covered with 3 patches
with coefficients in integer ring
given as the formal sum of
  1 * sheaf of ideals

```
"""
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

@doc raw"""
    algebraic_cycle(I::IdealSheaf) -> AlgebraicCycle

Return the `AlgebraicCycle` ``D = 1 ⋅ V(I)`` with coefficients
in ``ℤ`` for a sheaf of prime ideals ``I``.

# Examples
```jldoctest
julia> P, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> I = ideal([x^3-y^2*z]);

julia> Y = proj(P);

julia> II = IdealSheaf(Y, I);

julia> R = ZZ;

julia> algebraic_cycle(II, R)
Effective algebraic cycle
  on scheme over QQ covered with 3 patches
with coefficients in integer ring
given as the formal sum of
  1 * sheaf of ideals
```
"""
algebraic_cycle(I::IdealSheaf) = AlgebraicCycle(I)

### copy constructor
function copy(D::AlgebraicCycle) 
  new_dict = IdDict{IdealSheaf, elem_type(coefficient_ring_type(D))}()
  for I in keys(coefficient_dict(D))
    new_dict[I] = D[I]
  end
  return AlgebraicCycle(scheme(D), coefficient_ring(D), new_dict)
end

###############################################################################
#
#  Printing
#
###############################################################################

# Method quite nice: if it knowns something about the algebraic cycle, it will
# tell us! (Irreducibility, effectivity, zero cycle,...)
#
# As it is given as a formal sum, we want a nice printing where all the
# coefficients of the respective ideal sheaves are aligned on the right - one
# needs to take care about some left offsets.
function Base.show(io::IO, ::MIME"text/plain", D::AlgebraicCycle)
  io = pretty(io)
  X = scheme(D)
  # If the IO context knows about a covering to be used, we use this one.
  # Otherwise, we check whether X has a simplified covering. If not, we use the
  # default covering of X
  cov = Oscar._covering_for_printing(io, X)
  eff = all(i >= 0 for i in collect(values(D.coefficients)))
  if length(components(D)) == 0
    print(io, "Zero algebraic cycle")
  else
    if eff
      print(io, "Effective algebraic cycle")
    else
      print(io, "Algebraic cycle")
    end
    if has_name(D)
      print(io, " ", get_attribute(D, :name))
    end
    if has_attribute(D, :dim)
      print(io, " of dimension $(dim(D))")
    end
  end
  println(io)
  print(io, Indent(), "on ", Lowercase())
  show(IOContext(io, :covering => cov), X)
  println(io, Dedent())
  print(io, "with coefficients in ", Lowercase(), coefficient_ring(D))
  if length(components(D)) != 0
    println(io, Dedent())
    print(io, Dedent(), "given as the formal sum of")
    print(io, Indent())
    co_str = String["$(D[I])" for I in components(D)]
    k = max(length.(co_str)...)
    for i in 1:length(components(D))
      println(io)
      I = components(D)[i]
      kI = length(co_str[i])
      print(io, " "^(k-kI)*"$(D[I]) * ")
      print(io, Indent(), Lowercase())
      show(IOContext(io, :show_scheme => false), I)
      print(io, Dedent())
    end
  end
  print(io, Dedent())
end

function Base.show(io::IO, D::AlgebraicCycle)
  io = pretty(io)
  X = scheme(D)
  eff = all(i >= 0 for i in collect(values(D.coefficients)))
  if length(components(D)) == 1
    prim = D[components(D)[1]] == 1 ? true : false
  else
    prim = false
  end
  if has_name(D)
    print(io, name(D))
  elseif get(io, :supercompact, false)
    print(io, "Algebraic cycle")
  elseif length(components(D)) == 0
    print(io, "Zero algebraic cycle on ", Lowercase(), scheme(D))
  elseif eff
    if prim
      print(io, "Irreducible algebraic cycle on ", Lowercase(), scheme(D))
    else
      print(io, "Effective algebraic cycle on ", Lowercase(), scheme(D))
    end
  else
    print(io, "Algebraic cycle on ", Lowercase(), scheme(D))
  end
end


@attr function dim(D::AlgebraicCycle)
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

@doc raw"""
    irreducible_decomposition(D::AbsAlgebraicCycle)

Returns a divisor ``E`` equal to ``D`` but as a formal sum ``E = ∑ₖ aₖ ⋅ Iₖ``
where the `components` ``Iₖ`` of ``E`` are all sheaves of prime ideals.
"""
function irreducible_decomposition(D::AbsAlgebraicCycle)
  all(I->is_prime(I), keys(coefficient_dict(D))) && return D
  result = zero(D)
  for (I, a) in coefficient_dict(D)
    next_dict = IdDict{IdealSheaf, elem_type(coefficient_ring(D))}()
    decomp = maximal_associated_points(I)
    for P in decomp
      k = _colength_in_localization(I, P)
      next_dict[P] = coefficient_ring(D)(k)
    end
    result = result + a * AlgebraicCycle(scheme(D), coefficient_ring(D), next_dict, check=false)
  end
  return result
end

function _colength_in_localization(Q::IdealSheaf, P::IdealSheaf; covering=simplified_covering(scheme(P)))
  X = scheme(Q)
  X === scheme(P) || error("ideal sheaves do not live on the same scheme")
  n = minimum([ngens(OO(U)) for U in patches(covering) if !isone(P(U))])
  j = findfirst(U->(!isone(P(U)) && ngens(OO(U))==n), patches(covering))
  U = patches(covering)[j]
  QU = Q(U)
  PU = P(U)
  W, loc_map = localization(OO(U), complement_of_prime_ideal(saturated_ideal(P(U))))
  Q_loc = loc_map(QU)
  F = free_module(W, 1)
  M, _ = quo(F, sub(F, [g*F[1] for g in gens(Q_loc)])[1])
  return length(M)
end

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
    # Make sure all generators are actually prime so that they can be compared. 
    all(is_prime, keys(coefficient_dict(D))) || return irreducible_decomposition(D) == E
    all(is_prime, keys(coefficient_dict(E))) || return D == irreducible_decomposition(E)

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
    @check begin
      dim(I) == 0 || continue
    end
    result = result + W[I]*colength(I)
  end
  return result
end

