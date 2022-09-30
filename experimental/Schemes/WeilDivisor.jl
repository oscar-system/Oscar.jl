export WeilDivisor
export scheme_type, ideal_sheaf_type, coefficient_ring_type, coefficient_type
export scheme, components, coefficient_dict, coefficient_ring
export in_linear_system

export LinearSystem 
export divisor, find_subsystem

@Markdown.doc """
    WeilDivisor{
      CoveredSchemeType<:AbsCoveredScheme, 
      IdealSheafType<:IdealSheaf, 
      CoefficientRingType<:AbstractAlgebra.Ring, 
      CoefficientRingElemType<:AbstractAlgebra.RingElem
    }

A Weil divisor on a covered scheme ``X`` of type 
`CoveredSchemeType` as a formal linear combination of 
ideal sheaves of type `IdealSheafType` with coefficients 
in a ring ``R`` of type `CoefficientRingType`.
"""
@attributes mutable struct WeilDivisor{
    CoveredSchemeType<:AbsCoveredScheme, 
    CoefficientRingType<:AbstractAlgebra.Ring, 
    CoefficientRingElemType<:AbstractAlgebra.RingElem
   }
  X::CoveredSchemeType # the parent
  R::CoefficientRingType # the ring of coefficients
  coefficients::Dict{IdealSheaf, CoefficientRingElemType} # the formal linear combination

  function WeilDivisor(
      X::AbsCoveredScheme,
      R::CoefficientRingType, 
      coefficients::Dict{<:IdealSheaf, CoefficientRingElemType}
    ) where {CoefficientRingType, CoefficientRingElemType}
    # TODO: Do we want to require that the different effective divisors 
    # have the same underlying covering? Probably not.
    for D in keys(coefficients)
      scheme(D) == X || error("component of divisor does not lay in the given scheme")
      parent(coefficients[D]) == R || error("coefficients do not lay in the given ring")
    end
    return new{typeof(X), CoefficientRingType, CoefficientRingElemType}(X, R, coefficients)
  end
end

### type getters 
scheme_type(D::WeilDivisor{S, U, V}) where{S, T, U, V} = S
scheme_type(::Type{WeilDivisor{S, U, V}}) where{S, T, U, V} = S
coefficient_ring_type(D::WeilDivisor{S, U, V}) where{S, T, U, V} = U
coefficient_ring_type(::Type{WeilDivisor{S, U, V}}) where{S, T, U, V} = U
coefficient_type(D::WeilDivisor{S, U, V}) where{S, T, U, V} = V
coefficient_type(::Type{WeilDivisor{S, U, V}}) where{S, T, U, V} = V

### getter methods
scheme(D::WeilDivisor) = D.X
getindex(D::WeilDivisor, I::IdealSheaf) = (D.coefficients)[I]
components(D::WeilDivisor) = [ Z for Z in keys(D.coefficients)]
coefficient_dict(D::WeilDivisor) = D.coefficients
coefficient_ring(D::WeilDivisor) = D.R

set_name!(X::WeilDivisor, name::String) = set_attribute!(X, :name, name)
name_of(X::WeilDivisor) = get_attribute(X, :name)::String
has_name(X::WeilDivisor) = has_attribute(X, :name)

function setindex!(D::WeilDivisor, c::RingElem, I::IdealSheaf)
  parent(c) == coefficient_ring(D) || error("coefficient does not belong to the correct ring")
  coefficient_dict(D)[I] = c
end

function WeilDivisor(X::CoveredScheme, R::Ring)
  D = Dict{IdealSheaf, elem_type(R)}()
  return WeilDivisor(X, R, D)
end

function WeilDivisor(I::IdealSheaf, R::Ring)
  D = WeilDivisor(scheme(I), R)
  D[I] = one(R)
  return D
end

function WeilDivisor(I::IdealSheaf)
  D = WeilDivisor(scheme(I), ZZ)
  D[I] = one(ZZ)
  return D
end

### copy constructor
function copy(D::WeilDivisor) 
  new_dict = Dict{IdealSheaf, elem_type(coefficient_ring_type(D))}()
  for I in keys(coefficient_dict(D))
    new_dict[I] = D[I]
  end
  return WeilDivisor(scheme(D), coefficient_ring(D), new_dict)
end

function Base.show(io::IO, D::WeilDivisor)
  if has_name(D)
    print(io, name_of(D))
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

function +(D::T, E::T) where {T<:WeilDivisor}
  X = scheme(D)
  X == scheme(E) || error("divisors do not live on the same scheme")
  R = coefficient_ring(D)
  R == coefficient_ring(E) || error("coefficient rings do not coincide")
  C = copy(D)
  for I in components(E)
    if haskey(coefficient_dict(C), I)
      c = C[I] + E[I]
      if iszero(c) 
        delete!(coefficient_dict(C), I)
      else 
        C[I] = c
      end
    else
      coefficient_dict(C)[I] = E[I]
    end
  end
  return C
end

function -(D::T) where {T<:WeilDivisor}
  E = copy(D)
  for I in components(E)
    E[I] = -E[I]
  end
  return E
end

-(D::T, E::T) where {T<:WeilDivisor} = D + (-E)

function *(a::RingElem, E::WeilDivisor)
  c = coefficient_ring(E)()
  D = copy(E)
  parent(a) == coefficient_ring(E) ? (c = a)::elem_type(coefficient_ring(E)) : c = coefficient_ring(E)(a)
  for I in components(D)
    D[I] = c*D[I]
  end
  return D
end

*(a::Int, E::WeilDivisor) = coefficient_ring(E)(a)*E
*(a::Integer, E::WeilDivisor) = coefficient_ring(E)(a)*E

+(D::WeilDivisor, I::IdealSheaf) = D + WeilDivisor(I)

function ==(D::WeilDivisor, E::WeilDivisor) 
  keys(coefficient_dict(D)) == keys(coefficient_dict(E)) || return false
  for I in keys(coefficient_dict(D))
    D[I] == E[I] || return false
  end
  return true
end

#function intersection(D::T, E::T) where {T<:WeilDivisor}
#  X = scheme(D)
#  X == scheme(E) || error("divisors do not live on the same scheme")
#  R = coefficient_ring(D)
#  R == coefficient_ring(E) || error("divisors do not have the same coefficient ring")
#  # prepare a copy of the divisors
#  D_copy = WeilDivisor(X, R)
#  E_copy = WeilDivisor(X, R)
#  # check whether a common refinement of the covering is necessary
#  CD = covering(D)
#  CE = covering(E)
#  if CD != CE
#    CC, f, g = common_refinement(X, CD, CE)
#    D_copy = pullback(f, D)
#    E_copy = pullback(g, E)
#  else
#    D_copy = D
#    E_copy = E
#  end
#  # TODO: Work out the intersection
#end

@Markdown.doc """
    in_linear_system(f::VarietyFunctionFieldElem, D::WeilDivisor; check::Bool=true)

Returns `true` if the rational function `f` is in the linear system ``|D|``
and `false` otherwise.
"""
function in_linear_system(f::VarietyFunctionFieldElem, D::WeilDivisor; check::Bool=true)
  X = scheme(D) 
  X == variety(parent(f)) || error("schemes not compatible")
  C = default_covering(X)
  for I in components(D)
    order_on_divisor(f, I, check=check) >= -D[I] || return false
  end
  for U in patches(C)
    # we have to check that f[U] has no poles outside D[U]
  end
  return true
end

@Markdown.doc """
    LinearSystem{DivisorType<:WeilDivisor}

A linear system of a Weil divisor `D` on a scheme `X`, 
generated by rational functions ``f₁,…,fᵣ``.
"""
@attributes mutable struct LinearSystem{DivisorType<:WeilDivisor}
  D::DivisorType
  f::Vector{<:VarietyFunctionFieldElem}

  function LinearSystem(f::Vector, D::WeilDivisor; check::Bool=true)
    length(f) == 0 && return new{typeof(D)}(D, Vector{VarietyFunctionFieldElem}())
    KK = parent(f[1])
    all(g -> (parent(g) == KK), f[2:end]) || error("elements must have the same parent")
    X = scheme(D)
    X == variety(KK) || error("input not compatible")

    if check
      all(g->in_linear_system(g, D), f) || error("element not in linear system")
    end
    f = Vector{VarietyFunctionFieldElem}(f)
    return new{typeof(D)}(D, f)
  end
end

### essential getters 
@Markdown.doc """
    weil_divisor(L::LinearSystem)

Returns `D` on a linear system `L = |D|`.
"""
function weil_divisor(L::LinearSystem) 
  return L.D
end
gens(L::LinearSystem) = L.f
ngens(L::LinearSystem) = length(L.f)
@Markdown.doc """
    variety(L::LinearSystem)

Returns the variety on which `L` is defined.
"""
variety(L::LinearSystem) = scheme(weil_divisor(L))

@Markdown.doc """
    find_subsystem(L::LinearSystem, P::IdealSheaf, n::Int)

Given a linear system ``L = |D|``, a sheaf of prime ideals `P` 
and an integer `n`, this returns a pair ``(K, A)`` consisting 
of the subsystem of elements in ``|D + P|`` and the representing 
matrix ``A`` for its inclusion into ``L`` on the given set 
of generators.
"""
function find_subsystem(L::LinearSystem, P::IdealSheaf, n::Int)
  # find one chart in which P is supported
  # TODO: There might be preferred choices for charts with 
  # the least complexity.
  X = variety(L)
  X == scheme(P) || error("input incompatible")
  C = default_covering(X)
  U = first(patches(C))
  for V in patches(C)
    if !(one(OO(V)) in P[V])
      U = V
      break
    end
  end
  # Now U it is.

  # Assemble the local representatives
  R = ambient_ring(U)
  loc_rep = [g[U] for g in gens(L)]
  common_denominator = gcd([denominator(g) for g in loc_rep])
  numerators = [numerator(g)*divexact(common_denominator, denominator(g)) for g in loc_rep]
  RP, _ = Localization(R, complement_of_ideal(saturated_ideal(P[U])))
  PP = RP(prime_ideal(inverted_set(RP)))
  denom_mult = (_minimal_power_such_that(PP, I -> !(RP(common_denominator) in I))[1])-1
  w = n + denom_mult # Ajust!
  pPw = saturated_ideal(PP^w) # the symbolic power
  images = elem_type(R)[]
  for a in numerators
    push!(images, normal_form(a, pPw))
  end

  # collect a monomial basis in which to represent the results
  all_mons = elem_type(R)[]
  for b in images
    all_mons = vcat(all_mons, [m for m in monomials(b) if !(b in all_mons)])
  end

  kk = base_ring(X)
  A = zero(MatrixSpace(kk, ngens(L), length(all_mons)))
  for i in 1:ngens(L)
    for (c, m) in zip(coefficients(images[i]), monomials(images[i]))
      k = findfirst(x->(x==m), all_mons)
      A[i, k] = c
    end
  end

  r, K = left_kernel(A)
  new_gens = [sum([K[i,j]*gens(L)[j] for j in 1:ncols(K)]) for i in 1:nrows(K)]
  return LinearSystem(new_gens, weil_divisor(L) + n*WeilDivisor(P)), K
end

