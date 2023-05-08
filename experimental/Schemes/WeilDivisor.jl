export LinearSystem
export WeilDivisor
export coefficient_dict
export coefficient_ring
export coefficient_ring_type
export coefficient_type
export components
export divisor
export ideal_sheaf_type
export in_linear_system
export linear_system
export scheme
export scheme_type
export subsystem
export weil_divisor

@doc raw"""
    WeilDivisor

A Weil divisor on an integral separated `AbsCoveredScheme` ``X``; 
stored as a formal linear combination over some ring ``R`` of 
(prime) ideal sheaves on ``X``.
"""
@attributes mutable struct WeilDivisor{
    CoveredSchemeType<:AbsCoveredScheme, 
    CoefficientRingType<:AbstractAlgebra.Ring, 
    CoefficientRingElemType<:AbstractAlgebra.RingElem
   } <: AbsAlgebraicCycle{CoveredSchemeType, CoefficientRingType}
  C::AlgebraicCycle{CoveredSchemeType, CoefficientRingType, CoefficientRingElemType}

  function WeilDivisor(
      X::AbsCoveredScheme,
      R::CoefficientRingType, 
      coefficients::IdDict{<:IdealSheaf, CoefficientRingElemType};
      check::Bool=true
    ) where {CoefficientRingType, CoefficientRingElemType}
    if check
      for D in keys(coefficients)
        isprime(D) || error("components of a divisor must be sheaves of prime ideals")
        dim(X) - dim(D) == 1 || error("components of a divisor must be of codimension one")
      end
    end
    return new{typeof(X), CoefficientRingType, CoefficientRingElemType}(AlgebraicCycle(X, R, coefficients, check=check))
  end

  function WeilDivisor(C::AlgebraicCycle; check::Bool=true)
    X = scheme(C)
    if check
      for D in keys(coefficient_dict(C))
        isprime(D) || error("components of a divisor must be sheaves of prime ideals")
        dim(X) - dim(D) == 1 || error("components of a divisor must be of codimension one")
      end
    end
    return new{typeof(X), coefficient_ring_type(C), coefficient_ring_elem_type(C)}(C)
  end
end

### forwarding of all essential functionality
underlying_cycle(D::WeilDivisor) = D.C

@attr function dim(I::IdealSheaf)
  dims = [dim(I(U)) for U in affine_charts(scheme(I))]
  return maximum(dims)
end

### type getters 
scheme_type(D::WeilDivisor{S, U, V}) where{S, U, V} = S
scheme_type(::Type{WeilDivisor{S, U, V}}) where{S, U, V} = S
coefficient_ring_type(D::WeilDivisor{S, U, V}) where{S, U, V} = U
coefficient_ring_type(::Type{WeilDivisor{S, U, V}}) where{S, U, V} = U
coefficient_type(D::WeilDivisor{S, U, V}) where{S, U, V} = V
coefficient_type(::Type{WeilDivisor{S, U, V}}) where{S, U, V} = V


@doc raw"""
    WeilDivisor(X::CoveredScheme, R::Ring)

Return the zero `WeilDivisor` over `X` with coefficients 
in `R`.
"""
function WeilDivisor(X::AbsCoveredScheme, R::Ring)
  D = IdDict{IdealSheaf, elem_type(R)}()
  return WeilDivisor(X, R, D)
end

function zero(W::WeilDivisor)
  return WeilDivisor(scheme(W), coefficient_ring(W))
end

# provide non-camelcase methods
@doc raw"""
    weil_divisor(X::AbsCoveredScheme, R::Ring)

See the documentation for `WeilDivisor`.
"""
weil_divisor(X::AbsCoveredScheme, R::Ring) = WeilDivisor(X, R)

@doc raw"""
    WeilDivisor(I::IdealSheaf, R::Ring)

Return the `WeilDivisor` ``D = 1 ⋅ V(I)`` with coefficients 
in ``R`` for a sheaf of prime ideals ``I``.
"""
function WeilDivisor(I::IdealSheaf, R::Ring)
  D = WeilDivisor(space(I), R)
  D[I] = one(R)
  return D
end

weil_divisor(I::IdealSheaf, R::Ring) = WeilDivisor(I, R)

@doc raw"""
    WeilDivisor(I::IdealSheaf)

Return the `WeilDivisor` ``D = 1 ⋅ V(I)`` with coefficients
in ``ℤ`` for a sheaf of prime ideals ``I``.
"""
function WeilDivisor(I::IdealSheaf)
  D = WeilDivisor(space(I), ZZ)
  D[I] = one(ZZ)
  return D
end

weil_divisor(I::IdealSheaf) = WeilDivisor(I)

### copy constructor
function copy(D::WeilDivisor) 
  new_dict = IdDict{IdealSheaf, elem_type(coefficient_ring_type(D))}()
  for I in keys(coefficient_dict(D))
    new_dict[I] = D[I]
  end
  return WeilDivisor(scheme(D), coefficient_ring(D), new_dict)
end

function Base.show(io::IO, D::WeilDivisor)
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

function +(D::T, E::T) where {T<:WeilDivisor}
  return WeilDivisor(underlying_cycle(D) + underlying_cycle(E), check=false)
end

function -(D::T) where {T<:WeilDivisor}
  return WeilDivisor(-underlying_cycle(D), check=false)
end

-(D::T, E::T) where {T<:WeilDivisor} = D + (-E)

function *(a::RingElem, E::WeilDivisor)
  return WeilDivisor(a*underlying_cycle(E), check=false)
end

*(a::Int, E::WeilDivisor) = coefficient_ring(E)(a)*E
*(a::Integer, E::WeilDivisor) = coefficient_ring(E)(a)*E

+(D::WeilDivisor, I::IdealSheaf) = D + WeilDivisor(I)

function ==(D::WeilDivisor, E::WeilDivisor) 
  return underlying_cycle(D) == underlying_cycle(E)
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

@doc raw"""
    in_linear_system(f::VarietyFunctionFieldElem, D::WeilDivisor; check::Bool=true) -> Bool

Check if the rational function `f` is in the linear system ``|D|``.
"""
function in_linear_system(f::VarietyFunctionFieldElem, D::WeilDivisor; check::Bool=true)
  X = scheme(D) 
  X === variety(parent(f)) || error("schemes not compatible")
  C = default_covering(X)
  for I in components(D)
    order_on_divisor(f, I, check=check) >= -D[I] || return false
  end
  for U in patches(C)
    # we have to check that f[U] has no poles outside the support of D[U]
    g = numerator(f[U])
    h = denominator(f[U])
    J = prod([J(U) for J in components(D)])
    incH = ClosedEmbedding(U, J)
    W = complement(incH) # This is a SpecOpen
    is_regular(f, W) || return false
  end
  return true
end

@doc raw"""
    LinearSystem

A linear system of a Weil divisor `D` on a variety `X`, 
generated by rational functions ``f₁,…,fᵣ ∈ K(X)``.
"""
@attributes mutable struct LinearSystem{DivisorType<:WeilDivisor}
  D::DivisorType
  f::Vector{<:VarietyFunctionFieldElem}

  function LinearSystem(f::Vector, D::WeilDivisor; check::Bool=true)
    length(f) == 0 && return new{typeof(D)}(D, Vector{VarietyFunctionFieldElem}())
    KK = parent(f[1])
    all(g -> (parent(g) === KK), f[2:end]) || error("elements must have the same parent")
    X = scheme(D)
    X === variety(KK) || error("input not compatible")

    if check
      all(g->in_linear_system(g, D), f) || error("element not in linear system")
    end
    f = Vector{VarietyFunctionFieldElem}(f)
    return new{typeof(D)}(D, f)
  end
end
  
linear_system(f::Vector, D::WeilDivisor; check::Bool=true) = LinearSystem(f, D, check=check)

### essential getters 
@doc raw"""
    weil_divisor(L::LinearSystem)

Return the divisor `D` of the linear system `L = |D|`.
"""
function weil_divisor(L::LinearSystem) 
  return L.D
end
gens(L::LinearSystem) = L.f
ngens(L::LinearSystem) = length(L.f)
gen(L::LinearSystem,i::Int) = L.f[i]

@doc raw"""
    variety(L::LinearSystem)

Return the variety on which `L` is defined.
"""
variety(L::LinearSystem) = scheme(weil_divisor(L))
# an alias for the user's convenience 
scheme(L::LinearSystem) = variety(L)

@doc raw"""
    subsystem(L::LinearSystem, P::IdealSheaf, n::Int) -> LinearSystem

Given a linear system ``L = |D|``, a sheaf of prime ideals `P` 
and an integer `n`, return a pair ``(K, A)`` consisting
of the subsystem of elements in ``|D + P|`` and the representing 
matrix ``A`` for its inclusion into ``L`` on the given set 
of generators.
"""
function subsystem(L::LinearSystem, P::IdealSheaf, n::Int)
  # find one chart in which P is supported
  # TODO: There might be preferred choices for charts with 
  # the least complexity.
  X = variety(L)
  X === space(P) || error("input incompatible")
  C = default_covering(X)
  U = first(patches(C))
  for V in patches(C)
    if !(one(OO(V)) in P(V))
      U = V
      break
    end
  end
  # Now U it is.

  # Assemble the local representatives
  R = ambient_coordinate_ring(U)
  loc_rep = [g[U] for g in gens(L)]
  common_denominator = gcd([denominator(g) for g in loc_rep])
  numerators = [numerator(g)*divexact(common_denominator, denominator(g)) for g in loc_rep]
  RP, _ = Localization(R, complement_of_prime_ideal(saturated_ideal(P(U))))
  PP = RP(prime_ideal(inverted_set(RP)))
  denom_mult = (_minimal_power_such_that(PP, I -> !(RP(common_denominator) in I))[1])-1
  w = n + denom_mult # Adjust!
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
  A = zero_matrix(kk, ngens(L), length(all_mons))
  for i in 1:ngens(L)
    for (c, m) in zip(coefficients(images[i]), monomials(images[i]))
      k = findfirst(x->(x==m), all_mons)
      A[i, k] = c
    end
  end

  r, K = left_kernel(A)
  new_gens = [sum([K[i,j]*gen(L, j) for j in 1:ncols(K)]) for i in 1:nrows(K)]
  return LinearSystem(new_gens, weil_divisor(L) + n*WeilDivisor(P)), K
end

