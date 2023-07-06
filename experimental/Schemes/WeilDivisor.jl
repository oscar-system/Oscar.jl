export LinearSystem
export WeilDivisor
export coefficient_dict
export coefficient_ring
export coefficient_ring_type
export coefficient_type
export components
export divisor
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
    @check begin
      for D in keys(coefficients)
        isprime(D) || error("components of a divisor must be sheaves of prime ideals")
        dim(X) - dim(D) == 1 || error("components of a divisor must be of codimension one")
      end
    end
    return new{typeof(X), CoefficientRingType, CoefficientRingElemType}(AlgebraicCycle(X, R, coefficients, check=check))
  end

  function WeilDivisor(C::AlgebraicCycle; check::Bool=true)
    X = scheme(C)
    @check begin
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
  return WeilDivisor(X, R, D, check=false)
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
    WeilDivisor(I::IdealSheaf, R::Ring=ZZ)

Return the `WeilDivisor` ``D = 1 ⋅ V(I)`` with coefficients 
in ``R`` for a sheaf of prime ideals ``I``.
"""
function WeilDivisor(I::IdealSheaf, R::Ring=ZZ; check::Bool=true)
  D = WeilDivisor(space(I), R)
  @check isprime(I) "ideal sheaf must be prime"
  @check dim(space(I)) - dim(I) == 1 "components of a divisor must be of codimension one"
  coefficient_dict(D)[I] = one(R)
  return D
end

weil_divisor(I::IdealSheaf, R::Ring; check::Bool=true) = WeilDivisor(I, R, check=check)

weil_divisor(I::IdealSheaf; check::Bool=true) = WeilDivisor(I, check=check)

### copy constructor
function copy(D::WeilDivisor) 
  new_dict = IdDict{IdealSheaf, elem_type(coefficient_ring_type(D))}()
  for I in keys(coefficient_dict(D))
    new_dict[I] = D[I]
  end
  return WeilDivisor(scheme(D), coefficient_ring(D), new_dict, check=false)
end

function Base.show(io::IO, D::WeilDivisor)
  if has_name(D)
    print(io, name(D))
    return
  end
  if length(components(D)) == 0
    print(io, "the zero Weil divisor on $(scheme(D))")
    return
  end
  println(io, "Weil divisor on $(scheme(D)) given as the formal sum:")
  comp = ["$(D[I]) ⋅ $(I)" for I in components(D)]
  join(io, comp, " + ")
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

@doc raw"""
    intersect(D::WeilDivisor, E::WeilDivisor)

For two `WeilDivisor`s on a complete smooth surface the intersection number is defined 
as in Hartshorne's "Algebraic Geometry". This computes this intersection number.
"""
function intersect(D::WeilDivisor, E::WeilDivisor)
  X = scheme(D)
  @assert dim(X) == 2 "intersection of Weil divisors is only implemented for surfaces."
  X === scheme(E) || error("divisors do not live on the same scheme")
  R = coefficient_ring(D)
  R === coefficient_ring(E) || error("divisors do not have the same coefficient ring")
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
  # TODO: Work out the intersection
  result = zero(R)
  for c1 in components(D)
    a1 = D[c1]
    for c2 in components(E)
      a2 = E[c2]
      I = c1 + c2
      @assert dim(I) <= 0 "divisors have nontrivial self intersection"
      result = result + a1 * a2 * colength(I)
    end
  end
  return result
end

function colength(I::IdealSheaf; covering::Covering=default_covering(scheme(I)))
  X = scheme(I)
  patches_todo = copy(patches(covering))
  patches_done = AbsSpec[]
  result = 0
  while length(patches_todo) != 0
    U = pop!(patches_todo)
    J = I(U)
    if has_decomposition_info(covering)
      h = decomposition_info(covering)[U]
      # The elements in h indicate where components must 
      # be located so that they can not be spotted in other charts.
      # We iteratively single out these components by adding a sufficiently high 
      # power of the equation to the ideal.
      for f in h
        g = f
        while !(g in ideal(OO(U), g*f) + J)
          g = g * g
        end
        J = J + ideal(OO(U), g)
        isone(J) && break
      end
    else
      # To avoid overcounting, throw away all components that 
      # were already visible in other charts.
      for V in patches_done
        if !haskey(glueings(covering), (U, V))
          continue
        end
        G = covering[U, V]
        (UV, VU) = glueing_domains(G)
        UV isa PrincipalOpenSubset || error("method is only implemented for simple glueings")
        f = complement_equation(UV)
        # Find a sufficiently high power of f such that it throws
        # away all components away from the horizon, but does not affect
        # those on the horizon itself.
        g = f
        while !(g in ideal(OO(U), g*f) + J)
          g = g * g
        end
        J = J + ideal(OO(U), g)
        isone(J) && break
      end
    end
    if !isone(J)
      JJ = leading_ideal(saturated_ideal(J))
      A, _ = quo(base_ring(JJ), JJ)
      result = result + ngens(vector_space(coefficient_ring(base_ring(A)), A)[1])
    end
    push!(patches_done, U)
  end
  return result
end


@doc raw"""
    in_linear_system(f::VarietyFunctionFieldElem, D::WeilDivisor; regular_on_complement::Bool=true) -> Bool

Check if the rational function `f` is in the linear system ``|D|``.
"""
function in_linear_system(f::VarietyFunctionFieldElem, D::WeilDivisor; regular_on_complement::Bool=false)
  X = scheme(D) 
  X === variety(parent(f)) || error("schemes not compatible")
  C = simplified_covering(X)
  for I in components(D)
    # no check needed because the components of a prime divisor are prime anyways
    order_on_divisor(f, I, check=false) >= -D[I] || return false
  end
  regular_on_complement && return true
  for U in patches(C)
    # we have to check that f[U] has no poles outside the support of D[U]
    J = intersect([J(U) for J in components(D)])
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

@doc raw"""
    linear_system(f::Vector, D::WeilDivisor; check::Bool=true)

Return the linear system ``L`` generated by rational functions ``f₁,…,fᵣ ∈ K(X)``
with $L \subseteq |D|$ for `D` on a variety `X`. If `check` is set,
confirm that $L \subseteq |D|$.
"""
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
    subsystem(L::LinearSystem, D::WeilDivisor) -> LinearSystem, MatElem

Given a linear system $L = |E|$ and a divisor $D \leq E$ compute $|D|$
and the matrix representing the inclusion $|D| \hookrightarrow |E|$
with respect to the given bases of both systems.
"""
function subsystem(L::LinearSystem, D::WeilDivisor)
  E = weil_divisor(L)
  @req D <= E "input does not define a subsystem"
  Lnew = L
  T = identity_matrix(base_ring(scheme(L)), ngens(L))
  for P in components(E)
    Lnew, Tnew = _subsystem(Lnew, P, -coeff(D,P))
    T = Tnew*T
  end
  return Lnew, T
end


@doc raw"""
    _subsystem(L::LinearSystem, P::IdealSheaf, n) -> LinearSystem

Given a linear system ``L = |D|``, a sheaf of prime ideals `P` 
and an integer `n`, return a pair ``(K, A)`` consisting
of the subsystem of elements in ``|D|`` that vanish to order at least n at ``P``.
The matrix ``A`` for its inclusion into ``L`` on the given set
of generators.
"""
function _subsystem(L::LinearSystem, P::IdealSheaf, n)
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
  if length(loc_rep) == 0
    common_denominator = R(1)
  else
    common_denominator = lcm([denominator(g) for g in loc_rep])
  end
  numerators = [numerator(g)*divexact(common_denominator, denominator(g)) for g in loc_rep]

  # compute a symbolic power
  RP, _ = Localization(OO(U), complement_of_prime_ideal(saturated_ideal(P(U))))
  PP = RP(prime_ideal(inverted_set(RP)))
  K = function_field(X)
  denom_mult = order_on_divisor(K(common_denominator), P, check=false)
  #denom_mult = (_minimal_power_such_that(PP, I -> !(RP(common_denominator) in I))[1])-1
  w = n + denom_mult # Adjust!
  if w < 0
    pPw = ideal(R, one(R))
  else
    pPw = saturated_ideal(PP^w) # the symbolic power
  end
  # reduce the numerators modulo P^(w)
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
  new_gens = [sum([K[i,j]*gen(L, j) for j in 1:ncols(K)]) for i in 1:r]
  return LinearSystem(new_gens, weil_divisor(L) + n*WeilDivisor(P), check=false), K
end

@attr Bool function is_prime(D::WeilDivisor)
  if length(components(D))!=1
    return false
  end
  c = components(D)[1]
  return coefficient_dict(D)[c] == 1
end

is_irreducible(D::WeilDivisor) = is_prime(D)

function order_on_divisor(
    f::VarietyFunctionFieldElem,
    D::WeilDivisor;
    check::Bool=true
  )
  @check is_prime(D) || error("divisor must be prime")
  I = components(D)[1]
  return order_on_divisor(f, I, check=false)
end

