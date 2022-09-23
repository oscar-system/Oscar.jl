export AffineGroupRepresentation

export affine_group_scheme, reynolds_operator, set_reynolds_operator!, vector_space_poly_ring

export canonical_representation, induced_representation_on_symmetric_power

export omega_process, nullcone_ideal, reynolds_operator_from_omega_process, check_invariance_function, invariant_ring, lift_to_invariant_polynomial_func

########################################################################
# 
# A representation of an affine group ``G ⊂ 𝔸ⁿ`` on a finite 
# dimensional vector space ``V``. 
#
# The information on the representation is specified and stored 
# via a matrix ``A ∈ 𝒪(G)ʳˣʳ`` describing the morphism 
# ``G → Aut(V)`` with ``r = dim V``; the so-called `coordinate_matrix`
# of the representation.
#
########################################################################
@attributes mutable struct AffineGroupRepresentation{BRT<:Ring, BRET<:RingElem, 
                                                     GroupType<:AbsAffineGroupScheme, 
                                                     MatrixType<:MatrixElem
                                                    }
  G::GroupType
  M::MatrixType
  S::MPolyRing # the polynomial ring for the vector space

  function AffineGroupRepresentation(
      G::AbsAffineGroupScheme, M::MatrixElem; 
      var_names::Vector{String}=["v$i" for i in 1:nrows(M)],
      check::Bool=true
    )
    R = base_ring(M)
    kk = ground_field(G)
    Mconv = map_entries(OO(G), M) # will complain if matrix is not compatible
    S, _ = PolynomialRing(kk, var_names)


    if check
      # check compatibility of group law with matrix multiplication, etc.
    end

    return new{typeof(kk), elem_type(kk), typeof(G), typeof(Mconv)}(G, Mconv, S)
  end
end

affine_group_scheme(rho::AffineGroupRepresentation) = rho.G
coordinate_matrix(rho::AffineGroupRepresentation) = rho.M
vector_space_poly_ring(rho::AffineGroupRepresentation) = rho.S

@attr Union{Hecke.Map, Nothing} function reynolds_operator(rho::AffineGroupRepresentation) 
  S = vector_space_poly_ring(rho)
  return MapFromFunc(reynolds_operator_from_omega_process(rho), S, S)
end

function set_reynolds_operator!(rho::AffineGroupScheme, R::Hecke.Map)
  set_attribute!(rho, :reynolds_operator, R)
  return rho
end

function (rho::AffineGroupRepresentation{<:Any, <:Any, <:AffineMatrixGroup, <:Any})(A::MatrixElem)
  a = coordinates(affine_group_scheme(rho), A)
  return map_entries(f->(evaluate(f, a)), G.M)
end

###
# For a given representation ``ρ : G → Aut(V)`` this returns a 
# function `check(f::MPolyElem)` that takes a polynomial ``f`` on 
# the vector space ``V`` and returns `true` if ``f`` is ``ρ``-invariant
# and `false` otherwise.

@attr function check_invariance_function(rho::AffineGroupRepresentation)
  S = vector_space_poly_ring(rho)
  A = coordinate_matrix(rho)
  P = base_ring(A)
  PS, V = PolynomialRing(P, symbols(S))
  StoPS = hom(S, PS, V)
  APS = map_entries(x->PS(x), A)
  W = MatrixSpace(PS, 1, length(V))(V)*APS
  w = [W[1, i] for i in 1:length(W)]

  function func(f::MPolyElem)
    parent(f) == S || error("polynomial does not belong to the correct ring")
    F = StoPS(f)
    result = (F == evaluate(F, w))
    return result
  end

  return func
end

### A constructor for the canonical representation of a matrix group
canonical_representation(G::AffineMatrixGroup) = AffineGroupRepresentation(G, coordinate_matrix(G))

### 
# For a monomial x this returns the index i for 
# which it appears in the iteration over all monomials 
# in the base ring of the same degree as x.
function _linear_index_in_all_monomials(x::MPolyElem)
  d = total_degree(x)
  R = parent(x)
  length(coefficients(x)) == 1 || error("polynomial must be monomial")
  isone(first(coefficients(x))) || error("polynomial must be monomial")
  mons = Oscar.all_monomials(R, d)
  for (i, y) in zip(1:length(mons), mons)
    x == y && return i
  end
  error("monomial could not be found")
end

### Note: It would be nice if the user was given more control on the order 
# in which the monomials appear in the symmetric power of the original vector 
# space. By now, this is determined by the iterator Oscar.all_monomials, 
# but the whole construction should eventually be replaced by calls to 
# a function like `symmetric_power(V::VectorSpace)` and related mappings.
@Markdown.doc """
    induced_representation_on_symmetric_power(G::AffineMatrixGroup, p::Int)

For a matrix group ``G ⊂ 𝕜ᵐˣᵐ`` defined on a vector space ``V ≅ 𝕜ᵐ`` 
and an integer ``p`` this returns the induced representation on 
``Symᵖ Vᵛ``, the ``p``-th symmetric power of the dual of ``V``.
"""
function induced_representation_on_symmetric_power(G::AffineMatrixGroup, p::Int)
  M = coordinate_matrix(G)
  m = ncols(M)
  kk = ground_field(G)
  S, v = PolynomialRing(kk, ["v$i" for i in 1:m])
  mons = Oscar.all_monomials(S, p)
  n = length(mons)
  
  R = base_ring(M)
  RS, _ = PolynomialRing(R, symbols(S))
  StoRS = hom(S, RS, gens(RS))

  MRS = map_entries(RS, M)
  V = MatrixSpace(RS, 1, m)([StoRS(v) for v in gens(S)])
  C = V*MRS
  D = [C[1, i] for i in 1:ncols(C)]

  A = zero(MatrixSpace(R, n, n)) # the matrix for the induced representation
  for y in mons
    p = evaluate(y, D)
    for (c, z) in zip(coefficients(p), monomials(p))
      A[_linear_index_in_all_monomials(y), _linear_index_in_all_monomials(z)] = c
    end
  end

  return AffineGroupRepresentation(G, A)
end

### 
# This function takes a multivariate polynomial ∑ₐ cₐ ⋅ xᵃ 
# and turns it into the differential operator 
#
#   D = ∑ₐ cₐ ⋅ (∂/∂ x₁)ᵃ¹ ⋯ (∂/∂ xₙ)ᵃⁿ
#
# where a = (a₁,…,aₙ) is a multiindex.
function as_constant_differential_operator(f::MPolyElem)
  R = parent(f)
  function p(g::MPolyElem)
    parent(g) == R || error("element belongs to the wrong ring")
    is_zero(g) && return g
    result = zero(g)
    for (c, e) in zip(coefficients(f), exponent_vectors(f))
      h = copy(g)
      for i in 1:ngens(R)
        for j in 1:e[i]
          h = derivative(h, i)
          iszero(h) && continue
        end
        iszero(h) && continue
      end
      result = result + c*h
    end
    return result
  end
  return MapFromFunc(p, R, R)
end

function omega_process(G::AffineMatrixGroup)
  return as_constant_differential_operator(det(coordinate_matrix(G)))
end

@Markdown.doc """
    nullcone_ideal(rho::AffineGroupRepresentation)

For a representation `rho` of a linearly reductive affine matrix group 
``G`` on a vector space ``V`` this computes the ideal ``I ⊂ S`` 
in the polynomial ring ``S`` on ``V`` which describes the zero fiber 
of the projection to the categorical quotient ``V → V//G``.
"""
@attr MPolyIdeal function nullcone_ideal(
    rho::AffineGroupRepresentation{<:Any, <:Any, <:AffineGroupType, <:Any}
  ) where {
           AffineGroupType <: AffineMatrixGroup
          }

  # Prepare the coordinate matrix for Derksen's algorithm; 
  # it needs to be defined over an honest affine algebra.
  A = coordinate_matrix(rho)
  OOG = base_ring(A)
  R, I, f, phi, theta = as_affine_algebra(OOG)
  Q, pr = quo(R, I)
  OOGtoQ = hom(OOG, Q, gens(Q)[2:end], check=false)
  QtoOOG = hom(Q, OOG, vcat([inv(OOG(f))], gens(OOG)), check=false)
  A_lift = map_entries(x->lift(OOGtoQ(x)), A)

  # This code now follows [Derksen: Computation of Invariants 
  # for Reductive Groups, Advances in Mathematics 141, pp. 366--384,
  # 1999], Section 5.
  m = ncols(A)
  kk = coefficient_ring(R)
  Rext, _ = PolynomialRing(kk, vcat(symbols(R), Symbol.(["u$i" for i in 1:m]), Symbol.(["v$i" for i in 1:m])))
  z = gens(Rext)[1:ngens(R)]
  x = gens(Rext)[ngens(R)+1:ngens(R)+m]
  y = gens(Rext)[ngens(R)+m+1:ngens(R)+2*m]
  RtoRext = hom(R, Rext, z)
  M = map_entries(RtoRext, A_lift)
  X = MatrixSpace(Rext, 1, m)(x)
  Y = MatrixSpace(Rext, 1, m)(y)
  J = ideal(Rext, minors(Y - X*M, 1))
  J = J + ideal(Rext, RtoRext.(gens(I)))

  K = eliminate(J, z)
  
  S = vector_space_poly_ring(rho)
  v_ext = vcat([zero(S) for i in 1:ngens(R)], gens(S), [zero(S) for i in 1:m])
  g = (x->(evaluate(x, v_ext))).(gens(K))
  rop = reynolds_operator(rho)
  g = rop.(g)
  return ideal(S, g)
end

function reynolds_operator_from_omega_process(
    rho::AffineGroupRepresentation{<:Any, <:Any, <:AffineGroupType, <:Any}
  ) where {
           AffineGroupType <: AffineMatrixGroup
          }
  G = affine_group_scheme(rho)
  A = coordinate_matrix(rho)
  M = coordinate_matrix(G)
  all(x->(isone(lifted_denominator(x))), A) || error("reynolds operator from omega process not implemented for non-polynomial representations; please make sure that the coordinate matrix for your representation has trivial denominators")
  A_lift = map_entries(x->(lifted_numerator(x)), A)

  R = base_ring(M)
  m = ncols(M)
  kk = coefficient_ring(R)
  n = ncols(A_lift)
  S = vector_space_poly_ring(rho)
  RS, _ = PolynomialRing(kk, vcat(["t_$(i)_$(j)" for i in 1:m for j in 1:m], ["v$i" for i in 1:n]))
  t = gens(RS)[1:m^2]
  v = gens(RS)[m^2+1:end]
  RtoRS = hom(R, RS, t)
  T = map_entries(RtoRS, A_lift)
  V = MatrixSpace(RS, 1, n)(v)
  W = V*T
  w = [W[1, i] for i in 1:ncols(W)]
  det_T = det(T)
  det_M = RtoRS(det(M))
  Omega_t = as_constant_differential_operator(det_M)

  StoRS = hom(S, RS, v)

# function I_func(p::Int, q::Int, f::MPolyElem)
#   parent(f) == S || error("polynomial does not belong to the correct ring")
#   p < q && return zero(RS)
#   h = det_M^q*evaluate(f, [W[1, i] for i in 1:ncols(W)])
#   for i in 1:p
#     h = Omega_t(h)
#   end
#   return h
# end
#
# function I_func(g::Int, f::MPolyElem) 
#   return evaluate(I_func(g, 0, f), vcat([zero(S) for i in 1:ngens(R)], gens(S)))
# end
#
# o = lex(t)*lex(v)
#
  function reynolds(f::MPolyElem; raw::Bool=false)
    parent(f) == S || error("polynomial does not belong to the correct ring")
    is_zero(f) && return f
    h = evaluate(f, w)
    x_ext = vcat([zero(S) for i in 1:ngens(R)], gens(S))
    g1 = evaluate(h, x_ext)
    while iszero(g1)
      h = Omega_t(h)
      iszero(h) && return zero(f)
      g1 = evaluate(h, x_ext)
    end
    # Find out about the factor c
    # The raw version of the operator is c*projection.
    # In order to know c, we apply the raw version twice 
    # and compare.
    raw && return g1
    g2 = reynolds(g1, raw=true)
    check, c = divides(g2, g1)
    check || error("omega process does not give a projection")
    result =  div(g1, c)
    return result
  end

  return reynolds
end

@attr Tuple{<:MPolyQuo, <:Hecke.Map} function invariant_ring(rho::AffineGroupRepresentation)
  S = vector_space_poly_ring(rho)
  I = nullcone_ideal(rho)
  R = reynolds_operator(rho)
  inv = [R(g) for g in gens(I)]

  #check = check_invariance_function(rho)
  #@show [check(g) for g in inv]

  kk = coefficient_ring(S)
  A, y = PolynomialRing(kk, ["y$i" for i in 1:length(inv)])
  f = hom(A, S, inv)
  K = kernel(f)
  Q, p = quo(A, K)
  return Q, hom(Q, S, inv, check=false)
end

@attr function lift_to_invariant_polynomial_func(rho::AffineGroupRepresentation)
  S = vector_space_poly_ring(rho)
  A, pr = invariant_ring(rho)
  I = nullcone_ideal(rho)
  rop = reynolds_operator(rho)
  check_func = check_invariance_function(rho)

  Sgr, _ = grade(S)

  function is_homog(f::MPolyElem)
    return is_homogeneous(Sgr(f))
  end

  function homog_comp(f::MPolyElem)
    return [a.f for a in collect(values(homogeneous_components(Sgr(f))))]
  end

  function func(f::MPolyElem; check::Bool=false)
    parent(f) == S || error("polynomial does not belong to the correct ring")
    is_zero(f) && return zero(A)
    if !is_homog(f) 
      return sum([func(a) for a in homog_comp(f)])
    end
    check && (check_func(f) || error("the given polynomial is not an invariant"))

    is_constant(f) && return A(first(coefficients(f)))
    c = coordinates(f, I)
    c = rop.(c)
    rec_res = [func(a, check=false) for a in c]
    res = sum([A[i]*rec_res[i] for i in 1:length(c)])
    return res
  end

  return MapFromFunc(func, S, A)
end

