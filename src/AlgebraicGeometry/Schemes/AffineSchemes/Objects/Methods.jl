

####################################################################################
# (1) Equality
####################################################################################

# TODO: Add further cross-type comparison methods as needed.

function ==(X::T, Y::T) where {T<:AbsAffineScheme{<:Ring, <:MPolyRing}}
  return OO(X) === OO(Y)
end


function ==(X::AbsAffineScheme, Y::AbsAffineScheme)
  X === Y && return true
  return is_subscheme(X, Y) && is_subscheme(Y, X)
end


function ==(X::AbsAffineScheme, Y::EmptyScheme)
  return is_subscheme(X, Y)
end


==(X::EmptyScheme, Y::AbsAffineScheme) = (Y == X)


########################################################
# (2) Display
########################################################

# We show a detailed version of the coordinate ring since they are all the
# details we can get.. Otherwise our detailed printing is quite poor and
# "useless".
function Base.show(io::IO, ::MIME"text/plain", X::AbsAffineScheme)
  io = pretty(io)
  println(io, "Spectrum")
  print(io, Indent(), "of ", Lowercase())
  show(io, MIME"text/plain"(), OO(X))
  print(io, Dedent())
end

function Base.show(io::IO, X::AbsAffineScheme)
  if has_attribute(X, :name)
    print(io, name(X))
  elseif get(io, :supercompact, false)
    print(io, "Affine scheme")
  elseif get_attribute(X, :is_empty, false)
    print(io, "Empty affine scheme")
  else
    _show(io, X)
  end
end

function _show(io::IO, X::AbsAffineScheme)
  io = pretty(io)
  print(io, LowercaseOff(), "Spec of ")
  print(io, Lowercase(), OO(X))
end

function _show(io::IO, X::AbsAffineScheme{<:Any,<:MPolyRing})
  print(io, "Affine ", ngens(OO(X)), "-space")
end

function _show(io::IO, X::AbsAffineScheme{<:Any,<:MPolyQuoRing})
  io = pretty(io)
  print(io, "scheme(")
  I = modulus(OO(X))
  join(io, gens(I), ", ")
  print(io, ")")
end

function _show(io::IO, X::AbsAffineScheme{<:Any, <:MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}})
  io = pretty(io)
  print(io, "scheme(")
  I = modulus(OO(X))
  S = inverted_set(OO(X))
  join(io, gens(I), ", ")
  if is_empty(denominators(S))
    print(io, raw")")
  else
    AA = AbstractAlgebra
    PP = AbstractAlgebra.PrettyPrinting
    print(io, raw") \ ")
    AA.show_obj(io,
                MIME("text/plain"),
                PP.canonicalize(Expr(:call,
                                     :scheme,
                                     Expr(:call, :*, PP.expressify.(denominators(S))...)
                                    )
                               )
               )
  end
end

function _show(io::IO, X::AbsAffineScheme{<:Any, <:MPolyLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}})
  io = pretty(io)
  print(io, LowercaseOff(), "AA^", ngens(OO(X)))
  S = inverted_set(OO(X))
  if is_empty(denominators(S))
    # do nothing more
  else
    AA = AbstractAlgebra
    PP = AbstractAlgebra.PrettyPrinting
    print(io, raw" \ ")
    AA.show_obj(io,
                MIME("text/plain"),
                PP.canonicalize(Expr(:call,
                                     :scheme,
                                     Expr(:call, :*, PP.expressify.(denominators(S))...)
                                    )
                               )
               )
  end
end

########################################################
# (3) Check for zero divisors in rings
########################################################

@doc raw"""
    is_non_zero_divisor(f::RingElem, X::AbsAffineScheme)

Checks if a ring element is a non-zero divisor
in the coordinate ring of an affine scheme.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Affine space of dimension 3
  over rational field
with coordinates [x1, x2, x3]

julia> (x1, x2, x3) = gens(OO(X))
3-element Vector{QQMPolyRingElem}:
 x1
 x2
 x3

julia> is_non_zero_divisor(x1, X)
true

julia> is_non_zero_divisor(zero(OO(X)), X)
false
```
"""
function is_non_zero_divisor(f::RingElem, X::AbsAffineScheme)
  error("method not implemented for affine schemes of type $(typeof(X))")
end

function is_non_zero_divisor(f::RingElem, X::AbsAffineScheme{<:Ring, <:MPolyRing})
  return !iszero(OO(X)(f))
end

function is_non_zero_divisor(f::MPolyQuoRingElem, X::AbsAffineScheme{<:Ring, <:MPolyQuoRing})
  R = ambient_coordinate_ring(X)
  I = modulus(OO(X))
  J = ideal(R, lift(f))
  return I == quotient(I, J)
end

function is_non_zero_divisor(f::RingElem, X::AbsAffineScheme{<:Ring, <:MPolyLocRing})
  return !iszero(OO(X)(f))
end

function is_non_zero_divisor(f::RingElem, X::AbsAffineScheme{<:Ring, <:MPolyQuoLocRing})
  I = ideal(OO(X), [zero(OO(X))])
  zero_ideal = Oscar.pre_image_ideal(I)
  J = Oscar.pre_image_ideal(ideal(OO(X), [f]))
  Q = quotient(zero_ideal, J)
  return zero_ideal == Q
end

########################################################################
# Normalization                                                        #
########################################################################

## Example
# R, (x, y, z) = grade(QQ["x", "y", "z"][1])
# I = ideal(R, z*x^2 + y^3)
# X = covered_scheme(proj(R, I))
# C = X.coverings[1]
# gluing_example = X[1][1,2]
# Oscar.normalization(C)
# Oscar.normalization(X)

# Warning: assume scheme integral
function normalization(X::AbsCoveredScheme; check::Bool=true)
  orig_cov = default_covering(X)
  if has_attribute(X, :simplified_covering)
    orig_cov = simplified_covering(X)
  end
  X_norm_cov, pr_cov = normalization(orig_cov; check)
  Y = CoveredScheme(X_norm_cov)
  pr = CoveredSchemeMorphism(Y, X, pr_cov)
  return Y, pr
end

# Warning: assume scheme integral
function normalization(C::Covering; check::Bool=true)
  X_i_norm_outputs = [_normalization(U; check) for U in patches(C)]
  C_norm = Covering(first.(first.(X_i_norm_outputs)))
  n = n_patches(C_norm)
  for i in 1:n
    for j in (i+1):n
      gluing_i_j = _normalization(C[i,j], X_i_norm_outputs[i][1], X_i_norm_outputs[j][1])
      add_gluing!(C_norm, gluing_i_j)
    end
  end
  second(xs) = xs[2]
  morphisms = second.(first.(X_i_norm_outputs))
  covering_morphism_dict = IdDict(zip(domain.(morphisms), morphisms))
  normalization_map = CoveringMorphism(C_norm, C, covering_morphism_dict)
  return C_norm, normalization_map
end

struct NormalizationGluingData
  original_gluing::AbsGluing
  norm_output1 # TODO: Add some types here?
  norm_output2
  check::Bool
end

function _normalization(
    G::AbsGluing,
    X_1_norm_output, # output of _normalization
    X_2_norm_output; # output of _normalization
    check::Bool=true
  )
  data = NormalizationGluingData(G, X_1_norm_output, X_2_norm_output, check)
  X, Y = patches(G)
  return LazyGluing(X, Y, _compute_normalization, data)
end

# Warning: assume patches irreducible
# Input is gluing of X_1 and X_2 and the output of the _normalization
# function of X_1 and X_2
function _compute_normalization(
    data::NormalizationGluingData
  )
  # Initialize the variables
  G = data.original_gluing
  X_1_norm_output = data.norm_output1
  X_2_norm_output = data.norm_output2
  check = data.check

  A, B = patches(G)
  @check is_integral(A) && is_integral(B) "schemes must be integral"
  (X_1_norm, F_1, hom_to_K_1) = X_1_norm_output
  (X_2_norm, F_2, hom_to_K_2) = X_2_norm_output
  X_1, X_2 = patches(G)
  U_1, U_2 = gluing_domains(G)
  g_1, g_2 = gluing_morphisms(G)
  codomain(F_1) === X_1 || error("Codomain of F_1=$(F_1) should be X_1=$(X_1).")
  codomain(F_2) === X_2 || error("Codomain of F_2=$(F_2) should be X_2=$(X_2).")

  # We assume X_i_norm, F_i, hom_to_K_i are defined as below
  # X_1_norm, F_1, hom_to_K_1 = _normalization(X_1)[1]
  # X_2_norm, F_2, hom_to_K_2 = _normalization(X_2)[1]

  U_1_norm = PrincipalOpenSubset(X_1_norm, pullback(F_1)(complement_equation(U_1)))
  U_2_norm = PrincipalOpenSubset(X_2_norm, pullback(F_2)(complement_equation(U_2)))

  # TODO: disable check
  U_1_norm_mor = restrict(F_1, U_1_norm, U_1)
  U_2_norm_mor = restrict(F_2, U_2_norm, U_2)

  hom_to_U_2_norm_coordinates = elem_type(OO(U_1_norm))[]
  for a in hom_to_K_2.(coordinates(X_2_norm))
    b_num = pullback(U_1_norm_mor)(pullback(g_1)(OO(U_2)(numerator(a))))
    b_denom = pullback(U_1_norm_mor)(pullback(g_1)(OO(U_2)(denominator(a))))
    # TODO remove check that b_denom is invertible
    push!(hom_to_U_2_norm_coordinates, b_num * inv(b_denom))
  end
  hom_to_U_2_norm = morphism(U_1_norm, U_2_norm, hom_to_U_2_norm_coordinates)

  hom_to_U_1_norm_coordinates = elem_type(OO(U_2_norm))[]
  for a in hom_to_K_1.(coordinates(X_1_norm))
    b_num = pullback(U_2_norm_mor)(pullback(g_2)(OO(U_1)(numerator(a))))
    b_denom = pullback(U_2_norm_mor)(pullback(g_2)(OO(U_1)(denominator(a))))
    # TODO remove check that b_denom is invertible
    push!(hom_to_U_1_norm_coordinates, b_num * inv(b_denom))
  end
  hom_to_U_1_norm = morphism(U_2_norm, U_1_norm, hom_to_U_1_norm_coordinates)

  return SimpleGluing(
    X_1_norm,
    X_2_norm,
    hom_to_U_2_norm,
    hom_to_U_1_norm;
    check=true
  )
end

# Todo: Normalizations of localizations
#function _normalization(X::AbsAffineScheme{<:Field, <:MPolyLocRing}; algorithm=:equidimDec)
#  L = OO(X)
#  R = base_ring(L)
#  Y, phi, loc = _normalization(spec(R))[1]
#end

function _normalization(X::AbsAffineScheme{<:Field, <:MPolyRing}; algorithm=:equidimDec, check::Bool=true)
  R = OO(X)
  K = total_ring_of_fractions(R)
  Y = spec(R)
  phi = morphism(Y, X, gens(OO(X)); check=false)
  loc = hom(R, K, [K(x, one(R), true) for x in gens(R)])
  return [(X, phi, loc)]
end

# Warning: assume scheme reduced
function _normalization(X::AbsAffineScheme{<:Field, <:MPolyQuoRing}; algorithm=:equidimDec, check::Bool=true)
  @check is_integral(X)
  A = OO(X)
  A_norm = normalization(A; algorithm)
  output = Tuple{<:AbsAffineScheme, <:AbsAffineSchemeMor, <:Map}[]
  for (A_k, f_k, a_k) in A_norm
    K = total_ring_of_fractions(X)
    d_k = base_ring(K)(a_k[1]) # An element so that 1/(d_k) J_k = A_k; J_k = d_k[2]
    # Workaround, normalization for rings buggy if already normal
    if is_one(d_k)
      X_k = X
      F_k = identity_map(X)
      homomorphism_to_K = hom(A_k, K, coordinates(X_k); check=true)
      push!(output, (X_k, F_k, homomorphism_to_K))
      continue
    end
    X_k = spec(A_k)
    F_k = morphism(X_k, X, f_k; check=true) # Set to false after testing
    homomorphism_to_K_coordinates = vcat(
      [ K(base_ring(K)(g_k), d_k, true) for g_k in gens(a_k[2])[1:end-1] ],
      gens(A),
    )
    homomorphism_to_K = hom(A_k, K, homomorphism_to_K_coordinates; check=true) # Set to false after testing
    push!(output, (X_k, F_k, homomorphism_to_K))
  end
  return output
end

########################################################################
# High level constructors of subschemes                                #
########################################################################

function sub(X::AbsAffineScheme, I::Ideal)
  inc = ClosedEmbedding(X, I) # Will complain if input is not compatible
  return domain(inc), inc
end

########################################################################
# Further intersections                                                #
########################################################################

function intersect(U::PrincipalOpenSubset, V::PrincipalOpenSubset)
  if !(ambient_scheme(U) === ambient_scheme(V))
    return intersect(underlying_scheme(U), underlying_scheme(V))
  end
  return PrincipalOpenSubset(ambient_scheme(U), complement_equation(U)*complement_equation(V))
end


# TODO: The method below is necessary for a temporary hotfix; see #1882.
# It seems likely that the implementation can be tuned so that, for
# instance, the massive use of primary decompositions can be avoided.
# This should eventually be addressed.
@doc raw"""
    connected_components(X::AbsAffineScheme)

Return a decomposition of ``X`` into its connected components
``X = U₁ ∪ U₂ ∪ … ∪ Uₙ`` with ``Uᵢ`` a `PrincipalOpenSubset` of ``X``.
"""
function connected_components(X::AbsAffineScheme)
  I = saturated_ideal(modulus(OO(X)))
  l = primary_decomposition(I)
  comp = [subscheme(X, Q) for (Q, _) in l]
  found_intersection = true
  while found_intersection
    found_intersection = false
    for i in 1:length(comp)-1
      for j in i+1:length(comp)
        if !isempty(intersect(comp[i], comp[j]))
          found_intersection = true
          C2 = popat!(comp, j)
          C1 = comp[i]
          C_new = subscheme(X, intersect(saturated_ideal(modulus(OO(C1))),
                                         saturated_ideal(modulus(OO(C2)))
                                        ))
          comp[i] = C_new
          break
        end
      end
      found_intersection && break
    end
  end

  # We need to reproduce the components of X as `PrincipalOpenSubset`s.
  # To this end, we first determine separating equations between pairs
  # of components.

  n = length(comp)
  separating_equations = Array{RingElem}(undef, n, n)
  for i in 1:length(comp)
    C1 = comp[i]
    for j in i:length(comp)
      C2 = comp[j]
      if i == j
        separating_equations[i, j] = one(OO(X))
        continue
      end

      v = OO(X).(gens(saturated_ideal(modulus(OO(C1)))))
      w = OO(X).(gens(saturated_ideal(modulus(OO(C2)))))

      J = ideal(OO(X), vcat(v, w))
      # The idea is to write 1 = f + g with f ∈ I(C1) and g ∈ I(C2).
      # Then 1 - f = g vanishes identically on C2, but is a unit in OO(C1)
      # and vice versa.
      c = coordinates(one(OO(X)), J)

      # Conversion to assure compatibility.
      # This can be removed, once the ideal interface is streamlined.
      if c isa MatrixElem
        c = [c[1, i] for i in 1:ncols(c)]::Vector
      end

      a = c[1:length(v)]
      b = c[length(v)+1:length(v)+length(w)]

      f = sum([v[i]*a[i] for i in 1:length(v)], init=zero(OO(X)))
      g = sum([w[i]*b[i] for i in 1:length(w)], init=zero(OO(X)))
      separating_equations[i, j] = 1 - f
      separating_equations[j, i] = 1 - g
    end
  end
  result = Vector{PrincipalOpenSubset}()
  for i in 1:n
    f = prod(separating_equations[i, :], init=one(OO(X)))
    push!(result, PrincipalOpenSubset(X, f))
  end
  return result
end

########################################################################
# Base change and reduction modulo p
########################################################################

@doc raw"""
    base_change(phi::Any, X::AbsAffineScheme)

For an affine scheme `X` over a `base_ring` ``𝕜`` and a morphism
``φ : 𝕜 → 𝕂`` this computes ``Y = X × Spec(𝕂)`` and returns a pair
`(Y, psi)` where `psi` is the canonical map ``Y → X``.
"""
function base_change(phi::Any, X::AbsAffineScheme)
  kk = base_ring(X)
  kk_red = parent(phi(zero(kk)))
  R = OO(X)
  R_red, Phi = _change_base_ring(phi, R)
  Y = spec(R_red)
  return Y, morphism(Y, X, Phi)
end

### Some helper functions
function _change_base_ring(phi::Any, R::MPolyRing)
  K = coefficient_ring(R)
  kk = parent(phi(zero(K)))
  P, _ = polynomial_ring(kk, symbols(R))
  Phi = hom(R, P, phi, gens(P))
  return P, Phi
end

function _change_base_ring(phi::Any, A::MPolyQuoRing)
  R = base_ring(A)
  I = modulus(A)
  P, Phi = _change_base_ring(phi, R)
  I_red = ideal(P, Phi.(gens(I)))
  Q, pr = quo(P, I_red)
  Phi_bar = hom(A, Q, phi, gens(Q), check=false)
  return Q, Phi_bar
end

function _change_base_ring(phi::Any,
    W::MPolyLocRing{<:Any, <:Any, <:Any, <:Any,
                    <:MPolyPowersOfElement}
  )
  R = base_ring(W)
  P, Phi = _change_base_ring(phi, R)
  U = inverted_set(W)
  U_red = MPolyPowersOfElement(P, Phi.(denominators(U)))
  W_red, loc_map = localization(P, U_red)
  return W_red, hom(W, W_red, compose(Phi, loc_map), check=false)
end

function _change_base_ring(phi::Any,
    L::MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any,
                       <:MPolyPowersOfElement}
  )
  R = base_ring(L)
  W = localized_ring(L)
  W_red, Phi_W = _change_base_ring(phi, W)
  I = modulus(L)
  I_red = ideal(W_red, Phi_W.(gens(I)))
  L_red, pr = quo(W_red, I_red)
  res = compose(restricted_map(Phi_W), pr)
  return L_red, hom(L, L_red, res, check=false)
end

