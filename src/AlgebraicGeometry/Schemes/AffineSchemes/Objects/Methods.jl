

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
  elseif is_terse(io)
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

Check if a ring element is a non-zero-divisor
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

@doc raw"""
    is_normal(X::AbsAffineScheme; check::Bool=true) -> Bool

# Input:
- a reduced scheme ``X``,
- if `check` is `true`, then confirm that ``X`` is reduced; this is expensive.

# Output:
Return whether the scheme ``X`` is normal.

# Examples
```jldoctest
julia> R, (x, y, z) = QQ[:x, :y, :z];

julia> X = spec(R);

julia> is_normal(X)
true
```
"""
function is_normal(X::AbsAffineScheme; check::Bool=true)
  !check || is_reduced(X) || return false
  R = coordinate_ring(X)
  return is_normal(R; check=check)
end

# Todo: Normalizations of localizations
#function _normalization(X::AbsAffineScheme{<:Field, <:MPolyLocRing}; algorithm=:equidimDec)
#  L = OO(X)
#  R = base_ring(L)
#  Y, phi, loc = _normalization(spec(R))[1]
#end

"""
    _normalization

Return a vector of triples
`(Xnorm_k, f_k, phi_k)`
where the normalization is the disjoint union of the `Xnorm_k`
`f_k: Xnorm_k -> X` is the restriction of the normalization map.
Let `X_k` be the image of `f_k`.
Then `phi_k` maps the coordinate ring of `Xnorm_k` to the
total field of fraction of `OO(X_k)`.
"""
_normalization(X::AbsAffineScheme; algorithm, check::Bool=true) = error("not implemented")



function _normalization(X::AbsAffineScheme{<:Field, <:MPolyRing}; algorithm=:equidimDec, check::Bool=true)
  R = OO(X)
  K = total_ring_of_fractions(R)
  Y = spec(R)
  phi = morphism(Y, X, gens(OO(X)); check=false)
  loc = hom(R, K, [K(x, one(R), false) for x in gens(R)]; check=false)
  return [(X, phi, loc)]
end

function _normalization(X::AbsAffineScheme{<:Field, <:MPolyQuoRing}; algorithm=:equidimDec, check::Bool=true)
  @check is_reduced(X) "The scheme X=$(X) needs to be reduced."
  A = OO(X)
  A_norm = normalization(A; algorithm)
  output = Tuple{<:AbsAffineScheme, <:AbsAffineSchemeMor, <:Map}[]
  F = ambient_embedding(X)
  if length(A_norm) == 1 && is_one(A_norm[1][3][1])
    # Workaround, normalization for rings buggy if already normal
    Xnorm_k = X
    X_k = X
    F_k = identity_map(X)
    K_k = total_ring_of_fractions(X_k)
    A_k = OO(Xnorm_k)
    A_k_to_K_k = hom(A_k, K_k, K_k.(gens(A_k)); check=false)
    push!(output, (Xnorm_k, F_k, A_k_to_K_k))
    return output
  end
  for (A_k, f_k, a_k) in A_norm
    d_k = a_k[1] # An element so that 1/(d_k) J_k = A_k; J_k = d_k[2]
    Xnorm_k = spec(A_k)
    F_k = morphism(Xnorm_k, X, f_k; check=false)
    if length(A_norm) == 1
      # X = X_k is integral
      X_k = X
      inc_k = identity_map(X)
    else
      I = kernel(pullback(F_k))
      X_k, inc_k = sub(X, I) # the component of X lying over Xnorm_k
    end
    K_k = total_ring_of_fractions(X_k)
    d_k = pullback(inc_k)(d_k)
    gens_to_K_k_coordinates = vcat(
      [ K_k(pullback(inc_k)(g_k), d_k, false) for g_k in gens(a_k[2])[1:end-1] ],
      gens(OO(X_k))
    )
    A_k_to_K_k = hom(A_k, K_k, gens_to_K_k_coordinates; check=false)

    push!(output, (Xnorm_k, F_k, A_k_to_K_k))
  end
  return output
end

# further documented for `Scheme`
"""
    normalization(X::AbsAffineScheme) -> Vector{Tuple{AbsAffineScheme, AbsAffineSchemeMor}}

Return the normalization of the reduced affine scheme ``X``.

# Input:
- A reduced affine scheme ``X``
- if `check` is `true` confirm that ``X`` is reduced; this is expensive
- the keyword argument `algorithm` is passed on to [`normalization(::MPolyQuoRing)`](@ref)

# Output:
A list of pairs ``(Y_i, f_i)`` where ``Y_i`` is a normal scheme and
``f_i`` is a morphism from ``Y_i`` to ``X``.
The disjoint union of the ``Y_i`` is the normalization of ``X``
and the ``f_i`` are the restrictions of the normalization morphism to ``Y_i``.
"""
function normalization(X::AbsAffineScheme; check::Bool=true, algorithm=:equidimDec)
  @check is_reduced(X) "only reduced schemes can be normalized"
  norm_outputs = _normalization(X; algorithm)
  return [i[1:2] for i in norm_outputs]
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
  return Y, morphism(Y, X, Phi; check=false)
end

### Some helper functions
function _change_base_ring(phi::Any, R::MPolyRing)
  K = coefficient_ring(R)
  kk = parent(phi(zero(K)))
  P, _ = polynomial_ring(kk, symbols(R); cached=false)
  Phi = hom(R, P, phi, gens(P); check=false)
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
  @assert _has_coefficient_map(Phi)
  U = inverted_set(W)
  U_red = MPolyPowersOfElement(P, Phi.(denominators(U)))
  W_red, loc_map = localization(P, U_red)
  comp = hom(R, W_red, phi, gens(W_red); check=false)
  @assert _has_coefficient_map(comp)
  res_map = hom(W, W_red, comp, check=false)
  @assert _has_coefficient_map(res_map)
  return W_red, res_map
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
  @assert _has_coefficient_map(res)
  res_map = hom(L, L_red, res, check=false)
  @assert _has_coefficient_map(res_map)
  return L_red, res_map
end

function Base.hash(X::Scheme, u::UInt)
  return u
end

