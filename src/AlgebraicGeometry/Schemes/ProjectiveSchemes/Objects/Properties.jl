#######################################################################
# Basic properties
#######################################################################
@attr Bool function is_empty(P::AbsProjectiveScheme{<:Field, <:MPolyQuoRing})
  I = defining_ideal(P)
  R = base_ring(I)
  J = ideal(R, gens(R))
  return is_one(saturation(I, J))
end

@attr Bool function is_empty(P::AbsProjectiveScheme{<:Ring, <:MPolyDecRing})
  return isone(zero(base_ring(P)))
end

@attr Bool function is_empty(P::AbsProjectiveScheme{<:Ring, <:MPolyQuoRing})
  return all(x->radical_membership(x, defining_ideal(P)), gens(ambient_coordinate_ring(P)))
end

@doc raw"""
    is_smooth(P::AbsProjectiveScheme; algorithm::Symbol=:default) -> Bool

Check whether the scheme `P` is smooth.

# Algorithms

There are three possible algorithms for checking smoothness, determined
by the value of the keyword argument `algorithm`:
  * `:projective_jacobian` - uses the Jacobian criterion for projective
    schemes, see Exercise 4.2.10 of [Liu06](@cite),
  * `:covered_jacobian` - uses covered version of the Jacobian criterion,
  * `:affine_cone` - checks that the affine cone is smooth outside the origin.

The `:projective_jacobian` and the `:covered` algorithms only work for equidimensional schemes. The algorithms first check for equidimensionality, which can be expensive.
If you already know that the scheme is equidimensional, then you can
avoid recomputing that by writing
  `set_attribute!(P, :is_equidimensional, true)`
before checking for smoothness.

The algorithms `:covered_jacobian` and `:affine_cone` only work when the base ring is a field.

The default algorithm is `:projective_jacobian` if the scheme is equidimensional, otherwise it is `:affine_cone`.

# Examples
```jldoctest
julia> A, (x, y, z) = grade(QQ[:x, :y, :z][1]);

julia> B, _ = quo(A, ideal(A, [x^2 + y^2]));

julia> C = proj(B)
Projective scheme
  over rational field
defined by ideal (x^2 + y^2)

julia> is_smooth(C)
false

julia> is_smooth(C; algorithm=:covered_jacobian)
false
```
"""
is_smooth(P::AbsProjectiveScheme; algorithm::Symbol=:default)

@attr Bool function is_smooth(P::AbsProjectiveScheme{<:Ring, <:MPolyQuoRing}; algorithm::Symbol=:default)
  if is_empty(P)
    return true
  end

  if algorithm == :default
    if (has_attribute(P, :is_equidimensional) && is_equidimensional(P))
      algorithm = :projective_jacobian
    elseif base_ring(P) isa Field
      algorithm = :affine_cone
    else
      if !(base_ring(P) isa Field)
        throw(NotImplementedError(
          :is_smooth,
          "is_smooth only implemented when the scheme is equidimensional or the base ring is a field"
        ))
      end
    end
  end

  algorithms = [
    :projective_jacobian,
    :covered_jacobian,
    :affine_cone,
  ]
  if !(algorithm in algorithms)
    throw(ArgumentError(
      "the optional argument to the function is_smooth can only be one"
        * " of the following: " * join(algorithms, ", ") * "."
    ))
  end

  if algorithm == :covered_jacobian
    if !(base_ring(P) isa Field)
      throw(NotImplementedError(
        :is_smooth,
        "Algorithm `:covered_jacobian` only implemented when the base ring is a field"
        # because this algorithm uses `is_smooth` for affine schemes, and `is_smooth` is not implemented for affine schemes over a non-field base ring
      ))
    end
    return _jacobian_criterion(covered_scheme(P))
  elseif algorithm == :affine_cone
    # Compute the singular locus of the affine cone (an affine scheme)
    # and test whether the affine dimension of this locus is less than 1.
    # The dimension check avoids saturation of the result by the origin.
    # Compute_radical is switched off, as we only want the dimension.
    SL = _affine_cone_singular_locus_ideal(P; compute_radical=false, saturate=false)
    return dim(SL) < 1
  elseif algorithm == :projective_jacobian
    return is_one(_projective_jacobian_singular_locus_ideal(P))
  end
end

is_smooth(P::AbsProjectiveScheme{<:Ring, <:MPolyRing}; algorithm::Symbol=:default) = true
is_smooth(P::AbsProjectiveScheme{<:Ring, <:MPolyLocRing}; algorithm::Symbol=:default) = true

function _affine_cone_singular_locus_ideal(P::AbsProjectiveScheme{<:Ring, <:MPolyQuoRing}; saturate::Bool=true, compute_radical::Bool=false)
  if !(base_ring(P) isa Field)
    throw(NotImplementedError(
      :is_smooth,
      "affine Jacobian criterion only implemented when the base ring is a field"
      # because the underlying algorithms are not implemented for
      # affine schemes over a non-field base ring
    ))
  end
  # TODO: `non_smooth_locus` and `non_regular_locus` for affine schemes over 
  # non-field baserings. Then, this algorithm would work for arbitrary schemes. 
  # We explain why the algorithm of `:affine_cone` works for arbitrary schemes over arbitrary base schemes. By Remark 13.38(1) of [GW20](@cite), the morphism from the pointed affine cone to $P$ is locally the morphism $\mathbb{A}_U^1 \setminus \{0\} \to U$, where $U$ is an affine open of $P$ and $\mathbb{A}_U^1$ is the relative affine 1-space over $U$. By Definition 6.14(1) of [GW20](@cite), the Jacobian matrix for $U$ differs from the Jacobian matrix for $P$ only by a column containing zeros, implying that the ranks of the Jacobian matrices are the same. Therefore, $P$ is smooth if and only if the affine cone is smooth outside the origin.
  aff, _ = affine_cone(P)
  sing, _ = singular_locus(aff; compute_radical)
  !saturate && return saturated_ideal(defining_ideal(sing))
  origin = ideal(ambient_coordinate_ring(sing),gens(ambient_coordinate_ring(sing)))
  return saturation(saturated_ideal(defining_ideal(sing)), origin)
end

function _projective_jacobian_singular_locus_ideal(P::AbsProjectiveScheme{<:Ring, <:MPolyQuoRing}; saturate::Bool=true)
  if !(is_equidimensional(P))
    throw(NotImplementedError(
      :is_smooth,
      "projective Jacobian criterion only implemented when scheme is equidimensional"
    ))
  end
  R = ambient_coordinate_ring(P)
  I = defining_ideal(P)
  mat = jacobian_matrix(R, gens(I))
  sing_locus = ideal(R, minors(mat, codim(P))) + I
  !saturate && return sing_locus
  irrelevant_ideal = ideal(R, gens(R))
  return saturation(sing_locus, irrelevant_ideal)
end

@doc raw"""
    singular_locus(P::AbsProjectiveScheme; algorithm::Symbol=:default, saturate::Bool=true) -> AbsProjectiveScheme

Given an equidimensional projective scheme, returns the possibly
nonreduced subscheme of projective space describing the singular locus.

If the boolean keyword argument `saturate` is true (which is the default
value), then we saturate the ideal with respect to the irrelevant ideal.

# Algorithms

There are two possible algorithms for computing the singular locus, determined
by the value of the keyword argument `algorithm`:
  * `:projective_jacobian` - uses the Jacobian criterion for projective
    schemes, see Exercise 4.2.10 of [Liu06](@cite),
  * `:affine_cone` - uses the singular locus of the affine cone.

The `:projective_jacobian` algorithm only works for equidimensional schemes.
The algorithm first checks for equidimensionality, which can be expensive.
If you already know that the scheme is equidimensional, then you can
avoid recomputing that by writing
  `set_attribute!(P, :is_equidimensional, true)`

The algorithm `:affine_cone` only works when the base ring is a field.

The default algorithm is `:projective_jacobian` if the scheme is
equidimensional, otherwise it is `:affine_cone`.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> I = ideal(R, x^2 + y^2)
Ideal generated by
  x^2 + y^2

julia> X = proj(I)
Projective scheme
  over rational field
defined by ideal (x^2 + y^2)

julia> singular_locus(X)
Projective scheme
  over rational field
defined by ideal (y, x)
```
"""
singular_locus(P::AbsProjectiveScheme; algorithm::Symbol=:default, saturate::Bool=true)

function singular_locus(P::AbsProjectiveScheme{<:Ring, <:MPolyQuoRing}; algorithm::Symbol=:default, saturate::Bool=true)
  if is_empty(P)
    return P
  end

  if algorithm == :default
    if (has_attribute(P,:is_equidimensional) && is_equidimensional(P))
      algorithm = :projective_jacobian
    elseif base_ring(P) isa Field
      algorithm = :affine_cone
    else
      if !(base_ring(P) isa Field)
        throw(NotImplementedError(
          :singular_locus,
          "singular_locus only implemented when the scheme is equidimensional or the base ring is a field"
        ))
      end
    end
  end

  algorithms = [
    :projective_jacobian,
    :affine_cone,
  ]
  if !(algorithm in algorithms)
    throw(ArgumentError(
      "the optional argument to the function is_smooth can only be one"
        * " of the following: " * join(algorithms, ", ") * "."
    ))
  end

  if algorithm == :affine_cone
    I_aff = _affine_cone_singular_locus_ideal(P; saturate)
    R = ambient_coordinate_ring(P)
    I = ideal(R, map(R, gens(I_aff)))
    return subscheme(ambient_space(P), I)
  elseif algorithm == :projective_jacobian
    return subscheme(ambient_space(P), _projective_jacobian_singular_locus_ideal(P; saturate))
  end
end

singular_locus(P::AbsProjectiveScheme{<:Ring, <:MPolyRing}; algorithm::Symbol=:default, saturate::Bool=true) = subscheme(P, ideal(ambient_coordinate_ring(P), [1]))
singular_locus(P::AbsProjectiveScheme{<:Ring, <:MPolyLocRing}; algorithm::Symbol=:default, saturate::Bool=true) = subscheme(P, ideal(ambient_coordinate_ring(P), [1]))

@attr Bool function is_irreducible(P::AbsProjectiveScheme)
  is_empty(P) && return false  # the empty set is not irreducible
  I = defining_ideal(P)
  return is_prime(radical(I)) # topological property
end

@attr Bool function is_reduced(P::AbsProjectiveScheme)
  I = defining_ideal(P)
  return is_radical(I)
end

@attr Bool function is_geometrically_reduced(P::AbsProjectiveScheme{<:Field})
  if is_perfect(base_ring(P))
    return is_reduced(P)
  end
  throw(NotImplementedError(:is_geometrically_reduced,"currently we can decide this only over a perfect base field"))
end

@attr Bool function is_geometrically_irreducible(P::AbsProjectiveScheme{<:Field, <:MPolyQuoRing{<:MPolyDecRingElem}})
  is_empty(P) && return false # the empty set is not irreducible
  I = defining_ideal(P)
  is_prime(I) || return false
  AI = absolute_primary_decomposition(I)
  return length(AI)==1
end

@attr Bool function is_geometrically_irreducible(X::AbsProjectiveScheme{<:Field, <:MPolyDecRing})
  return true # projective space is geometrically irreducible
end

@attr Bool function is_integral(X::AbsProjectiveScheme{<:Field, <:MPolyAnyRing})
  !is_empty(X) || return false
  if has_attribute(X,:is_reduced) && has_attribute(X,:is_irreducible)
     return get_attribute(X,:is_reduced) && get_attribute(X,:is_irreducible)
  end
  I = defining_ideal(X)
  return is_prime(I)
end

@attr Bool function is_geometrically_integral(X::AbsProjectiveScheme{<:Field, <:MPolyDecRing})
  return true # projective space is geometrically integral
end

@attr Bool function is_geometrically_integral(X::AbsProjectiveScheme{<:Field})
  is_integral(X) || return false
  throw(NotImplementedError(:is_geometrically_integral, "absolute primary decomposition is currently only available over the rationals"))
end

@attr Bool function is_geometrically_integral(X::AbsProjectiveScheme{<:QQField, <:MPolyQuoRing{<:MPolyDecRingElem}})
  is_integral(X) || return false
  I = defining_ideal(X)
  AI = absolute_primary_decomposition(I)
  @assert length(AI)==1 # it is prime since X is integral
  return AI[1][4]==1
end

@attr Bool function is_equidimensional(P::AbsProjectiveScheme)
  I = defining_ideal(P)
  return is_equidimensional(I)
end
