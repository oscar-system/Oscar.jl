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
    is_smooth(P::AbsProjectiveScheme; algorithm=:default) -> Bool

Check whether the scheme `P` is smooth.

# Algorithms

There are three possible algorithms for checking smoothness, determined
by the value of the keyword argument `algorithm`:
  * `:projective_jacobian` - uses the Jacobian criterion for projective
    varieties, see Exercise 4.2.10 of [Liu06](@cite),
  * `:covered_jacobian` - uses covered version of the Jacobian criterion,
  * `:affine_cone` - checks that the affine cone is smooth outside the origin.

The `:projective_jacobian` and the `:covered` algorithms only work for equidimensional schemes. The algorithms first check for equidimensionality, which can be expensive.

By default, if the base ring is a field, then we first compute whether the scheme is equidimensional. If yes, then the `:projective_jacobian` algorithm is used. Otherwise, the `:affine_cone` algorithm is used.

If you already know that the scheme is equidimensional, then you can
avoid recomputing that by writing
  `set_attribute!(P, :is_equidimensional, true)`
before checking for smoothness.

We explain why the algorithm of `:affine_cone` works for arbitrary schemes over arbitrary base schemes. By Remark 13.38(1) of [GW20](@cite), the morphism from the pointed affine cone to $P$ is locally the morphism $\mathbb{A}_U^1 \setminus \{0\} \to U$, where $U$ is an affine open of $P$ and $\mathbb{A}_U^1$ is the relative affine 1-space over $U$. By Definition 6.14(1) of [GW20](@cite), the Jacobian matrix for $U$ differs from the Jacobian matrix for $P$ only by a column containing zeros, implying that the ranks of the Jacobian matrices are the same. Therefore, $P$ is smooth if and only if the affine cone is smooth outside the origin.

# Examples
```jldoctest
julia> A, (x, y, z) = grade(QQ["x", "y", "z"][1]);

julia> B, _ = quo(A, ideal(A, [x^2 + y^2]));

julia> C = projective_scheme(B)
Projective scheme
  over rational field
defined by ideal (x^2 + y^2)

julia> is_smooth(C)
false

julia> is_smooth(C; algorithm=:covered_jacobian)
false
```
"""
function is_smooth(P::AbsProjectiveScheme; algorithm=:default)
  get_attribute!(P, :is_smooth) do
    if is_empty(P)
      return true
    end

    if algorithm == :default
      if (
        base_ring(P) isa Field
        && is_equidimensional(P)
      )
        algorithm = :projective_jacobian
      else
        algorithm = :affine_cone
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
      return _jacobian_criterion(covered_scheme(P))
    elseif algorithm == :affine_cone
      aff, _ = affine_cone(P)
      sing, _ = singular_locus(aff)
      origin = ideal(gens(ambient_coordinate_ring(sing)))
      return isone(saturation(saturated_ideal(defining_ideal(sing)), origin))
    elseif algorithm == :projective_jacobian
      return _projective_jacobian_criterion(P)
    end
  end
end

function _projective_jacobian_criterion(P::AbsProjectiveScheme)
  if !(base_ring(P) isa Field)
    throw(NotImplementedError(
      :is_smooth,
      "projective Jacobian criterion only implemented when base ring is a field"
    ))
  end
  if !(is_equidimensional(P))
    throw(NotImplementedError(
      :is_smooth,
      "projective Jacobian criterion only implemented when scheme is equidimensional"
    ))
  end
  R = base_ring(homogeneous_coordinate_ring(P))
  I = defining_ideal(P)
  mat = jacobian_matrix(R, gens(I))
  sing_locus = ideal(R, minors(mat, codim(P)))
  return dim(sing_locus + I) <= 0
end

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
