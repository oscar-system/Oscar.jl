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
    function is_smooth(P::AbsProjectiveScheme; algorithm=:default) -> Bool

Check whether the scheme `P` is smooth.

For schemes that are not equidimensional, use the optional argument
  `algorithm=:affine_cone`.
If you already know that the scheme is equidimensional and wish to avoid
computing that, write
  !set_attribute(P, :is_equidimensional, true)
before checking for smoothness.

# Algorithms

There are three possible algorithms for checking smoothness, determined
by the value of the keyword argument `algorithm`:
  * `:projective_jacobian` - uses the Jacobian criterion for projective
    varieties, see Exercise 4.2.10 of [Liu06](@cite),
  * `:covered` - converts to a covered scheme,
  * `:affine_cone` - checks that the affine cone is smooth outside the origin.

The `:projective_jacobian` and the `:covered` algorithms only work for equidimensional schemes. The algorithms first check for equidimensionality, which can be expensive.

By default, if the base ring is a field, then we first compute whether the scheme is equidimensional. If yes, then the `:projective_jacobian` algorithm is used. Otherwise, the `:affine_cone` algorithm is used.

If you already know that the scheme is equidimensional, then you can
avoid recomputing that by writing
  `set_attribute!(P, :is_equidimensional, true)`
before checking for smoothness.

# Examples
```jldoctest
julia> A, (x, y, z) = grade(QQ["x", "y", "z"][1]);

julia> B, _ = quo(A, ideal(A, [x^2 + y^2]));

julia> C = projective_scheme(B)
Projective scheme
  over rational field
defined by ideal(x^2 + y^2)

julia> is_smooth(C)
false

julia> is_smooth(C; algorithm=:covered_jacobian)
false
```
"""
function is_smooth(P::AbsProjectiveScheme; algorithm=:default)
  get_attribute!(P, :is_smooth) do
    if algorithm == :default
      if (
        base_ring(P) isa Field
        && has_attribute(P, :is_equidimensional)
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
         + " of the following: " + join(algorithms, ", ") + "."
      ))
    end

    if algorithm == :covered_jacobian
      return _jacobian_criterion(covered_scheme(P))
    elseif algorithm == :affine_cone
      aff = affine_cone(P)[1]
      singular_locus(aff)[1]
      return dim(singular_locus(aff)[1]) <= 0
    elseif algorithm == :projective_jacobian
      return _projective_jacobian_criterion(P)
    end
  end
end

# Use decomposition_info
function _projective_jacobian_criterion(P::AbsProjectiveScheme)
  if !(base_ring(P) isa Field)
    throw(NotImplementedError(
      :is_smooth,
      "projective jacobian criterion not implemented when base ring not a field"
    ))
  end
  if !(is_equidimensional(P))
    throw(NotImplementedError(
      :is_smooth,
      "projective jacobian criterion not implemented when scheme not equidimensional"
    ))
  end
  R = base_ring(homogeneous_coordinate_ring(P))
  I = defining_ideal(P)
  mat = jacobi_matrix(R, gens(I))
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
