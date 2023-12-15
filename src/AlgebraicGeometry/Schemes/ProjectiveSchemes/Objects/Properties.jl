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

Check whether the scheme `X` is smooth.

# Algorithms

There are three possible algorithms for checking smoothness, determined
by the value of the keyword argument `algorithm`:
  * `:covered` - converts to a covered scheme,
  * `:jacobian` - uses the Jacobian criterion for projective varieties,
    see Exercise 4.2.10 of [Liu06](@cite),
  * `:affine_cone` - checks that the affine cone is smooth outside the origin.

The `:jacobian` algorithm first checks that the scheme is integral,
which can be expensive.

By default, if the base ring is a field and the scheme has already been
computed to be integral, then the `:jacobian` algorithm is used.
Otherwise, the `:affine_cone` algorithm is used.

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

julia> is_smooth(C; algorithm=:jacobian)
false
```
"""
function is_smooth(P::AbsProjectiveScheme; algorithm=:default)
  get_attribute!(P, :is_smooth) do
    if algorithm == :default
      if base_ring(P) isa Field && has_attribute(P, :is_integral) && is_integral(P)
        algorithm = :jacobian
      else
        algorithm = :affine_cone
      end
    end
    if algorithm == :covered
      return is_smooth(covered_scheme(P))
    elseif algorithm == :affine_cone
      aff = affine_cone(P)[1]
      origin = subscheme(aff, coordinates(aff))
      return singular_locus_reduced(aff)[1] == origin
    elseif algorithm == :jacobian
      if !(base_ring(P) isa Field)
        throw(NotImplementedError(:is_smooth, "jacobi criterion not implemented when base ring not a field"))
      end
      if !(is_integral(P))
        throw(NotImplementedError(:is_smooth, "jacobi criterion not implemented when scheme not integral"))
      end
      R = base_ring(homogeneous_coordinate_ring(P))
      I = defining_ideal(P)
      mat = jacobi_matrix(R, gens(I))
      sing_locus = I + ideal(R, minors(mat, dim(ambient_space(P)) - dim(P)))
      sing_subscheme = subscheme(ambient_space(P), sing_locus)
      return isempty(sing_subscheme)
    end
  end
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
