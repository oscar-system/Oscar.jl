export is_equidimensional
export singular_locus
export reduced_scheme
export is_smooth

function singular_locus(X::AbsSpec{<:Ring, <:MPolyQuoLocalizedRing})
  comp = _singular_locus_with_decomposition(X)
  if length(comp) == 0 
    return subscheme(X, one(OO(X)))
  end
  R = base_ring(OO(X))
  return Spec(R, prod([modulus(quotient_ring(OO(Y))) for Y in comp]), inverted_set(OO(X)))
end

function _singular_locus_with_decomposition(X::AbsSpec{<:Ring, <:MPolyQuoLocalizedRing})
  I = modulus(OO(X))
  result = typeof(X)[]

  P = []
  if has_attribute(X, :is_equidimensional) && is_equidimensional(X) 
    push!(P, saturated_ideal(I))
  else 
    ## TODO: This is the scheme case, for varieties case use ..._radical
    P = equidimensional_decomposition_weak(saturated_ideal(I))
  end

  if length(P)==1
    d = dim(X)
    R = base_ring(OO(X))
    n = ngens(R)
    J = saturated_ideal(I)
    M = jacobi_matrix(gens(J))
    K = ideal(R, minors(M, n-d))
    one(OO(X)) in OO(X)(K) && return result
    return [subscheme(X, K)]
  else
    components = [subscheme(X, J) for J in P]
    for Y in components
      for Z in components
        if !(Y == Z)
          W = intersect(Y, Z)
          if !isempty(W) 
            push!(result, W)
          end
        end
      end
    end
    for Y in components
      result = vcat(result, singular_locus(Y))
    end
  end
  return result
end

@attr Bool function is_equidimensional(X::AbsSpec{<:Ring, <:MPolyQuoLocalizedRing})
  I = modulus(OO(X))
  P = equidimensional_decomposition_weak(saturated_ideal(I))
  length(P) < 2 && return true
  return false
end

@attr function reduced_scheme(X::AbsSpec{<:Ring, <:MPolyQuoLocalizedRing})
  I = modulus(OO(X))
  J = radical(saturated_ideal(I))
  return Spec(base_ring(J), J, inverted_set(OO(X)))
end

@attr function is_reduced(X::AbsSpec{<:Ring, <:MPolyQuoLocalizedRing})
  I = modulus(OO(X))
  J = saturated_ideal(I)
  return is_reduced(quo(base_ring(OO(X)), J)[1])
end


########################################################################
# Additional functionality for fractions                               #
########################################################################

function derivative(f::MPolyLocalizedRingElem, i::Int)
  num = derivative(numerator(f), i)*denominator(f) - derivative(denominator(f), i)*numerator(f)
  den = denominator(f)^2
  g = gcd(num, den)
  return parent(f)(divexact(num, g), divexact(den, g), check=false)
end

function derivative(f::MPolyQuoLocalizedRingElem, i::Int)
  num = derivative(lifted_numerator(f), i)*lifted_denominator(f) - derivative(lifted_denominator(f), i)*lifted_numerator(f)
  den = lifted_denominator(f)^2
  g = gcd(num, den)
  return parent(f)(divexact(num, g), divexact(den, g), check=false)
end

function jacobi_matrix(f::MPolyLocalizedRingElem)
  L = parent(f)
  n = nvars(base_ring(L))
  return matrix(L, n, 1, [derivative(f, i) for i=1:n])
end

function jacobi_matrix(g::Vector{<:MPolyLocalizedRingElem})
  R = parent(g[1])
  n = nvars(base_ring(R))
  @assert all(x->parent(x) == R, g)
  return matrix(R, n, length(g), [derivative(x, i) for i=1:n for x = g])
end

function jacobi_matrix(f::MPolyQuoLocalizedRingElem)
  L = parent(f)
  n = nvars(base_ring(L))
  return matrix(L, n, 1, [derivative(f, i) for i=1:n])
end

function jacobi_matrix(g::Vector{<:MPolyQuoLocalizedRingElem})
  L = parent(g[1])
  n = nvars(base_ring(L))
  @assert all(x->parent(x) == L, g)
  return matrix(L, n, length(g), [derivative(x, i) for i=1:n for x = g])
end

########################################################################
# Smoothness test based on projective modules.
#
# The routine checks whether the module for the cotangent sheaf Ω¹(X)
# is locally free over 𝒪(X) and returns `true` if this is the case. 
########################################################################

@attr Bool function is_smooth(X::AbsSpec{<:Field, <:MPolyQuoLocalizedRing})
  R = base_ring(OO(X))
  characteristic(base_ring(X)) == 0 || error("method not implemented in positive characteristic")
  L = localized_ring(OO(X))
  I = modulus(OO(X))
  f = gens(Oscar.pre_saturated_ideal(I))
  Df = jacobi_matrix(f)

  A = map_entries(x->OO(X)(x), Df)
  success, _, _ = Oscar._is_projective_without_denominators(A)
  return success
end

is_smooth(X::AbsSpec{<:Field, <:MPolyRing}) = true
is_smooth(X::AbsSpec{<:Field, <:MPolyLocalizedRing}) = true

@attr Bool function is_smooth(X::AbsSpec{<:Field, <:MPolyQuo})
  R = base_ring(OO(X))
  characteristic(base_ring(X)) == 0 || error("method not implemented in positive characteristic")
  I = modulus(OO(X))
  f = gens(I)
  Df = jacobi_matrix(f)

  A = map_entries(x->OO(X)(x), Df)
  success, _, _ = Oscar._is_projective_without_denominators(A)
  return success
end

