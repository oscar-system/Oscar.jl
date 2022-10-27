export is_equidimensional
export singular_locus, singular_locus_reduced
export reduced_scheme
export is_smooth

QuoRing = Union{MPolyQuoLocalizedRing, 
                MPolyQuo
               }

NonQuoRing = Union{MPolyRing, MPolyLocalizedRing
                  }

### TODO: The following two functions need to be made type-sensitive 
###       and reduced=true needs to be set automatically for varieties 
###       as soon as not only schemes, but also varieties as special 
###       schemes have been introduced in OSCAR
function singular_locus(X::AbsSpec{<:Ring, <:QuoRing})
  comp = _singular_locus_with_decomposition(X,false)
  if length(comp) == 0 
    return subscheme(X, ideal(OO(X),one(OO(X))))
  end
  R = base_ring(OO(X))
  I= prod([modulus(quotient_ring(OO(Y))) for Y in comp])
#  return Spec(R, I, inverted_set(OO(X)))
  return subscheme(X,I)
end

function singular_locus_reduced(X::AbsSpec{<:Ring, <:QuoRing})
  comp =  _singular_locus_with_decomposition(X, true)
  I= ideal(localized_ring(OO(X)),one(localized_ring(OO(X))))
  for Z in comp
    I = intersect(I, modulus(OO(Z)))
  end     
  return subscheme(X,I)
end

function _singular_locus_with_decomposition(X::AbsSpec{<:Ring, <:QuoRing}, reduced::Bool=true)
  I = saturated_ideal(modulus(OO(X)))
  result = typeof(X)[]

  P = []
  if has_attribute(X, :is_equidimensional) && is_equidimensional(X) 
    push!(P, I)
  else 
    if reduced
      P = equidimensional_decomposition_radical(I)
    else
      P = equidimensional_decomposition_weak(I)
    end
  end

  if length(P)==1
    d = dim(X)
    R = base_ring(I)
    n = nvars(R) 
    M = _jacobi_matrix_modulus(X)
    minvec = minors(M, n-d)
    J = ideal(R, minvec)
    one(OO(X)) in OO(X)(J) && return result
    return [subscheme(X, J)]
  else
    components = [subscheme(X, J) for J in P]
    for i in 1:length(components)
      for j in (i+1):length(components)
        W = intersect(components[i],components[j])
        if !isempty(W) 
          push!(result, W)
        end
      end
    end
    for Y in components
      result = vcat(result, singular_locus(Y))
    end
  end
  return result
end

### only for users' convenience: appropriate return value for NonQuoRings

function singular_locus(X::AbsSpec{<:Ring, <:NonQuoRing})
  return subscheme(X,ideal(OO(X),one(OO(X))))
end

function singular_locus_reduced(X::AbsSpec{<:Ring, <:NonQuoRing})
  return subscheme(X,ideal(OO(X),one(OO(X))))
end

### to keep number of cases down, create a dummy for 
### saturated_ideal for the affine case
function saturated_ideal(I::MPolyIdeal)
  return(I)
end

@attr Bool function is_equidimensional(X::AbsSpec{<:Ring, <:QuoRing})
  I = modulus(OO(X))
  P = equidimensional_decomposition_radical(saturated_ideal(I))
  length(P) < 2 && return true
  return false
end

@attr function reduced_scheme(X::AbsSpec{<:Ring, <:QuoRing})
  I = modulus(OO(X))
  J = radical(pre_saturated_ideal(I))
  return Spec(base_ring(J), J, inverted_set(OO(X)))
end

@attr function is_reduced(X::AbsSpec{<:Ring, <:QuoRing})
  I = saturated_ideal(modulus(OO(X)))
  return is_reduced(quo(base_ring(I), I)[1])
end

### only for users' convenience: appropriate return value for NonQuoRings

@attr Bool function is_equidimensional(X::AbsSpec{<:Ring, <:NonQuoRing})
  return true
end

@attr AbsSpec function reduced_scheme(X::AbsSpec{<:Ring, <:NonQuoRing})
  return X
end

@attr Bool function is_reduced(X::AbsSpec{<:Ring, <:NonQuoRing})
  return true
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

function _jacobi_matrix_modulus(X::AbsSpec{<:Ring, <:MPolyQuoLocalizedRing})
  g = gens(modulus(quotient_ring(OO(X))))
  L = base_ring(quotient_ring(OO(X)))
  n = nvars(L)
  M = matrix(L, n, length(g),[derivative(f,i) for i=1:n for f in g])
  return M
end

function _jacobi_matrix_modulus(X::AbsSpec{<:Ring, <:MPolyQuo})
  g = gens(modulus(OO(X)))
  L = base_ring(OO(X))
  n = nvars(L)
  M = matrix(L, n , length(g), [derivative(f,i) for i=1:n for f in g])
  return M
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

