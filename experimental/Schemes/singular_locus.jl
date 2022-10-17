export is_equidimensional
export singular_locus
export reduced_scheme
export is_smooth

QuoRing = Union{MPolyQuoLocalizedRing, 
                MPolyQuo
               }

### TODO: The following two functions need to be made type-sensitive 
###       and reduced=true needs to be set for varieties 
###       as soon as not only schemes, but also varieties as special 
###       schemes have been introduced
function singular_locus(X::AbsSpec{<:Ring, <:QuoRing})
  comp = _singular_locus_with_decomposition(X,reduced=false)
  if length(comp) == 0 
    return subscheme(X, one(OO(X)))
  end
  R = base_ring(OO(X))
  I= prod([modulus(quotient_ring(OO(Y))) for Y in comp])
  return Spec(R, I, inverted_set(OO(X)))
end

function singular_locus_reduced(X::AbsSpec{<:Ring, <:QuoRing})
  comp =  _singular_locus_with_decomposition(X, reduced=true)
  I= ideal(ambient_ring(X),one(ambient_ring(X)))
  for Z in comp
    I = intersect(I, modulus(quotient_ring(Z)))
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
    R = base_ring(OO(X))
    n = ngens(R)
    M = jacobi_matrix(gens(I))
    J = ideal(R, minors(M, n-d))
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

#####################################################################
### TODO: create a 'dummy' saturated_ideal for the global case
###       which is just a NOOP
#####################################################################
@attr Bool function is_equidimensional(X::AbsSpec{<:Ring, <:QuoRing})
  I = modulus(OO(X))
  P = equidimensional_decomposition_radical(saturated_ideal(I))
  length(P) < 2 && return true
  return false
end

@attr function reduced_scheme(X::AbsSpec{<:Ring, <:QuoRing})
  I = modulus(OO(X))
  J = radical(saturated_ideal(I))
  return Spec(base_ring(J), J, inverted_set(OO(X)))
end

@attr function is_reduced(X::AbsSpec{<:Ring, <:QuoRing})
  I = modulus(OO(X))
  J = saturated_ideal(I)
  return is_reduced(quo(base_ring(J), J)[1])
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

