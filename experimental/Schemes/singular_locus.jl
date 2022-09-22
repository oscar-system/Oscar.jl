export is_equidimensional
export singular_locus_with_decomposition
export reduced_scheme
export relative_polar_loci, relative_polar_multiplicities
export polar_multiplicities, polar_loci
export random_linear_forms

function singular_locus(X::Spec)
  comp = singular_locus_with_decomposition(X)
  if length(comp) == 0 
    return subscheme(X, one(OO(X)))
  end
  R = base_ring(OO(X))
  return Spec(R, prod([modulus(OO(Y)) for Y in comp]), inverted_set(OO(X)))
end

function singular_locus_with_decomposition(X::Spec)
  I = localized_modulus(OO(X))
  result = typeof(X)[]
  if is_equidimensional(X)
    d = dim(X)
    R = base_ring(OO(X))
    n = ngens(R)
    J = saturated_ideal(I)
    M = jacobi_matrix(gens(J))
    K = ideal(R, minors(M, n-d))
    one(OO(X)) in OO(X)(K) && return result
    return [subscheme(X, K)]
  else
    P = primary_decomposition(saturated_ideal(I))
    components = [subscheme(X, J[2]) for J in P]
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

@attr Bool function is_equidimensional(X::Spec)
  I = localized_modulus(OO(X))
  P = primary_decomposition(saturated_ideal(I))
  length(P) == 0 && return true
  d = dim(P[1][1])
  all(x->dim(x[1])==d, P[2:end]) && return true
  return false
end

@attr function reduced_scheme(X::Spec)
  I = localized_modulus(OO(X))
  J = radical(saturated_ideal(I))
  return Spec(base_ring(J), J, inverted_set(OO(X)))
end

@attr function is_reduced(X::Spec)
  I = localized_modulus(OO(X))
  J = saturated_ideal(I)
  return is_reduced(quo(base_ring(OO(X)), J)[1])
end

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

function random_linear_forms(R::MPolyRing, d::Int)
  n = ngens(R)
  kk = coefficient_ring(R)
  A = MatrixSpace(kk, d, n)
  M = zero(A)
  for i in 1:nrows(A)
    for j in 1:ncols(A)
      M[i,j] = rand(kk, -10:10)
    end
  end
  b = [sum([R[j]*M[i,j] for j in 1:n]) for i in 1:d]
  return b
end

function polar_loci(X::Spec, a::RingElem)
  R = base_ring(OO(X))
  kk = coefficient_ring(R)
  d = dim(X)
  n = ngens(R)
  A = MatrixSpace(kk, d, n)
  M = zero(A)
  for i in 1:nrows(A)
    for j in 1:ncols(A)
      M[i,j] = rand(kk, -10:10)
    end
  end
  b = [sum([R[j]*M[i,j] for j in 1:n]) for i in 1:d]
  return relative_polar_loci(X, a, b)
end

function polar_multiplicities(X::Spec, a::RingElem)
  R = base_ring(OO(X))
  kk = coefficient_ring(R)
  d = dim(X)
  n = ngens(R)
  A = MatrixSpace(kk, d, n)
  M = zero(A)
  for i in 1:nrows(A)
    for j in 1:ncols(A)
      M[i,j] = rand(kk, -10:10)
    end
  end
  b = [sum([R[j]*M[i,j] for j in 1:n]) for i in 1:d]
  return relative_polar_multiplicities(X, a, b)
end

function relative_polar_loci(X::Spec, a::RingElem, b::Vector{<:RingElem})
  L = localized_ring(OO(X))
  f = L(a)
  l = L.(b)
  g = gens(localized_modulus(OO(X)))
  is_equidimensional(X) && is_reduced(X) || error("method not implemented for $X unreduced or of mixed dimension.")
  Z = singular_locus(X)
  Dl = jacobi_matrix(l)
  Df = jacobi_matrix(f)
  Dg = jacobi_matrix(g)
  n = ngens(base_ring(OO(X)))
  d = dim(X)
  k = n-d
  result = typeof(X)[]
  for i in 1:d
    A = hcat(Df, Dl[:,1:i], Dg)
    I = ideal(L, minors(A, k+i+1))
    P = subscheme(X, I)
    push!(result, P)
  end
  return result
end

function relative_polar_multiplicities(X::Spec, a::RingElem, b::Vector{<:RingElem})
  P = relative_polar_loci(X, a, b)
  PI = [localized_modulus(OO(Y))+ideal(localized_ring(OO(Y)), b[1:i]) for (Y, i) in zip(P, 1:length(b))]
  m = multiplicity.(PI)
  return reverse(m)
end

function length(X::Spec)
  R = base_ring(OO(X))
  I = saturated_ideal(localized_modulus(OO(X)))
  leadI = leading_ideal(I)
  dim(leadI) <= 0 || error("scheme is not zero dimensional")
  d=0
  result = 0
  mons = [x for x in all_monomials(R, d) if !(x in leadI)]
  while length(mons) > 0
    result = result + length(mons)
    d = d+1
    mons = [x for x in all_monomials(R, d) if !(x in leadI)]
  end
  return result
end



