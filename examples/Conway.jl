module ConwayDB

using Nemo, Markdown

export conway_polynomial, hasconway_polynomial, conway_polynomials

@doc Markdown.doc"""
    conway_polynomial(p::Int, n::Int)
Return the Conway polynomial of degree `n` for the prime `p` - if known.
"""
function conway_polynomial(p::Int, n::Int)
  c = cglobal((:flint_conway_polynomials, Nemo.libflint), Cint)
  i = 1
  while unsafe_load(c, i) != 0 && unsafe_load(c, i) < p
    i += unsafe_load(c, i+1)+3 #prime, degree, coeffs (deg+1 many)
  end
  while unsafe_load(c, i) != 0 && unsafe_load(c, i+1) < n
    i += unsafe_load(c, i+1)+3
  end
  if unsafe_load(c, i) == 0 ||
     unsafe_load(c, i+1) != n
     error("Conway polynomial not known")
  end
  return polynomial(GF(p), [unsafe_load(c, j) for j=i+2:i+2+unsafe_load(c, i+1)+1])
end

@doc Markdown.doc"""
    hasconway_polynomial(p::Int, n::Int)
Test if the Conway polynomial of degree `n` is known for the prime `p`.
"""
function hasconway_polynomial(p::Int, n::Int)
  c = cglobal((:flint_conway_polynomials, Nemo.libflint), Cint)
  i = 1
  while unsafe_load(c, i) != 0 && unsafe_load(c, i) < p
    i += unsafe_load(c, i+1)+3 #prime, degree, coeffs (deg+1 many)
  end
  while unsafe_load(c, i) != 0 && unsafe_load(c, i+1) < n
    i += unsafe_load(c, i+1)+3
  end
  if unsafe_load(c, i) == 0 ||
     unsafe_load(c, i+1) != n
     return false
  end
  return true
end

@doc Markdown.doc"""
    conway_polynomials(p::Int)
Return an array with all Conway polynomials for the given prime `p`.
"""
function conway_polynomials(p::Int)
  c = cglobal((:flint_conway_polynomials, Nemo.libflint), Cint)
  i = 1
  while unsafe_load(c, i) != 0 && unsafe_load(c, i) < p
    i += unsafe_load(c, i+1)+3 #prime, degree, coeffs (deg+1 many)
  end
  if unsafe_load(c, i) == 0 
     error("no conway polynomial known")
  end

  k = GF(p)
  kt, t = PolynomialRing(k, cached = false)
  pols = typeof(t)[]

  while unsafe_load(c, i) == p
    push!(pols, kt([unsafe_load(c, j) for j=i+2:i+2+unsafe_load(c, i+1)]))

    i += unsafe_load(c, i+1)+3
  end
  return pols
end

end

using .ConwayDB
export conway_polynomials, conway_polynomial, hasconway_polynomial
