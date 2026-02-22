"""
    bernstein_basis_polynomial(R, i, n)

Return the Bernstein basis polynomial B_{i,n}(t) = binomial(n,i) * t^i * (1-t)^(n-i).

- `R`: univariate polynomial ring
- `i`: basis index (0 ≤ i ≤ n)
- `n`: degree (n ≥ 0)
"""
function bernstein_basis_polynomial(R::PolyRing, i::Int, n::Int)
  @assert 0 ≤ i ≤ n "i must be in 0..n"
  t = gen(R)
  binomial(n,i) * (t^i) * ((1 - t)^(n-i))
end

"""
    bernstein_basis_polynomial(i, n)

Return the Bernstein basis polynomial B_{i,n}(t) = binomial(n,i) * t^i * (1-t)^(n-i).

- `i`: basis index (0 ≤ i ≤ n)
- `n`: degree (n ≥ 0)
"""
function bernstein_basis_polynomial(i::Int, n::Int)
  R, _ = polynomial_ring(ZZ, :t)
  bernstein_basis_polynomial(R, i, n)
end

"""
    bernstein_polynomial(R, coeffs)

Return the Bernstein polynomial which is the linear combination of the Bernstein basis polynomials with the given coefficients.

- `R`: univariate polynomial ring
- `coeffs`: vector of rational coefficients
"""
function bernstein_polynomial(R::QQPolyRing, coeffs::Vector{QQFieldElem})
  n = length(coeffs)-1
  sum([ coeffs[i+1] * bernstein_basis_polynomial(R, i, n) for i in 0:n ])
end

"""
    bernstein_polynomial(coeffs)

Return the Bernstein polynomial which is the linear combination of the Bernstein basis polynomials with the given coefficients.

- `coeffs`: vector of rational coefficients
"""
function bernstein_polynomial(coeffs::Vector{QQFieldElem})
  R, _ = polynomial_ring(QQ, :t)
  bernstein_polynomial(R, coeffs)
end

"""
    bernstein_01_approximation(f, n)

Return a Bernstein polynomial of degree `n` which approximates the polynomial `f` on the interval [0,1].
The sequence of these polynomials converges uniformly to `f` on [0,1] if `n` goes to infinity.
However, the convergence is rather slow.

- `f`: univariate polynomial with rational coefficients
- `n`: degree
"""
function bernstein_01_approximation(f::QQPolyRingElem, n::Int)
  @assert 0 < n "n must be positive"
  R = parent(f)
  g = bernstein_polynomial(R, [ f(i//n) for i in 0:n ])
  return g, f-g
end

function bezier_curve(p::Vector{Vector{QQFieldElem}})
  n = length(p) # deg+1
  d = length(p[1])
  for k in 1:d
    @show bernstein_polynomial([ p[i][k] for i in 1:n ] )
  end
end

function _draw_edge_sequence_bernstein(io, pts::Vector{_Point}, scale; color::String="black")
  # Use formula from
  # https://web.archive.org/web/20131225210855/http://people.sc.fsu.edu/~jburkardt/html/bezier_interpolation.html
  # linked on Wikipedia to achieve curve through points.
  ptsmod2 = 1//6*(-5*pts[1]+18*pts[1]-9*pts[3]+2*pts[4])
  ptsmod3 = 1//6*(-5*pts[4]+18*pts[3]-9*pts[2]+2*pts[1])
  ptsmod = [pts[1], ptsmod2, ptsmod3, pts[4]]
   f1 = bernstein_polynomial([pt.xcoord for pt in pts])
   f2 = bernstein_polynomial([pt.ycoord for pt in pts])
   pts100 = [_Point(Oscar.evaluate(f1, i//100), Oscar.evaluate(f2, i//100)) for i in 1:100]
   for i in 1:length(pts100)-1
      Oscar._draw_edge_tikz(io, pts100[i], pts100[i+1], scale; color="blue")
   end
end
