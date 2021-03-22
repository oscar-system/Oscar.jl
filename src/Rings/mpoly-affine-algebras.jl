
export normalize_with_delta
export noether_normalization, normalization
export isreduced, subalgebra_membership
export hilbert_series, hilbert_series_reduced, hilbert_series_expanded, hilbert_function, hilbert_polynomial, degree
export issurjective, isinjective, isbijective, inverse, preimage, isfinite

##############################################################################
#
# Data associated to affine algebras
#
##############################################################################

@doc Markdown.doc"""
    dim(A::MPolyQuo)

Returns the dimension of `A`.
"""
function dim(A::MPolyQuo) 
  I = A.I
  return dim(I)
end

##############################################################################
#
# Data associated to affine algebras
#
##############################################################################

@doc Markdown.doc"""
    hilbert_series(A::MPolyQuo)

Given a graded affine algebra $A$ over a field $K$, return a pair $(p,q)$, say, 
of univariate polynomials in $t$ with integer coefficients
such that: If $A = R/I$, where $R$ is a multivariate polynomial ring in $n$
variables over $K$, with positive integer weights $w_1, \dots, w_n$ assigned to the variables,
and $I$ is a homogeneous ideal of $R$ when we grade $R$ according to the corresponding 
weighted degree, then $p/q$ represents the Hilbert series of $A$ as a rational function
with denominator $q = (1-t^{w_1})\cdots (1-t^{w_n})$. See also 
`hilbert_series_reduced(A::MPolyQuo)`.
"""
function hilbert_series(A::MPolyQuo)
   H = HilbertData(A.I)
   return hilbert_series(H,1)
end

@doc Markdown.doc"""
    hilbert_series_reduced(A::MPolyQuo)

Given a graded affine algebra $A$ over a field $K$, return a pair $(p,q)$, say, 
of univariate polynomials in $t$ with integer coefficients such that: 
If $A = R/I$, where $R$ is a multivariate polynomial ring in $n$
variables over $K$, with positive integer weights $w_1, \dots, w_n$ assigned to the variables,
and $I$ is a homogeneous ideal of $R$ when we grade $R$ according to the corresponding 
weighted degree, then $p/q$ represents the Hilbert series of $A$ as a rational function
written in lowest terms. See also `hilbert_series(A::MPolyQuo)`.
"""
function hilbert_series_reduced(A::MPolyQuo)
   H = HilbertData(A.I)
   return hilbert_series(H,2)
end

@doc Markdown.doc"""
    hilbert_series_expanded(A::MPolyQuo, d::Int)

Given an affine algebra $A = R/I$ over a field $K$ and an integer $d\geq 0$, return the
Hilbert series of $A$ to precision $d$. 
"""
function hilbert_series_expanded(A::MPolyQuo, d::Int)
   H = HilbertData(A.I)  
   return hilbert_series_expanded(H, d)
end

@doc Markdown.doc"""
    hilbert_function(A::MPolyQuo, d::Int)

Given an affine algebra $A = R/I$ over a field $K$ and an integer $d\geq 0$, return the value
$H(A, d)$, where $H(A, \underline{\phantom{d}}): \N \rightarrow \N, d \mapsto \dim_K A_d$ is 
the Hilbert function of $A$.
"""
function hilbert_function(A::MPolyQuo, d::Int)
   H = HilbertData(A.I)
   return hilbert_function(H, d)
end

@doc Markdown.doc"""
    hilbert_polynomial(A::MPolyQuo)
"""
function hilbert_polynomial(A::MPolyQuo)
   H = HilbertData(A.I)
   return hilbert_polynomial(H)
end

@doc Markdown.doc"""
    degree(A::MPolyQuo)
"""
function degree(A::MPolyQuo)
   H = HilbertData(A.I)
   return degree(H)
end

##############################################################################
#
# Properties of affine algebras
#
##############################################################################

@doc Markdown.doc"""
    isreduced(A::MPolyQuo)

Returns `true` if `A` is reduced, `false` otherwise.
"""
function isreduced(A::MPolyQuo) 
  I = A.I
  return I == radical(I)
end

@doc Markdown.doc"""
    isnormal(A::MPolyQuo)
"""
function isnormal(A::MPolyQuo)
  _, _, d = normalize_with_delta(A)
  return d == 1
end

##############################################################################
#
# Algebra Containment
#
##############################################################################

# helper function for the containement problem, surjectivity and preimage
function _containement_helper(R::MPolyRing, N::Int, M::Int, I::MPolyIdeal, W::Vector)
   T = PolynomialRing(base_ring(R), N + M, ordering = :lex)[1]
   phi = hom(R, T, [gen(T, i) for i in 1:M])

   # Groebner computation
   A = phi.(gens(I))
   B = [phi(W[i])-gen(T, M + i) for i in 1:N]
   J = ideal(T, vcat(A, B))
   G = groebner_basis(J, complete_reduction = true, ord = :lex)
   return (T, phi, G)
end

# helper function to obtain information about qring status and 
# convert elements, if needed
function _ring_helper(r, f, V)
   if isdefined(r, :I)
      R = r.R
      I = r.I
      W = [v.f for v in V]
      F = f.f
   else
      R = r
      I = ideal(R, [R(0)])
      W = [v for v in V]
      F = f
   end
   return (R, I, W, F)
end


@doc Markdown.doc"""
    subalgebra_membership(f::S, v::Vector{S}) where S <: Union{MPolyElem{U}, MPolyQuoElem{U}} where U

Returns a tuple `(true`, h)` if `f` is contained in the subalgebra of `A` generated by the elements contained in `v` with relation defined by `h`, `(false, 0)` otherwise.
"""
function subalgebra_membership(f::S, v::Vector{S}) where S <: Union{MPolyElem{U}, MPolyQuoElem{U}} where U
   @assert !isempty(v)
   r = parent(f)
   (R, I, W, F) = _ring_helper(r, f, v)
   n = length(W)
   m = ngens(R)

   # Build auxilliary objects
   (T, phi, G) = _containement_helper(R, n, m, I, W)
   TT, _ = PolynomialRing(base_ring(T), ["t_$i" for i in 1:n])
   
   # Check containement
   D = divrem(phi(F), G)
   if leading_monomial(D[2]) < gen(T, m)
      A = [zero(TT) for i in 1:m]
      B = [gen(TT, i) for i in 1:n]
      psi = hom(T, TT, vcat(A, B))
      return (true, psi(D[2]))
   else
      return (false, TT(0))
   end
end

##############################################################################
#
# Properties of maps of affine algebras
#
##############################################################################

@doc Markdown.doc"""
    isinjective(F::AlgHom)

Return `true` if `F` is injective, `false` otherwise.
"""
function isinjective(F::AlgHom)
   G = gens(kernel(F))
   for i in 1:length(G)
      !iszero(G[i]) && return false
   end
   return true
end

# Helper function related to the computation of surjectivity, preimage etc.
# Stores the necessary data in F.surj_help
function _surj_helper(F::AlgHom)
   r = domain(F)
   s = codomain(F)
   n = ngens(r)
   m = ngens(s)

   if !isdefined(F, :surj_helper)
      (S, I, W, _) = _ring_helper(s, zero(s), F.image)
      # Build auxilliary objects
      (T, inc, G) = _containement_helper(S, n, m, I, W)
      D = divrem([gen(T, i) for i in 1:m], G)
      A = [zero(r) for i in 1:m]
      B = [gen(r, i) for i in 1:n]
      pr = hom(T, r, vcat(A, B))
      F.surj_helper = (T, inc, pr, G, D)
   end
   return F.surj_helper
end


@doc Markdown.doc"""
    issurjective(F::AlgHom)

Return `true` if `F` is surjective, `false` otherwise.
"""
function issurjective(F::AlgHom)

   # Compute data necessary for computation
   r = domain(F)
   s = codomain(F)
   n = ngens(r)
   m = ngens(s)
   (T, _, _, _, D) = _surj_helper(F)

   # Check if map is surjective

   for i in 1:m
      if !(leading_monomial(D[i][2]) < gen(T, m))
         return false
      end
   end
   return true
end

@doc Markdown.doc"""
    isbijective(F::AlgHom)

Return `true` if `F` is bijective, `false` otherwise.
"""
function isbijective(F::AlgHom)
  return isinjective(F) && issurjective(F)
end

@doc Markdown.doc"""
    isfinite(F::AlgHom)

Return `true` if `F` is finite, `false` otherwise.
"""
function isfinite(F::AlgHom)
  (T, _, _, G, _) = _surj_helper(F)
  # Find all elements with leading monomial which contains the 
  # variables x_i.
  s = codomain(F)
  m = ngens(s)
  L = leading_monomial.(G)

  # Check if for all i, powers of x_i occur as leading monomials
  N = Vector{Int}()
  for i in 1:length(L)
     exp = exponent_vector(L[i], 1)
     f = findall(x->x!=0, exp) 
     length(f) == 1 && f[1] <= m && union!(N, f)
  end

  return length(N) == m
end

##############################################################################
#
# Inverse of maps of affine algebras and preimages of elements
#
##############################################################################

function inverse(F::AlgHom)
   !isinjective(F) && error("Homomorphism is not injective")
   !issurjective(F) && error("Homomorphism is not surjective")

   # Compute inverse map via preimages of algebra generators
   r = domain(F)
   s = codomain(F)
   n = ngens(r)
   m = ngens(s)

   (T, _, pr, _, D) = _surj_helper(F)
   psi = hom(s, r, [pr(D[i][2]) for i in 1:m])
   psi.kernel = ideal(s, [zero(s)])
   return psi
end

function preimage(F::AlgHom, f::Union{MPolyElem, MPolyQuoElem})
   @assert parent(f) == codomain(F)
   r = domain(F)
   s = codomain(F)
   n = ngens(r)
   m = ngens(s)

   (S, _, _, g) = _ring_helper(s, f, [zero(s)])
   (T, inc, pr, G, _) = _surj_helper(F)
   D = divrem(inc(g), G)
   !(leading_monomial(D[2]) < gen(T, m)) && error("Element not contained in image")
   return (pr(D[2]), kernel(F))
end

##############################################################################
#
# Normalization
#
##############################################################################

function _conv_normalize_alg(alg::Symbol)
  if alg == :primeDec
    return "prim"
  elseif alg == :equidimDec
    return "equidim"
  else
    error("algorithm invalid")
  end
end

function _conv_normalize_data(A, l, br)
  return [
    begin
      newSR = l[1][i][1]::Singular.PolyRing
      newOR, _ = PolynomialRing(br, [string(x) for x in gens(newSR)])
      newA, newAmap = quo(newOR, ideal(newOR, l[1][i][2][:norid]))
      newgens = newOR.(gens(l[1][i][2][:normap]))
      hom = AlgebraHomomorphism(A, newA, newA.(newgens))
      idgens = A.R.(gens(l[2][i]))
      (newA, hom, (A(idgens[end]), ideal(A, idgens)))
    end
    for i in 1:length(l[1])]
end

@doc Markdown.doc"""
    normalize(A::MPolyQuo; alg=:equidimDec)

Finds the normalization of a reduced affine algebra over a perfect field $K$:
Given the quotient $A=R/I$ of a multivariate polynomial ring $R$ over $K$
modulo a radical ideal $I$, compute the integral closure $\overline{A}$ 
of $A$ in its total ring of fractions $Q(A)$, together with the embedding 
$f: A \rightarrow \overline{A}$. The function relies on the algorithm 
of Greuel, Laplagne, and Seelisch which proceeds by finding a suitable decomposition 
$I=I_1\cap\dots\cap I_r$ into radical ideals $I_k$, together with
the normalization maps $f_k: R/I_k \rightarrow A_k=\overline{R/I_k}$, such that

$f=f_1\times \dots\times f_r: A \rightarrow A_1\times \dots\times A_r=\overline{A}$

is the normalization map of $A$. For each $k$, the function specifies two representations
of $A_k$: It returns an array of triples $(A_k, f_k, \mathfrak a_k)$,
where $A_k$ is represented as an affine $K$-algebra, and $f_k$ as a map of affine $K$-algebras.
The third entry $\mathfrak a_k$ is a tuple $(d_k, J_k)$, consisting of an element
$d_k\in A$ and an ideal $J_k\subset A$, such that $\frac{1}{d_k}J_k = A_k$ 
as $A$-submodules of the total ring of fractions of $A$.

By default, as a first step on its way to find the decomposition $I=I_1\cap\dots\cap I_r$, 
the algorithm computes an equidimensional decomposition of the radical ideal $I$.
Alternatively, if specified by `alg=:primeDec`, the algorithm computes $I=I_1\cap\dots\cap I_r$
as the prime decomposition of the radical ideal $I$.

CAVEAT: The function does not check whether $A$ is reduced. Use `isreduced(A)` in case 
you are unsure (this may take some time).
"""
function normalize(A::MPolyQuo; alg=:equidimDec)
  I = A.I
  br = base_ring(A.R)
  singular_assure(I)
  l = Singular.LibNormal.normal(I.gens.S, _conv_normalize_alg(alg))
  return _conv_normalize_data(A, l, br)
end

@doc Markdown.doc"""
    normalize_with_delta(A::MPolyQuo; alg=:equidimDec)

Finds the normalization

$f=f_1\times \dots\times f_r: A \rightarrow A_1\times \dots\times A_r=\overline{A}$

of $A$ as does `normalize(A)`, but computes additionally the delta invariant of $A$,
that is, the dimension $\dim_K(\overline{A}/A)$. More precisely, it returns a
tuple whose first element is `normalize(A)`, whose second element is an array
containing the delta invariants of the $A_k$, and whose third element is the
(total) delta invariant of $A$. The return value -1 in the third element
indicates that the delta invariant is infinite.
"""
function normalize_with_delta(A::MPolyQuo; alg=:equidimDec)
  I = A.I
  br = base_ring(A.R)
  singular_assure(I)
  l = Singular.LibNormal.normal(I.gens.S, _conv_normalize_alg(alg), "withDelta")
  return (_conv_normalize_data(A, l, br), l[3][1]::Vector{Int}, l[3][2]::Int)
end


##############################################################################
#
# Noether Normalization
#
##############################################################################

@doc Markdown.doc"""
    noether_normalization(A::MPolyQuo)

Given an affine algebra $A=R/I$ over a field $K$, the function returns  a triple $(L,f,g)$, such that:
$L$ is an array of $d=\dim A$ elements of $A$, all represented by linear forms in $R$, and
such that $K[L]\rightarrow A$ is a Noether normalization for $A$; $f: A \rightarrow A$ is an
automorphism of $A$, induced by a linear change of coordinates of $R$, and mapping the
$f_i$ to the last $d$ variables of $A$; and $g = f^{-1}$.

CAVEAT: The algorithm may not terminate over a small finite field. If it terminates,
the result is correct.
"""
function noether_normalization(A::MPolyQuo)
 I = A.I
 R = base_ring(I)
 singular_assure(I)
 l = Singular.LibAlgebra.noetherNormal(I.gens.S)
 i1 = [R(x) for x = gens(l[1])]
 i2 = [R(x) for x = gens(l[2])]
 m = matrix([[coeff(x, y) for y = gens(R)] for x = i1])
 mi = inv(m)
 h1 = AlgebraHomomorphism(A, A, map(A, i1))
 h2 = AlgebraHomomorphism(A, A, map(A, collect(matrix([gens(R)])'*map_entries(R, mi))))
 return map(x->h2(A(x)), i2), h1, h2
end


##############################################################################
#
# Integral bases
#
##############################################################################

@doc Markdown.doc"""
    integral_basis(f::MPolyElem, i::Int)

Given a polynomial $f$ in two variables with rational coefficients and an
integer $i\in\{1,2\}$ specifying one of the variables, $f$ must be irreducible
and monic in the specified variable: Say, $f\in\mathbb Q[x,y]$ is monic in $y$.
Then the normalization of $A = Q[x,y]/\langle f \rangle$, that is, the
integral closure $\overline{A}$ of $A$ in its quotient field, is a free
module over $K[x]$ of finite rank, and any set of free generators for
$\overline{A}$ over $K[x]$ is called an `integral basis` for $\overline{A}$
over $K[x]$. Relying on the algorithm by B\"ohm, Decker, Laplagne, and Pfister,
the function returns a pair $(d, V)$, where $d$ is an element of $A$,
and $V$ is a vector of elements in $A$, such that the fractions $v/d, v\in V$,
form an `integral basis` for $\overline{A}$
over $K[x]$.

NOTE: The conditions on $f$ are automatically checked.
"""
function integral_basis(f::MPolyElem, i::Int)
 aaa
end
