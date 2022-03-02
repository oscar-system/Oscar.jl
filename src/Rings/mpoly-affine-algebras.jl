
export normalization_with_delta
export noether_normalization, normalization, integral_basis
export isreduced, subalgebra_membership
export hilbert_series, hilbert_series_reduced, hilbert_series_expanded, hilbert_function, hilbert_polynomial, degree
export issurjective, isinjective, isbijective, inverse, preimage, isfinite
export _multi_hilbert_series, _multi_hilbert_series_reduced

##############################################################################
#
# Data associated to affine algebras
#
##############################################################################

@doc Markdown.doc"""
    dim(A::MPolyQuo)

Return the Krull dimension of `A`.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> A, _ = quo(R, ideal(R, [y-x^2, z-x^3]))
(Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(-x^2 + y, -x^3 + z), Map from
Multivariate Polynomial Ring in x, y, z over Rational Field to Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(-x^2 + y, -x^3 + z) defined by a julia-function with inverse)

julia> dim(A)
1
```
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
with denominator $q = (1-t^{w_1})\cdots (1-t^{w_n})$. 

See also `hilbert_series_reduced`.
"""
function hilbert_series(A:: MPolyQuo)
   if iszero(A.I)
      Zt, t = ZZ["t"]
      return (one(parent(t)), (1-t)^(ngens(A)))
   end
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
written in lowest terms. 

See also `hilbert_series`.
"""
function hilbert_series_reduced(A::MPolyQuo)
   if iszero(A.I)
      Zt, t = ZZ["t"]
      return (one(parent(t)), (1-t)^(ngens(A)))
   end
   H = HilbertData(A.I)
   return hilbert_series(H,2)
end

@doc Markdown.doc"""
    hilbert_series_expanded(A::MPolyQuo, d::Int)

Given a graded affine algebra $A = R/I$ over a field $K$ and an integer $d\geq 0$, return the
Hilbert series of $A$ to precision $d$. 
"""
function hilbert_series_expanded(A::MPolyQuo, d::Int)
   if iszero(A.I)
      P, t = PowerSeriesRing(QQ, d, "t")   
      b = zero(parent(t))
      n = ngens(A)
      for i=0:d
           b = b + binomial(n-1+i, i)*t^i
        end
      return b	   
     end
   H = HilbertData(A.I)  
   return hilbert_series_expanded(H, d)
end

@doc Markdown.doc"""
    hilbert_function(A::MPolyQuo, d::Int)

Given a graded affine algebra $A = R/I$ over a field $K$ and an integer $d\geq 0$, return the value
$H(A, d)$, where $H(A, \underline{\phantom{d}}): \N \rightarrow \N, d \mapsto \dim_K A_d$ is 
the Hilbert function of $A$.
"""
function hilbert_function(A::MPolyQuo, d::Int)
   if iszero(A.I)
       n = ngens(A)
       return binomial(n-1+d, n-1)
     end
   H = HilbertData(A.I)
   return hilbert_function(H, d)
end
   
@doc Markdown.doc"""
     hilbert_polynomial(A::MPolyQuo)

Given a graded affine algebra $A = R/I$ over a field $K$ such that the weights on the variables are all 1,
return the Hilbert polynomial of $A$.
"""
function hilbert_polynomial(A::MPolyQuo)::fmpq_poly
   if iszero(A.I)
       n = ngens(A)
       Qt, t = QQ["t"]
       b = one(parent(t))
       for i=1:(n-1)
           b = b * (t+i)
        end
       b = b//factorial(n-1)
       return b	   
     end
   H = HilbertData(A.I)
   return hilbert_polynomial(H)
end

@doc Markdown.doc"""
    degree(A::MPolyQuo)

Given a graded affine algebra $A = R/I$ over a field $K$ such that the weights on the variables are all 1,
return the degree of $A$.
"""
function degree(A::MPolyQuo)
   if iszero(A.I)
       return 1
     end
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

Given an affine algebra `A`, return `true` if `A` is reduced, `false` otherwise.

!!! warning
    The function computes the radical of the modulus of `A`. This may take some time.

# Examples
```jldoctest
julia> R, (x,) = PolynomialRing(QQ, ["x"])
(Multivariate Polynomial Ring in x over Rational Field, fmpq_mpoly[x])

julia> A, _ = quo(R, ideal(R, [x^4]))
(Quotient of Multivariate Polynomial Ring in x over Rational Field by ideal(x^4), Map from
Multivariate Polynomial Ring in x over Rational Field to Quotient of Multivariate Polynomial Ring in x over Rational Field by ideal(x^4) defined by a julia-function with inverse)

julia> isreduced(A)
false
```
"""
function isreduced(A::MPolyQuo) 
  I = A.I
  return I == radical(I)
end

@doc Markdown.doc"""
    isnormal(A::MPolyQuo)

Given an affine algebra `A` over a perfect field,
return `true` if `A` is normal, `false` otherwise.

!!! warning
    The function computes the normalization of `A`. This may take some time.
"""
function isnormal(A::MPolyQuo)
  _, _, d = normalization_with_delta(A)
  return d == 0
end

##############################################################################
#
# Algebra Containment
#
##############################################################################

# helper function for the containement problem, surjectivity and preimage
function _containement_helper(R::MPolyRing, N::Int, M::Int, I::MPolyIdeal, W::Vector, ord::Symbol)
   T, _ = PolynomialRing(base_ring(R), N + M, ordering = :lex) # :lex is needed in further computation,
                                                               # since we do not have block orderings in Oscar
   phi = hom(R, T, [gen(T, i) for i in 1:M])

   # Groebner computation
   A = phi.(gens(I))
   B = [phi(W[i])-gen(T, M + i) for i in 1:N]
   J = ideal(T, vcat(A, B))
   GG = gens(T)

   if ord == :lex
      #G = groebner_basis(J, complete_reduction = true)
      #return (T, phi, G)
      O = degrevlex(GG[1:M])*lex(GG[M+1:M+N])
   else ## ord == :degrevlex
      O = degrevlex(GG[1:M])*degrevlex(GG[M+1:M+N])
   end
   groebner_assure(J, O, complete_reduction = true)
   return (T, phi, J)
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
    subalgebra_membership(f::T, V::Vector{T}) where T <: Union{MPolyElem, MPolyQuoElem}

If  `f`  is contained in the subalgebra of `A` generated by the entries of `V`, return `(true, h)`, where `h` is giving the polynomial relation.
Otherwise, return `(false, 0)`.

# Examples
```jldoctest
julia> R, x = PolynomialRing(QQ, "x" => 1:3)
(Multivariate Polynomial Ring in x[1], x[2], x[3] over Rational Field, fmpq_mpoly[x[1], x[2], x[3]])

julia> f = x[1]^6*x[2]^6-x[1]^6*x[3]^6;

julia> V = [x[1]^3*x[2]^3-x[1]^3*x[3]^3, x[1]^3*x[2]^3+x[1]^3*x[3]^3]
2-element Vector{fmpq_mpoly}:
 x[1]^3*x[2]^3 - x[1]^3*x[3]^3
 x[1]^3*x[2]^3 + x[1]^3*x[3]^3

julia> subalgebra_membership(f, V)
(true, t_1*t_2)
```
"""
function subalgebra_membership(f::S, v::Vector{S}) where S <: Union{MPolyElem, MPolyQuoElem}
   @assert !isempty(v)
   r = parent(f)
   (R, I, W, F) = _ring_helper(r, f, v)
   n = length(W)
   m = ngens(R)

   # Build auxilliary objects
   (T, phi, J) = _containement_helper(R, n, m, I, W, :degrevlex)
   TT, _ = PolynomialRing(base_ring(T), ["t_$i" for i in 1:n], ordering = :lex)
   
   # Check containement
   D = normal_form([phi(F)], J)
   if leading_monomial(D[1]) < gen(T, m)
      A = [zero(TT) for i in 1:m]
      B = [gen(TT, i) for i in 1:n]
      psi = hom(T, TT, vcat(A, B))
      return (true, psi(D[1]))
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
# Stores the necessary data in groebner_data, resp. groebner_data_lex
function groebner_data(F::AlgHom, ord::Symbol)
   r = domain(F)
   s = codomain(F)
   n = ngens(r)
   m = ngens(s)

   if ord == :lex
         if !isdefined(F, :groebner_data_lex)
            (S, I, W, _) = _ring_helper(s, zero(s), F.image)
            # Build auxilliary objects
            (T, inc, J) = _containement_helper(S, n, m, I, W, ord)
            D = normal_form([gen(T, i) for i in 1:m], J)
            A = [zero(r) for i in 1:m]
            B = [gen(r, i) for i in 1:n]
            pr = hom(T, r, vcat(A, B))
            F.groebner_data_lex = (T, inc, pr, J, D)
          end
      return F.groebner_data_lex
   else ##ord == :degrevlex
         if !isdefined(F, :groebner_data)
            (S, I, W, _) = _ring_helper(s, zero(s), F.image)
            # Build auxilliary objects
            (T, inc, J) = _containement_helper(S, n, m, I, W, ord)
            D = normal_form([gen(T, i) for i in 1:m], J)
            A = [zero(r) for i in 1:m]
            B = [gen(r, i) for i in 1:n]
            pr = hom(T, r, vcat(A, B))
            F.groebner_data = (T, inc, pr, J, D)
         end
      return F.groebner_data
   end
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
   (T, _, _, _, D) = groebner_data(F, :degrevlex)

   # Check if map is surjective

   for i in 1:m
      if !(leading_monomial(D[i]) < gen(T, m))
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
  (T, _, _, J, _) = groebner_data(F, :lex)
  G = collect(J.gb)
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

   (T, _, pr, _, D) = groebner_data(F, :degrevlex)
   psi = hom(s, r, [pr(D[i]) for i in 1:m])
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
   (T, inc, pr, J, o) = groebner_data(F, :degrevlex)
   D = normal_form([inc(g)], J)
   !(leading_monomial(D[1]) < gen(T, m)) && error("Element not contained in image")
   return (pr(D[1]), kernel(F))
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
    normalization(A::MPolyQuo{<:MPolyElem{<:FieldElem}}; alg = :equidimDec)

Find the normalization of a reduced affine algebra over a perfect field $K$.
That is, given the quotient $A=R/I$ of a multivariate polynomial ring $R$ over $K$
modulo a radical ideal $I$, compute the integral closure $\overline{A}$ 
of $A$ in its total ring of fractions $Q(A)$, together with the embedding 
$f: A \rightarrow \overline{A}$. 

# Implemented Algorithms and how to Read the Output

The function relies on the algorithm 
of Greuel, Laplagne, and Seelisch which proceeds by finding a suitable decomposition 
$I=I_1\cap\dots\cap I_r$ into radical ideals $I_k$, together with
maps $A = R/I \rightarrow A_k=\overline{R/I_k}$ which give rise to the normalization map of $A$:

$A\hookrightarrow A_1\times \dots\times A_r=\overline{A}$

For each $k$, the function specifies two representations
of $A_k$: It returns an array of triples $(A_k, f_k, \mathfrak a_k)$,
where $A_k$ is represented as an affine $K$-algebra, and $f_k$ as a map of affine $K$-algebras.
The third entry $\mathfrak a_k$ is a tuple $(d_k, J_k)$, consisting of an element
$d_k\in A$ and an ideal $J_k\subset A$, such that $\frac{1}{d_k}J_k = A_k$ 
as $A$-submodules of the total ring of fractions of $A$.

By default (`alg = :equidimDec`), as a first step on its way to find the decomposition $I=I_1\cap\dots\cap I_r$, 
the algorithm computes an equidimensional decomposition of the radical ideal $I$.
Alternatively, if specified by `alg = :primeDec`, the algorithm computes $I=I_1\cap\dots\cap I_r$
as the prime decomposition of the radical ideal $I$.

See [GLS10](@cite).

!!! warning
    The function does not check whether $A$ is reduced. Use `isreduced(A)` in case you are unsure (this may take some time).

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"]);

julia> A, _ = quo(R, ideal(R, [(x^2-y^3)*(x^2+y^2)*x]));

julia> L = normalization(A);

julia> size(L)
(2,)

julia> LL = normalization(A, alg = :primeDec);

julia> size(LL)
(3,)

julia> LL[1][1]
Quotient of Multivariate Polynomial Ring in T(1), x, y over Rational Field by ideal(-T(1)*y + x, -T(1)*x + y^2, T(1)^2 - y, -x^2 + y^3)

julia> LL[1][2]
Algebra homomorphism with

domain: Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x^5 - x^3*y^3 + x^3*y^2 - x*y^5)

codomain: Quotient of Multivariate Polynomial Ring in T(1), x, y over Rational Field by ideal(-T(1)*y + x, -T(1)*x + y^2, T(1)^2 - y, -x^2 + y^3)

defining images of generators: MPolyQuoElem{fmpq_mpoly}[x, y]

julia> LL[1][3]
(y, ideal(x, y))
```
"""
function normalization(A::MPolyQuo{<:MPolyElem{<:FieldElem}}; alg=:equidimDec)
  I = A.I
  br = base_ring(A.R)
  singular_assure(I)
  l = Singular.LibNormal.normal(I.gens.S, _conv_normalize_alg(alg))
  return _conv_normalize_data(A, l, br)
end

@doc Markdown.doc"""
    normalization_with_delta(A::MPolyQuo{<:MPolyElem{<:FieldElem}}; alg = :equidimDec)

Compute the normalization

$A\hookrightarrow A_1\times \dots\times A_r=\overline{A}$

of $A$ as does `normalize(A)`, but return additionally the `delta invariant` of $A$,
that is, the dimension 

$\dim_K(\overline{A}/A)$. 

# How to Read the Output

The return value is a tuple whose first element is `normalize(A)`, whose second element is an array
containing the delta invariants of the $A_k$, and whose third element is the
(total) delta invariant of $A$. The return value -1 in the third element
indicates that the delta invariant is infinite.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"]);

julia> A, _ = quo(R, ideal(R, [(x^2-y^3)*(x^2+y^2)*x]));

julia> L = normalization_with_delta(A);

julia> L[2]
3-element Vector{Int64}:
 1
 1
 0

julia> L[3]
13

julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> A, _ = quo(R, ideal(R, [z^3-x*y^4]))
(Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(-x*y^4 + z^3), Map from
Multivariate Polynomial Ring in x, y, z over Rational Field to Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(-x*y^4 + z^3) defined by a julia-function with inverse)

julia> L = normalization_with_delta(A)
(Tuple{MPolyQuo{fmpq_mpoly}, AlgHom{fmpq}, Tuple{MPolyQuoElem{fmpq_mpoly}, MPolyQuoIdeal{fmpq_mpoly}}}[(Quotient of Multivariate Polynomial Ring in T(1), T(2), x, y, z over Rational Field by ideal(T(1)*y - T(2)*z, T(2)*y - z, -T(1)*z + x*y^2, T(1)^2 - x*z, T(1)*T(2) - x*y, -T(1) + T(2)^2, x*y^4 - z^3), Algebra homomorphism with

domain: Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(-x*y^4 + z^3)

codomain: Quotient of Multivariate Polynomial Ring in T(1), T(2), x, y, z over Rational Field by ideal(T(1)*y - T(2)*z, T(2)*y - z, -T(1)*z + x*y^2, T(1)^2 - x*z, T(1)*T(2) - x*y, -T(1) + T(2)^2, x*y^4 - z^3)

defining images of generators: MPolyQuoElem{fmpq_mpoly}[x, y, z]
, (z^2, ideal(x*y^2*z, x*y^3, z^2)))], [-1], -1)
```
"""
function normalization_with_delta(A::MPolyQuo{<:MPolyElem{<:FieldElem}}; alg=:equidimDec)
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
    noether_normalization(A::MPolyQuo{<:MPolyElem{<:FieldElem}})

Given an affine algebra $A=R/I$ over a field $K$, return a triple $(V,F,G)$ such that:
$V$ is a vector of $d=\dim A$ elements of $A$, represented by linear forms $l_i\in R$, and
such that $K[V]\hookrightarrow A$ is a Noether normalization for $A$; $F: A=R/I \rightarrow B = R/\phi(I)$ 
is an isomorphism, induced by a linear change $ \phi $ of coordinates of $R$ which maps the
$l_i$ to the the last $d$ variables of $R$; and $G = F^{-1}$.

!!! warning
    The algorithm may not terminate over a small finite field. If it terminates, the result is correct.

"""
function noether_normalization(A::MPolyQuo{<:MPolyElem{<:FieldElem}})
 I = A.I
 R = base_ring(I)
 singular_assure(I)
 l = Singular.LibAlgebra.noetherNormal(I.gens.S)
 i1 = [R(x) for x = gens(l[1])]
 i2 = [R(x) for x = gens(l[2])]
 m = matrix([[coeff(x, y) for y = gens(R)] for x = i1])
 mi = inv(m)
 ###mi_arr = [collect(matrix([gens(R)])*map_entries(R, mi))[i] for i in 1:ngens(R)]
 mi_arr = [collect(map_entries(R, mi)*gens(R))[i] for i in 1:ngens(R)]
 h = AlgebraHomomorphism(R, R, i1)
 V = map(x->h(x), gens(I))
 B, _ = quo(R, ideal(R, V))
 h1 = AlgebraHomomorphism(A, B, map(B, i1))
 h2 = AlgebraHomomorphism(B, A, map(A, mi_arr))
 return map(x->h2(B(x)), i2), h1, h2
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
$\overline{A}$ over $K[x]$ is called an *integral basis* for $\overline{A}$
over $K[x]$. Relying on the algorithm by [BDLP19](@cite),
the function returns a pair $(d, V)$, where $d$ is an element of $A$,
and $V$ is a vector of elements in $A$, such that the fractions $v/d, v\in V$,
form an integral basis for $\overline{A}$ over $K[x]$. 

!!! note
    The conditions on $f$ are automatically checked.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> f = (y^2-2)^2 + x^5
x^5 + y^4 - 4*y^2 + 4

julia> integral_basis(f, 2)
(x^2, MPolyQuoElem{fmpq_mpoly}[x^2, x^2*y, y^2 - 2, y^3 - 2*y])
```
"""
function integral_basis(f::MPolyElem, i::Int)
  R = parent(f)

  if typeof(R) <: MPolyRing_dec
    throw(ArgumentError("Not implemented for decorated rings."))
  end
  
  if !(nvars(R) == 2)
    throw(ArgumentError("The parent ring must be a polynomial ring in two variables."))
  end

  if !(i == 1 || i == 2)
    throw(ArgumentError("The index $i must be either 1 or 2, indicating the integral variable."))
  end

  if !(coefficient_ring(R) == QQ || base_ring(R) == Singular.QQ)
    throw(ArgumentError("The coefficient ring must be the rationals."))
  end

  if !isone(coeff(f, [i], [degree(f, i)]))
    throw(ArgumentError("The input polynomial must be monic as a polynomial in $(gen(R,i))"))
  end

  if !isirreducible(f)
    throw(ArgumentError("The input polynomial must be irreducible"))
  end

  SR = singular_ring(R)
  l = Singular.LibIntegralbasis.integralBasis(SR(f), i, "isIrred")
  A, p = quo(R, ideal(R, [f]))
  ###return (R(l[2]), R.(gens(l[1])))
  return (p(R(l[2])), [p(R(x)) for x = gens(l[1])])
end

