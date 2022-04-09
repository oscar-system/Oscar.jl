
export normalization_with_delta
export noether_normalization, normalization, integral_basis
export isreduced, subalgebra_membership, minimal_subalgebra_generators
export hilbert_series, hilbert_series_reduced, hilbert_series_expanded, hilbert_function, hilbert_polynomial, degree
export issurjective, isinjective, isbijective, inverse, preimage, isfinite
export multi_hilbert_series, multi_hilbert_series_reduced, multi_hilbert_function

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


##################################################################################
###    z-graded Hilbert series stuff using Singular for finding the Hilbert series
###    from mpoly-graded.jl
##################################################################################


@doc Markdown.doc"""
    hilbert_series(A::MPolyQuo)

Given a $\mathbb Z$-graded affine algebra $A = R/I$ over a field $K$, where the grading 
is inherited from a $\mathbb Z$-grading on the polynomial ring $R$ defined by assigning 
positive integer weights to the variables, and where $I$ is a homogeneous ideal of $R$ 
with respect to this grading, return a pair $(p,q)$, say, of univariate polynomials 
$p, q\in\mathbb Z[t]$ such that $p/q$ represents the Hilbert series of $A$ as 
a rational function with denominator 

$q = (1-t^{w_1})\cdots (1-t^{w_n}),$

where $n$ is the number of variables of $R$, and $w_1, \dots, w_n$ are the assigned weights.

See also `hilbert_series_reduced`.

# Examples
```jldoctest
julia> R, (w, x, y, z) = GradedPolynomialRing(QQ, ["w", "x", "y", "z"]);

julia> A, _ = quo(R, ideal(R, [w*y-x^2, w*z-x*y, x*z-y^2]));

julia> hilbert_series(A)
(2*t^3 - 3*t^2 + 1, t^4 - 4*t^3 + 6*t^2 - 4*t + 1)

julia> R, (x, y, z) = GradedPolynomialRing(QQ, ["x", "y", "z"], [1, 2, 3]);

julia> A, _ = quo(R, ideal(R, [x*y*z]));

julia> hilbert_series(A)
(-t^6 + 1, -t^6 + t^5 + t^4 - t^2 - t + 1)
```
"""
function hilbert_series(A:: MPolyQuo)
   if iszero(A.I)
      R = base_ring(A.I)
      W = R.d
      W = [Int(W[i][1]) for i = 1:ngens(R)]   
      Zt, t = ZZ["t"]
      den = prod([1-t^W[i] for i = 1:ngens(base_ring(A.I))])
      return (one(parent(t)), den)
   end
   H = HilbertData(A.I)
   return hilbert_series(H,1)
end


@doc Markdown.doc"""
    hilbert_series_reduced(A::MPolyQuo)

Given a $\mathbb Z$-graded affine algebra $A = R/I$ over a field $K$, where the grading 
is inherited from a $\mathbb Z$-grading on the polynomial ring $R$ defined by assigning 
positive integer weights to the variables, and where $I$ is a homogeneous ideal of $R$ 
with respect to this grading, return a pair $(p,q)$, say, of univariate polynomials 
$p, q\in\mathbb Z[t]$ such that $p/q$ represents the Hilbert series of $A$ as 
a rational function written in lowest terms. 

See also `hilbert_series`.

# Examples
```jldoctest
julia> R, (w, x, y, z) = GradedPolynomialRing(QQ, ["w", "x", "y", "z"]);

julia> A, _ = quo(R, ideal(R, [w*y-x^2, w*z-x*y, x*z-y^2]));

julia> hilbert_series_reduced(A)
(2*t + 1, t^2 - 2*t + 1)

julia> R, (x, y, z) = GradedPolynomialRing(QQ, ["x", "y", "z"], [1, 2, 3]);

julia> A, _ = quo(R, ideal(R, [x*y*z]));

julia> hilbert_series(A)
(-t^6 + 1, -t^6 + t^5 + t^4 - t^2 - t + 1)

julia> hilbert_series_reduced(A)
(t^2 - t + 1, t^2 - 2*t + 1)
```
"""
function hilbert_series_reduced(A::MPolyQuo)
   if iszero(A.I)
      return hilbert_series(A)
   end
   H = HilbertData(A.I)
   return hilbert_series(H,2)
end

@doc Markdown.doc"""
    hilbert_series_expanded(A::MPolyQuo, d::Int)

Given a $\mathbb Z$-graded affine algebra $A = R/I$ over a field $K$, where the grading 
is inherited from a $\mathbb Z$-grading on the polynomial ring $R$ defined by assigning 
positive integer weights to the variables, and where $I$ is a homogeneous ideal of $R$ 
with respect to this grading, return the Hilbert series of $A$ to precision $d$. 

# Examples
```jldoctest
julia> R, (w, x, y, z) = GradedPolynomialRing(QQ, ["w", "x", "y", "z"]);

julia> A, _ = quo(R, ideal(R, [w*y-x^2, w*z-x*y, x*z-y^2]));

julia> hilbert_series_expanded(A, 7)
1 + 4*t + 7*t^2 + 10*t^3 + 13*t^4 + 16*t^5 + 19*t^6 + 22*t^7 + O(t^8)

julia> R, (x, y, z) = GradedPolynomialRing(QQ, ["x", "y", "z"], [1, 2, 3]);

julia> A, _ = quo(R, ideal(R, [x*y*z]));

julia> hilbert_series_expanded(A, 5)
1 + t + 2*t^2 + 3*t^3 + 4*t^4 + 5*t^5 + O(t^6)
```
"""
function hilbert_series_expanded(A::MPolyQuo, d::Int)
   H = HilbertData(A.I)  
   return hilbert_series_expanded(H, d)
end

@doc Markdown.doc"""
    hilbert_function(A::MPolyQuo, d::Int)

Given a $\mathbb Z$-graded affine algebra $A = R/I$ over a field $K$, where the grading 
is inherited from a $\mathbb Z$-grading on the polynomial ring $R$ defined by assigning 
positive integer weights to the variables, and where $I$ is a homogeneous ideal of $R$ 
with respect to this grading, return the value $H(A, d)$, where 

$H(A, \underline{\phantom{d}}): \N \rightarrow \N, \; d  \mapsto \dim_K A_d,$ 

is the Hilbert function of $A$.

# Examples
```jldoctest
julia> R, (w, x, y, z) = GradedPolynomialRing(QQ, ["w", "x", "y", "z"]);

julia> A, _ = quo(R, ideal(R, [w*y-x^2, w*z-x*y, x*z-y^2]));

julia> hilbert_function(A,7)
22

julia> R, (x, y, z) = GradedPolynomialRing(QQ, ["x", "y", "z"], [1, 2, 3]);

julia> A, _ = quo(R, ideal(R, [x*y*z]));

julia> hilbert_function(A, 5)
5
```
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

Given a $\mathbb Z$-graded affine algebra $A = R/I$ over a field $K$, where the grading 
is inherited from the standard $\mathbb Z$-grading on the polynomial ring $R$,
return the Hilbert polynomial of $A$.

# Examples
```jldoctest
julia> R, (w, x, y, z) = GradedPolynomialRing(QQ, ["w", "x", "y", "z"]);

julia> A, _ = quo(R, ideal(R, [w*y-x^2, w*z-x*y, x*z-y^2]));

julia> hilbert_polynomial(A)
3*t + 1
```
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

Given a $\mathbb Z$-graded affine algebra $A = R/I$ over a field $K$, where the grading 
is inherited from the standard $\mathbb Z$-grading on the polynomial ring $R$,
return the degree of $A$.

# Examples
```jldoctest
julia> R, (w, x, y, z) = GradedPolynomialRing(QQ, ["w", "x", "y", "z"]);

julia> A, _ = quo(R, ideal(R, [w*y-x^2, w*z-x*y, x*z-y^2]));

julia> degree(A)
3
```
"""
function degree(A::MPolyQuo)
   if iszero(A.I)
       return 1
     end
   H = HilbertData(A.I)
   return degree(H)
end

###############################################################################
### General Hilbert series stuff using Singular for computing ideal quotients
###############################################################################

function transform_to_positive_orthant(rs::Matrix{Int})   
    C = positive_hull(rs)
    @assert isfulldimensional(C) "Cone spanned by generator degrees needs to be full-dimensional"
    F = linear_inequality_matrix(facets(C))
    
    # Find a simplicial cone containing C
    index = 2
    full_rank_subset = [1]
    full_rank = rank(F)
    current_rank = rank(F[full_rank_subset, :])
    while current_rank < full_rank
        for i in index:nrows(F)
            test = Vector{Int}(full_rank_subset)
            append!(test, i)
            testrank = rank(F[test,:])
            if testrank > current_rank
                index = i+1
                current_rank = testrank
                append!(full_rank_subset, i)
                break
            end
        end
    end
    Csimplicial = cone_from_inequalities(F[full_rank_subset,:])
    
    @assert Polymake.polytope.included_polyhedra(C.pm_cone, Csimplicial.pm_cone) "Cone containment violated"
    CsRays = Polymake.common.primitive(Csimplicial.pm_cone.RAYS)
    CsRays = matrix(ZZ, CsRays)
    nf = AbstractAlgebra.hnf_with_transform(transpose(CsRays))
    CsRays_transformed = transpose(nf[1])
    transformation = transpose(nf[2])
    @assert CsRays * transformation == CsRays_transformed "Maybe order of transformation is wrong?"
    original = matrix(ZZ, rs)
    return original * transformation, transformation
end

function _numerator_monomial_multi_hilbert_series(I::MPolyIdeal, S)
   ###for use in _multi_hilbert_series only
   ###if !ismonomial(I)
   ###      throw(ArgumentError("The ideal is not monomial"))
   ###end
   ### V = minimal_monomial_generators(I)  ### to be written
   V = gens(I)
   s = ngens(I)
   d = degree(Vector{Int}, V[s])
   B = MPolyBuildCtx(S)
   push_term!(B, 1, d)
   p = finish(B)
   if s == 1
      return 1-p
   end
   v = V[s]
   V = deleteat!(V, s)
   J = ideal(base_ring(I), V)
   p1 = _numerator_monomial_multi_hilbert_series(J, S)
   p2 = _numerator_monomial_multi_hilbert_series(J:v, S)
   ### TODO: Do I have the minimal set of monomial generators here?
   return p1-p*p2
end

@doc Markdown.doc"""
    multi_hilbert_series(A::MPolyQuo)

Return the Hilbert series of the positively graded affine algebra `A`.

# Examples
```jldoctest
julia> W = [1 1 1; 0 0 -1];

julia> R, x = GradedPolynomialRing(QQ, ["x[1]", "x[2]", "x[3]"], W)
(Multivariate Polynomial Ring in x[1], x[2], x[3] over Rational Field graded by
  x[1] -> [1 0]
  x[2] -> [1 0]
  x[3] -> [1 -1], MPolyElem_dec{fmpq, fmpq_mpoly}[x[1], x[2], x[3]])

julia> I = ideal(R, [x[1]^3*x[2], x[2]*x[3]^2, x[2]^2*x[3], x[3]^4]);

julia> A, _ = quo(R, I);

julia> H = multi_hilbert_series(A);

julia> H[1][1]
-t[1]^5*t[2]^2 + t[1]^5*t[2] + t[1]^4*t[2]^2 - t[1]^4 + t[1]^2*t[2]^2 - t[1]^2*t[2] + t[1]*t[2]^4 - t[1]*t[2]^2 - t[2]^4 + 1

julia> H[1][2]
-t[1]^2*t[2] + t[1]^2 + 2*t[1]*t[2] - 2*t[1] - t[2] + 1

julia> H[2]
[1    0]
[1   -1]

julia> G = abelian_group(fmpz_mat([1 -1]));

julia> g = gen(G, 1)
Element of
(General) abelian group with relation matrix
[1 -1]
with components [0 1]

julia> W = [g, g, g, g];

julia> R, (w, x, y, z) = GradedPolynomialRing(QQ, ["w", "x", "y", "z"], W);

julia> A, _ = quo(R, ideal(R, [w*y-x^2, w*z-x*y, x*z-y^2]));

julia> multi_hilbert_series(A)
((2*t^3 - 3*t^2 + 1, t^4 - 4*t^3 + 6*t^2 - 4*t + 1), [1])
```
"""
function multi_hilbert_series(A::MPolyQuo)
   R = A.R
   I = A.I
   if !(coefficient_ring(R) isa AbstractAlgebra.Field)
       throw(ArgumentError("The coefficient ring must be a field."))
   end
   if !(R isa MPolyRing_dec && isgraded(R) && is_positively_graded(R))
       throw(ArgumentError("The base ring must be positively graded."))
   end
   if !(is_zm_graded(R))
      G = grading_group(R)
      H, iso = snf(G)
      V = [preimage(iso, x) for x in gens(G)]
      isoinv = hom(G, H, V)
      W = R.d
      W = [isoinv(W[i]) for i = 1:length(W)]
      S, _ = GradedPolynomialRing(coefficient_ring(R), map(string, symbols(R)), W)
      change = hom(R, S, gens(S))
      I = change(A.I)
      R = S
   end
   m = ngens(grading_group(R))  
   n = ngens(R)
   W = R.d
   MI = Matrix{Int}(undef, n, m)
   for i=1:n
       for j=1:m
           MI[i, j] = Int(W[i][j])
       end
   end
   minMI = minimum(MI)
   if minMI<0
      MI, T = transform_to_positive_orthant(MI)     
   else
      T = identity_matrix(ZZ, m)
   end  
   if m == 1
      VAR = ["t"]
   else
      VAR = [_make_variable("t", i) for i = 1:m]
   end
   S, _ = PolynomialRing(ZZ, VAR) 
   q = one(S)
   for i = 1:n
      e = [Int(MI[i, :][j]) for j = 1:m]
      B = MPolyBuildCtx(S)
      push_term!(B, 1, e)
      q = q*(1-finish(B))
   end
   if iszero(I)
      p = one(S)
   else
      LI = leading_ideal(I, ordering=degrevlex(gens(R)))
      if minMI<0
         RNEW, _ = GradedPolynomialRing(coefficient_ring(R), [String(symbols(R)[i]) for i = 1:n], Matrix(transpose(MI)))
         LI = ideal(RNEW, [RNEW(LI[i]) for i = 1:ngens(LI)])
      end
      p = _numerator_monomial_multi_hilbert_series(LI, S)
   end
   return  (p, q), T
end

@doc Markdown.doc"""
    multi_hilbert_series_reduced(A::MPolyQuo)

Return the reduced Hilbert series of the positively graded affine algebra `A`.

# Examples
```jldoctest
julia> W = [1 1 1; 0 0 -1];

julia> R, x = GradedPolynomialRing(QQ, ["x[1]", "x[2]", "x[3]"], W)
(Multivariate Polynomial Ring in x[1], x[2], x[3] over Rational Field graded by
  x[1] -> [1 0]
  x[2] -> [1 0]
  x[3] -> [1 -1], MPolyElem_dec{fmpq, fmpq_mpoly}[x[1], x[2], x[3]])

julia> I = ideal(R, [x[1]^3*x[2], x[2]*x[3]^2, x[2]^2*x[3], x[3]^4]);

julia> A, _ = quo(R, I);

julia> H = multi_hilbert_series_reduced(A);

julia> H[1][1]
-t[1]^4*t[2] + t[1]^3 + t[1]^2 + t[1]*t[2] + t[1] + t[2]^3 + t[2]^2 + t[2] + 1

julia> H[1][2]
-t[1] + 1

julia> H[2]
[1    0]
[1   -1]

julia> G = abelian_group(fmpz_mat([1 -1]));

julia> g = gen(G, 1);

julia> W = [g, g, g, g];

julia> R, (w, x, y, z) = GradedPolynomialRing(QQ, ["w", "x", "y", "z"], W);

julia> A, _ = quo(R, ideal(R, [w*y-x^2, w*z-x*y, x*z-y^2]));

julia> multi_hilbert_series_reduced(A)
((2*t + 1, t^2 - 2*t + 1), [1])
```
"""
function multi_hilbert_series_reduced(A::MPolyQuo)
   (p, q), T = multi_hilbert_series(A::MPolyQuo)
   c = gcd(p, q)
   p = divexact(p, c)
   q = divexact(q, c)
   return (constant_coefficient(q)*p, constant_coefficient(q)*q), T
end

function _monomial_ideal_membership(m::MPolyElem, I::MPolyIdeal)
   ### for potential use in multi_hilbert_function
   ### I is supposed to be given by monomial generators, ordered by
   ### increasing (total) degree, m is supposed to be a monomial
        for i = 1:ngens(I)
	     if total_degree(gen(I, i))>total_degree(m)
	        break
	     end
	     if minimum(exponent_vector(m, 1)-exponent_vector(gen(I, i), 1))>=0
	         return true
             end
        end
	return false
end

@doc Markdown.doc"""
    multi_hilbert_function(A::MPolyQuo, g::GrpAbFinGenElem)

Given an element $g$ of a finitely generated Abelian group $G$, and given
an affine algebra $A = R/I$ over a field $K$, where $R$ is positively graded
by $G$, and where $I$ is homogeneous with respect to this grading, consider 
the induced grading on $A$, and return the value $H(A, g)$ of 
the Hilbert function

$H(A, \underline{\phantom{d}}): G \to \N, \; g\mapsto \dim_K(A_g).$

    multi_hilbert_function(A::MPolyQuo, g::Vector{<:IntegerUnion})

Given an affine algebra $A = R/I$ as above, where $G =\mathbb  Z^m$, 
and given a vector $g$ of $m$ integers, convert $g$ into an element 
of the group $\mathbb  Z^m$, and return the value $H(A, g)$
as above.

    multi_hilbert_function(A::MPolyQuo, g::IntegerUnion)

Given an affine algebra $A = R/I$ as above, where $G =\mathbb  Z$,
and given an integer $g$, convert $g$ into an element of the group 
$\mathbb  Z$, and return the value $H(A, g)$ as above.

# Examples
```jldoctest
julia> W = [1 1 1; 0 0 -1];

julia> R, x = GradedPolynomialRing(QQ, ["x[1]", "x[2]", "x[3]"], W)
(Multivariate Polynomial Ring in x[1], x[2], x[3] over Rational Field graded by 
  x[1] -> [1 0]
  x[2] -> [1 0]
  x[3] -> [1 -1], MPolyElem_dec{fmpq, fmpq_mpoly}[x[1], x[2], x[3]])

julia> I = ideal(R, [x[1]^3*x[2], x[2]*x[3]^2, x[2]^2*x[3], x[3]^4]);

julia> A, _ = quo(R, I);

julia> multi_hilbert_function(A::MPolyQuo, [1, 0])
2

julia> R, (w, x, y, z) = GradedPolynomialRing(QQ, ["w", "x", "y", "z"], [-1, -1, -1, -1]);

julia> A, _ = quo(R, ideal(R, [w*y-x^2, w*z-x*y, x*z-y^2]));

julia> multi_hilbert_function(A, -7)
22

julia> G = abelian_group(fmpz_mat([1 -1]));

julia> g = gen(G, 1);

julia> W = [g, g, g, g];

julia> R, (w, x, y, z) = GradedPolynomialRing(QQ, ["w", "x", "y", "z"], W);

julia> A, _ = quo(R, ideal(R, [w*y-x^2, w*z-x*y, x*z-y^2]));

julia> multi_hilbert_function(A, 7*g)
22
```
"""
function multi_hilbert_function(A::MPolyQuo, g::GrpAbFinGenElem)
    R = A.R
    if !(coefficient_ring(R) isa AbstractAlgebra.Field)
       throw(ArgumentError("The coefficient ring of must be a field."))
    end
    LI = leading_ideal(A.I, ordering=degrevlex(gens(R)))
    ### TODO: Decide whether we should check whether a GB with respect
    ### to another degree-compatible ordering is already available
    L = homogeneous_component(R, g);
    if rank(L[1]) == 0
       return 0
    end
    FG = gens(L[1]);
    EMB = L[2]
    cc = 0
    for i in 1:length(FG)
         if !(_monomial_ideal_membership(EMB(FG[i]), LI))
	 ### if !(EMB(FG[i]) in LI)
	    cc = cc +1
         end
    end
    return cc
end

function multi_hilbert_function(A::MPolyQuo, g::Vector{<:IntegerUnion})
  R = A.R
  @assert is_zm_graded(R)
  return multi_hilbert_function(A, grading_group(R)(g))
end

function multi_hilbert_function(A::MPolyQuo, g::IntegerUnion)
  R = A.R
  @assert is_z_graded(R)
  return multi_hilbert_function(A, grading_group(R)([g]))
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
   groebner_assure(J, O, true)
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
 
Given an element `f` of a graded multivariate polynomial ring over a field, or of a quotient ring of such a ring,
and given a vector `V` of elements in the same ring, consider the subalgebra generated by the entries of `V` in
that ring. If `f` is contained in the subalgebra, return `(true, h)`, where `h` is giving the polynomial relation.
Return, `(false, 0)`, otherwise.

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
   r = parent(f)
   if !(coefficient_ring(r) isa AbstractAlgebra.Field)
       throw(ArgumentError("The coefficient ring must be a field."))
   end
   @assert !isempty(v)
   @assert all(x->parent(x) == r, v)
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

################################################################################
#
#  Minimalizing a set of subalgebra generators in graded case
#
################################################################################
@doc Markdown.doc""" 
    minimal_subalgebra_generators(V::Vector{T}) where T <: Union{MPolyElem, MPolyQuoElem}

Given a vector `V` of homogeneous elements of a positively graded multivariate 
polynomial ring over a field, or of a quotient ring of such a ring, return a minimal
subset of the elements in `V` which, in the given ring, generate the same subalgebra
as all elements in `V`.

!!! note
    The conditions on `V` and the given ring are automatically checked.

# Examples
```jldoctest
julia> R, (x, y) = GradedPolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field graded by 
  x -> [1]
  y -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y])

julia> V = [x, y, x^2+y^2]
3-element Vector{MPolyElem_dec{fmpq, fmpq_mpoly}}:
 x
 y
 x^2 + y^2

julia> minimal_subalgebra_generators(V)
2-element Vector{MPolyElem_dec{fmpq, fmpq_mpoly}}:
 x
 y
```
"""
function minimal_subalgebra_generators(V::Vector{T}) where T <: Union{MPolyElem, MPolyQuoElem}
  p = parent(V[1])
  @assert all(x->parent(x) == p, V)
  if !(coefficient_ring(p) isa AbstractAlgebra.Field)
     throw(ArgumentError("The coefficient ring must be a field."))
  end
  if p isa MPolyRing
      p isa MPolyRing_dec && is_positively_graded(p) || throw(ArgumentError("The base ring must be positively graded"))
  else
      p.R isa MPolyRing_dec && isgraded(p.R) || throw(ArgumentError("The base ring must be graded"))
  end
  all(ishomogeneous, V) || throw(ArgumentError("The input data is not homogeneous"))
  # iterate over the generators, starting with those in lowest degree, then work up
  W = sort(V, by = x -> degree(x)[1])
  result = [ W[1] ]
  for elm in W
    if !subalgebra_membership(elm, result)[1]
       push!(result, elm)
    end
  end
  return result
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
      _hom = hom(A, newA, newA.(newgens))
      idgens = A.R.(gens(l[2][i]))
      (newA, _hom, (A(idgens[end]), ideal(A, idgens)))
    end
    for i in 1:length(l[1])]
end

@doc Markdown.doc"""
    normalization(A::MPolyQuo; alg = :equidimDec)

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
Map with following data
Domain:
=======
Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x^5 - x^3*y^3 + x^3*y^2 - x*y^5)
Codomain:
=========
Quotient of Multivariate Polynomial Ring in T(1), x, y over Rational Field by ideal(-T(1)*y + x, -T(1)*x + y^2, T(1)^2 - y, -x^2 + y^3)

julia> LL[1][3]
(y, ideal(x, y))
```
"""
function normalization(A::MPolyQuo; alg=:equidimDec)
  if !(coefficient_ring(A) isa AbstractAlgebra.Field)
       throw(ArgumentError("The coefficient ring must be a field."))
  end
  if A.R isa MPolyRing_dec
    throw(ArgumentError("Not implemented for quotients of decorated rings."))
  end
  I = A.I
  br = base_ring(A.R)
  singular_assure(I)
  l = Singular.LibNormal.normal(I.gens.S, _conv_normalize_alg(alg))
  return _conv_normalize_data(A, l, br)
end

@doc Markdown.doc"""
    normalization_with_delta(A::MPolyQuo; alg = :equidimDec)

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
(Tuple{MPolyQuo{fmpq_mpoly}, Oscar.MPolyAnyMap{MPolyQuo{fmpq_mpoly}, MPolyQuo{fmpq_mpoly}, Nothing, MPolyQuoElem{fmpq_mpoly}}, Tuple{MPolyQuoElem{fmpq_mpoly}, MPolyQuoIdeal{fmpq_mpoly}}}[(Quotient of Multivariate Polynomial Ring in T(1), T(2), x, y, z over Rational Field by ideal(T(1)*y - T(2)*z, T(2)*y - z, -T(1)*z + x*y^2, T(1)^2 - x*z, T(1)*T(2) - x*y, -T(1) + T(2)^2, x*y^4 - z^3), Map with following data
Domain:
=======
Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(-x*y^4 + z^3)
Codomain:
=========
Quotient of Multivariate Polynomial Ring in T(1), T(2), x, y, z over Rational Field by ideal(T(1)*y - T(2)*z, T(2)*y - z, -T(1)*z + x*y^2, T(1)^2 - x*z, T(1)*T(2) - x*y, -T(1) + T(2)^2, x*y^4 - z^3), (z^2, ideal(x*y^2*z, x*y^3, z^2)))], [-1], -1)
```
"""
function normalization_with_delta(A::MPolyQuo; alg=:equidimDec)
  if !(coefficient_ring(A) isa AbstractAlgebra.Field)
       throw(ArgumentError("The coefficient ring must be a field."))
  end
  if A.R isa MPolyRing_dec
    throw(ArgumentError("Not implemented for quotients of decorated rings."))
  end
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

Given an affine algebra $A=R/I$ over a field $K$, return a triple $(V,F,G)$ such that:
$V$ is a vector of $d=\dim A$ elements of $A$, represented by linear forms $l_i\in R$, and
such that $K[V]\hookrightarrow A$ is a Noether normalization for $A$; $F: A=R/I \rightarrow B = R/\phi(I)$ 
is an isomorphism, induced by a linear change $ \phi $ of coordinates of $R$ which maps the
$l_i$ to the the last $d$ variables of $R$; and $G = F^{-1}$.

!!! warning
    The algorithm may not terminate over a small finite field. If it terminates, the result is correct.

"""
function noether_normalization(A::MPolyQuo)
 if !(coefficient_ring(A) isa AbstractAlgebra.Field)
     throw(ArgumentError("The coefficient ring must be a field."))
 end
 if A.R isa MPolyRing_dec && !(is_standard_graded(A.R))
   throw(ArgumentError("If the base ring is decorated, it must be standard graded."))
 end
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
 h = hom(R, R, i1)
 V = map(x->h(x), gens(I))
 B, _ = quo(R, ideal(R, V))
 h1 = hom(A, B, map(B, i1))
 h2 = hom(B, A, map(A, mi_arr))
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

  if R isa MPolyRing_dec
    throw(ArgumentError("Not implemented for decorated rings."))
  end
  
  if nvars(R) != 2
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

  SR = singular_poly_ring(R)
  l = Singular.LibIntegralbasis.integralBasis(SR(f), i, "isIrred")
  A, p = quo(R, ideal(R, [f]))
  ###return (R(l[2]), R.(gens(l[1])))
  return (p(R(l[2])), [p(R(x)) for x = gens(l[1])])
end

