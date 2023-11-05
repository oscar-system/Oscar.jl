##############################################################################
#
# Data associated to affine algebras
#
##############################################################################

@doc raw"""
    dim(A::MPolyQuoRing)

Return the Krull dimension of `A`.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

julia> A, _ = quo(R, ideal(R, [y-x^2, z-x^3]));

julia> dim(A)
1
```
"""
function dim(A::MPolyQuoRing) 
  I = A.I
  return dim(I)
end


@doc raw"""
    vector_space_dimension(A::MPolyQuoRing)

If, say, `A = R/I`, where `R` is a multivariate polynomial ring over a field
`K`, and `I` is a zero-dimensional ideal of `R`, return the dimension of `A` 
as a `K`-vector space.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

julia> A, _ = quo(R, ideal(R, [x^3+y^3+z^3-1, x^2+y^2+z^2-1, x+y+z-1]));

julia> vector_space_dimension(A)
6

julia> I = modulus(A)
ideal(x^3 + y^3 + z^3 - 1, x^2 + y^2 + z^2 - 1, x + y + z - 1)

julia> groebner_basis(I, ordering = lex(base_ring(I)))
Gröbner basis with elements
1 -> z^3 - z^2
2 -> y^2 + y*z - y + z^2 - z
3 -> x + y + z - 1
with respect to the ordering
lex([x, y, z])
```
"""
function vector_space_dimension(A::MPolyQuoRing)
  if !isa(coefficient_ring(A), AbstractAlgebra.Field)
    error("vector_space_dimension requires a coefficient ring that is a field")
  end
  I = A.I
  G = standard_basis(I)
  @req dim(I) == 0 "The ideal must be zero-dimensional"
  return Singular.vdim(singular_generators(G, G.ord))
end

@doc raw"""
    monomial_basis(A::MPolyQuoRing)

If, say, `A = R/I`, where `R` is a multivariate polynomial ring over a field
`K`, and `I` is a zero-dimensional ideal of `R`, return a vector of monomials of `R` 
such that the residue classes of these monomials form a basis of `A` as a `K`-vector
space.

# Examples
```jldoctest
julia> R, (x, y) = graded_polynomial_ring(QQ, ["x", "y"]);

julia> I = ideal(R, [x^2, y^3])
ideal(x^2, y^3)

julia> A, _ = quo(R, I)
(Quotient of multivariate polynomial ring by ideal(x^2, y^3), Map: graded multivariate polynomial ring -> quotient of multivariate polynomial ring)

julia> L = monomial_basis(A)
6-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x*y^2
 y^2
 x*y
 y
 x
 1
```
"""
function monomial_basis(A::MPolyQuoRing)
  @req coefficient_ring(A) isa AbstractAlgebra.Field "The coefficient ring must be a field"
  I = A.I
  G = standard_basis(I)
  @req dim(I) == 0 "The ideal must be zero-dimensional"
  si = Singular.kbase(singular_generators(G, G.ord))
  return gens(MPolyIdeal(base_ring(I), si))
end


@doc raw"""
    monomial_basis(A::MPolyQuoRing, g::GrpAbFinGenElem)

Given an affine algebra `A` over a field which is graded by a free
group of type `GrpAbFinGen`, and given an element `g` of that group,
return a vector of monomials of `R` such that the residue classes of 
these monomials form a `K`-basis of the graded part of `A` of degree `g`.

    monomial_basis(A::MPolyQuoRing, W::Vector{<:IntegerUnion})

Given a $\mathbb  Z^m$-graded affine algebra `A` over a field and
a vector `W` of $m$ integers, convert `W` into an element `g` of the grading
group of `A` and proceed as above.

    monomial_basis(A::MPolyQuoRing, d::IntegerUnion)

Given a $\mathbb  Z$-graded  affine algebra `A` over a field and
an integer `d`, convert `d` into an element `g` of the grading
group of `A` and proceed as above.

!!! note
    If the component of the given degree is not finite dimensional, an error message will be thrown.

# Examples
```jldoctest
julia> R, (x, y) = graded_polynomial_ring(QQ, ["x", "y"]);

julia> I = ideal(R, [x^2])
ideal(x^2)

julia> A, _ = quo(R, I)
(Quotient of multivariate polynomial ring by ideal(x^2), Map: graded multivariate polynomial ring -> quotient of multivariate polynomial ring)

julia> L = monomial_basis(A, 3)
2-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 y^3
 x*y^2
```
"""
function monomial_basis(A::MPolyQuoRing, g::GrpAbFinGenElem)
  @req coefficient_ring(A) isa AbstractAlgebra.Field "The coefficient ring must be a field"
  R = base_ring(A)
  @req is_graded(R) "The ring must be graded"
  L = monomial_basis(R, g)
  LI = leading_ideal(A.I)
  ### TODO: Decide whether we should check whether a GB with respect
    ### to whatever <ordering is already available
  L = [x for x=L if !(x in LI)]
    return L
end

function monomial_basis(A::MPolyQuoRing, g::Vector{<:IntegerUnion})
  @assert is_zm_graded(A)
  return monomial_basis(A, grading_group(A)(g))
end

function monomial_basis(A::MPolyQuoRing, g::IntegerUnion)
  @assert is_z_graded(A)
  return monomial_basis(A, grading_group(A)([g]))
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


# TODO: The function below now also works for rings which are not standard graded 
# by virtue of Abbott's implementation. Clean up the docstring accordingly. 
@doc raw"""
    hilbert_series(A::MPolyQuoRing; backend::Symbol=:Singular, algorithm::Symbol=:BayerStillmanA)

Given a $\mathbb Z$-graded affine algebra $A = R/I$ over a field $K$, where the grading 
is inherited from a $\mathbb Z$-grading on the polynomial ring $R$ defined by assigning 
positive integer weights to the variables, return a pair $(p,q)$, say, of univariate 
polynomials $p, q\in\mathbb Z[t]$ such that $p/q$ represents the Hilbert series of $A$ as 
a rational function with denominator 

$q = (1-t^{w_1})\cdots (1-t^{w_n}),$

where $n$ is the number of variables of $R$, and $w_1, \dots, w_n$ are the assigned weights.

See also `hilbert_series_reduced`.

!!! note 
    The advanced user can select different backends for the computation (`:Singular` and 
    `:Abbott` for the moment), as well as different algorithms. The latter might be 
    ignored for certain backends. 

# Examples
```jldoctest
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, ["w", "x", "y", "z"]);

julia> A, _ = quo(R, ideal(R, [w*y-x^2, w*z-x*y, x*z-y^2]));

julia> hilbert_series(A)
(2*t^3 - 3*t^2 + 1, (-t + 1)^4)

julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"], [1, 2, 3]);

julia> A, _ = quo(R, ideal(R, [x*y*z]));

julia> hilbert_series(A)
(-t^6 + 1, (-t^2 + 1)^1*(-t + 1)^1*(-t^3 + 1)^1)
```
"""
function hilbert_series(A::MPolyQuoRing; #=backend::Symbol=:Singular, algorithm::Symbol=:BayerStillmanA,=# parent::Union{Nothing,Ring}=nothing)
  R = base_ring(A.I)
  @req is_z_graded(R) "ring must be graded by the integers"
  parent, t = (parent === nothing) ? polynomial_ring(ZZ, "t") : (parent, first(gens(parent)));
  W = R.d
  W = [Int(W[i][1]) for i = 1:ngens(R)]
  @req minimum(W) > 0 "The weights must be positive"
  # if iszero(A.I)
  #   den = prod([1-t^Int(w[1]) for w in R.d])
  #   return (one(parent(t)), den)
  # end

  (numer, denom), _ = multi_hilbert_series(A; parent=parent)
  return numer,denom
end

# TODO: The method below is missing. It should be made better and put to the correct place (AA).
ngens(S::AbstractAlgebra.Generic.LaurentMPolyWrapRing) = length(gens(S))
ngens(S::AbstractAlgebra.Generic.LaurentPolyWrapRing) = 1
ngens(P::PolyRing) = 1


@doc raw"""
    hilbert_series_reduced(A::MPolyQuoRing)

Given a $\mathbb Z$-graded affine algebra $A = R/I$ over a field $K$, where the grading 
is inherited from a $\mathbb Z$-grading on the polynomial ring $R$ defined by assigning 
positive integer weights to the variables, return a pair $(p,q)$, say, of univariate 
polynomials $p, q\in\mathbb Z[t]$ such that $p/q$ represents the Hilbert series of 
$A$ as a rational function written in lowest terms. 

See also `hilbert_series`.

# Examples
```jldoctest
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, ["w", "x", "y", "z"]);

julia> A, _ = quo(R, ideal(R, [w*y-x^2, w*z-x*y, x*z-y^2]));

julia> hilbert_series_reduced(A)
(2*t + 1, t^2 - 2*t + 1)

julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"], [1, 2, 3]);

julia> A, _ = quo(R, ideal(R, [x*y*z]));

julia> hilbert_series(A)
(-t^6 + 1, (-t^2 + 1)^1*(-t + 1)^1*(-t^3 + 1)^1)

julia> hilbert_series_reduced(A)
(t^2 - t + 1, t^2 - 2*t + 1)
```
"""
function hilbert_series_reduced(A::MPolyQuoRing)
   if iszero(A.I)
      return hilbert_series(A)
   end
   H = HilbertData(A.I)
   return hilbert_series_reduced(H)
end

@doc raw"""
    hilbert_series_expanded(A::MPolyQuoRing, d::Int)

Given a $\mathbb Z$-graded affine algebra $A = R/I$ over a field $K$, where the grading 
is inherited from a $\mathbb Z$-grading on the polynomial ring $R$ defined by assigning 
positive integer weights to the variables, return the Hilbert series of $A$ to precision $d$. 

# Examples
```jldoctest
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, ["w", "x", "y", "z"]);

julia> A, _ = quo(R, ideal(R, [w*y-x^2, w*z-x*y, x*z-y^2]));

julia> hilbert_series_expanded(A, 7)
1 + 4*t + 7*t^2 + 10*t^3 + 13*t^4 + 16*t^5 + 19*t^6 + 22*t^7 + O(t^8)

julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"], [1, 2, 3]);

julia> A, _ = quo(R, ideal(R, [x*y*z]));

julia> hilbert_series_expanded(A, 5)
1 + t + 2*t^2 + 3*t^3 + 4*t^4 + 5*t^5 + O(t^6)
```
"""
function hilbert_series_expanded(A::MPolyQuoRing, d::Int)
  if iszero(modulus(A))
    R = base_ring(A)
    @req is_z_graded(R) "The base ring must be ZZ-graded"
    W = R.d
    W = [Int(W[i][1]) for i = 1:ngens(R)]
    @req minimum(W) > 0 "The weights must be positive"
    num, denom = hilbert_series(A)
    T, t = power_series_ring(QQ, d+1, "t")
    return _rational_function_to_power_series(T, num, evaluate(denom))
  end
  H = HilbertData(A.I)  
  return hilbert_series_expanded(H, d)
end

@doc raw"""
    hilbert_function(A::MPolyQuoRing, d::Int)

Given a $\mathbb Z$-graded affine algebra $A = R/I$ over a field $K$, where the grading 
is inherited from a $\mathbb Z$-grading on the polynomial ring $R$ defined by assigning 
positive integer weights to the variables, return the value $H(A, d),$ where 

$H(A, \underline{\phantom{d}}): \N \to \N, \; d  \mapsto \dim_K A_d,$ 

is the Hilbert function of $A$.

# Examples
```jldoctest
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, ["w", "x", "y", "z"]);

julia> A, _ = quo(R, ideal(R, [w*y-x^2, w*z-x*y, x*z-y^2]));

julia> hilbert_function(A,7)
22

julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"], [1, 2, 3]);

julia> A, _ = quo(R, ideal(R, [x*y*z]));

julia> hilbert_function(A, 5)
5
```
"""
function hilbert_function(A::MPolyQuoRing, d::Int)
   if iszero(A.I)
       d < 0 && QQ(0)
       HS = hilbert_series_expanded(A, d)
       return coeff(HS, d)
     end
   H = HilbertData(A.I)
   return hilbert_function(H, d)
end
   
@doc raw"""
     hilbert_polynomial(A::MPolyQuoRing)

Given a $\mathbb Z$-graded affine algebra $A = R/I$ over a field $K$, where the grading 
is inherited from the standard $\mathbb Z$-grading on the polynomial ring $R$,
return the Hilbert polynomial of $A$.

# Examples
```jldoctest
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, ["w", "x", "y", "z"]);

julia> A, _ = quo(R, ideal(R, [w*y-x^2, w*z-x*y, x*z-y^2]));

julia> hilbert_polynomial(A)
3*t + 1
```
"""
function hilbert_polynomial(A::MPolyQuoRing)
   if iszero(A.I)
       R = base_ring(A.I)
       @req is_standard_graded(R) "The base ring must be standard ZZ-graded"
       n = ngens(A)
       Qt, t = QQ["t"]
       b = one(Qt)
       for i in QQ(1):QQ(n-1)
           b = b * (t+i)
       end
       b = b/QQ(factorial(n-1))
       return b
   end
   H = HilbertData(A.I)
   return hilbert_polynomial(H)::QQPolyRingElem
end

@doc raw"""
    degree(A::MPolyQuoRing)

Given a $\mathbb Z$-graded affine algebra $A = R/I$ over a field $K$, where the grading 
is inherited from the standard $\mathbb Z$-grading on the polynomial ring $R$,
return the degree of $A$.

# Examples
```jldoctest
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, ["w", "x", "y", "z"]);

julia> A, _ = quo(R, ideal(R, [w*y-x^2, w*z-x*y, x*z-y^2]));

julia> degree(A)
3
```
"""
function degree(A::MPolyQuoRing)
   if iszero(A.I)
       R = base_ring(A.I)
       @req is_standard_graded(R) "The base ring must be standard ZZ-graded"
       return ZZ(1)
     end
   H = HilbertData(A.I)
   return degree(H)
end

###############################################################################
### General Hilbert series stuff using Singular for computing ideal quotients
###############################################################################

### TODO: Originally meant to be used in multi_hilbert_series; might be useful elsewhere
function transform_to_positive_orthant(rs::Matrix{Int})   
    C = positive_hull(rs)
    @assert is_fulldimensional(C) "Cone spanned by generator degrees needs to be full-dimensional"
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

function _numerator_monomial_multi_hilbert_series(I::MPolyIdeal, S, m; algorithm::Symbol=:BayerStillmanA)
  x = gens(base_ring(I))
  W = [degree(Vector{Int}, x[i])[j] for j in 1:m, i in 1:length(x)]
  return _hilbert_numerator_from_leading_exponents([AbstractAlgebra.leading_exponent_vector(f) for f in gens(I)], W, S, :BayerStillmanA)
end


@doc raw"""
    multi_hilbert_series(A::MPolyQuoRing; algorithm::Symbol=:BayerStillmanA, parent::Union{Nothing,Ring}=nothing)

Return the Hilbert series of the graded affine algebra `A`.

!!! note 
    The advanced user can select an `algorithm` for the computation; 
    see the code for details.

# Examples
```jldoctest
julia> W = [1 1 1; 0 0 -1];

julia> R, x = graded_polynomial_ring(QQ, ["x[1]", "x[2]", "x[3]"], W)
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], x[3]])

julia> I = ideal(R, [x[1]^3*x[2], x[2]*x[3]^2, x[2]^2*x[3], x[3]^4]);

julia> A, _ = quo(R, I);

julia> H = multi_hilbert_series(A);

julia> H[1][1]
-t[1]^7*t[2]^-2 + t[1]^6*t[2]^-1 + t[1]^6*t[2]^-2 + t[1]^5*t[2]^-4 - t[1]^4 + t[1]^4*t[2]^-2 - t[1]^4*t[2]^-4 - t[1]^3*t[2]^-1 - t[1]^3*t[2]^-2 + 1

julia> H[1][2]
(-t[1] + 1)^2*(-t[1]*t[2]^-1 + 1)^1

julia> H[2][1]
GrpAb: Z^2

julia> H[2][2]
Identity map
  of GrpAb: Z^2

julia> G = abelian_group(ZZMatrix([1 -1]));

julia> g = gen(G, 1)
Element of G with components [0 1]

julia> W = [g, g, g, g];

julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, ["w", "x", "y", "z"], W);

julia> A, _ = quo(R, ideal(R, [w*y-x^2, w*z-x*y, x*z-y^2]));

julia> (num, den), (H, iso) = multi_hilbert_series(A);

julia> num
2*t^3 - 3*t^2 + 1

julia> den
(-t + 1)^4

julia> H
GrpAb: Z

julia> iso
Map
  from GrpAb: Z
  to (General) abelian group with relation matrix
  [1 -1]
  with structure of GrpAb: Z
```
"""
function multi_hilbert_series(
    A::MPolyQuoRing; 
    algorithm::Symbol=:BayerStillmanA, 
    backend::Symbol=:Abbott,
    parent::Union{Nothing, Ring}=nothing
  )
  R = base_ring(A)
  I = modulus(A)
  @req coefficient_ring(R) isa AbstractAlgebra.Field "The coefficient ring must be a field"
  @req is_positively_graded(A) "the ring must be positively graded"

  # Wrap the case where G is abstractly isomorphic to ℤᵐ, but not realized as a 
  # free Abelian group. 
  #
  # We use the Smith normal form to get there, recreate the graded ring with the 
  # free grading group, do the computation there and return the isomorphism for 
  # the grading. 
  G = grading_group(R)
  if !is_zm_graded(R)
    H, iso = snf(G)
    V = [preimage(iso, x) for x in gens(G)]
    isoinv = hom(G, H, V)
    W = [isoinv(R.d[i]) for i = 1:length(R.d)]
    S, _ = graded_polynomial_ring(coefficient_ring(R), symbols(R), W)
    map_into_S = hom(R, S, gens(S))
    J = map_into_S(I)
    AA, _ = quo(S, J)
    (numer, denom), _ = multi_hilbert_series(AA; algorithm, backend, parent)
    return (numer, denom), (H, iso)
  end

  # Now we may assume that the grading group is free Abelian.
  m = ngens(G)  
  n = ngens(R)
  HSRing = _hilbert_series_ring(parent, m)

  # Get the weights as Int values: W[k] contain the weight(s) of x[k]
  W = [[ Int(R.d[i][j])  for j in 1:m]  for i in 1:n]
  fac_denom = _hilbert_series_denominator(HSRing, W)

  # Old method below without factorization; left for debugging
  # q = one(parent)
  # for i = 1:n
  #    e = [Int(MI[i, :][j]) for j = 1:m]
  #    B = MPolyBuildCtx(parent)
  #    push_term!(B, 1, e)
  #    q = q*(1-finish(B))
  # end
  # @assert evaluate(fac_denom) == q

  # Shortcut for the trivial case
  iszero(I) && return (one(HSRing), fac_denom), (G, identity_map(G))

  # In general refer to internal methods for monomial ideals
  # TODO: Shouldn't the ordering be adapted to the grading in some sense?
  numer = one(HSRing)
  if backend == :Zach
    LI = leading_ideal(I; ordering=degrevlex(gens(R)))  # ??? better not to specify the grading ???
    numer = _numerator_monomial_multi_hilbert_series(LI, HSRing, m; algorithm=algorithm)
  elseif backend == :Abbott
    # TODO: Pass on the `algorithm` keyword argument also here.
    numer = HSNum_abbott(A, HSRing)
  else
    error("backend ($(backend)) not found")
  end
  return (numer, fac_denom), (G, identity_map(G))
end


### TODO: original version of multi_hilbert_series based on moving things to the positive orthant

#function multi_hilbert_series(A::MPolyQuoRing)
#   R = base_ring(A)
#   I = A.I
#   @req coefficient_ring(R) isa AbstractAlgebra.Field "The coefficient ring must be a field"
#   @req is_positively_graded(R) "The base ring must be positively graded"
#   if !(is_zm_graded(R))
#      G = grading_group(R)
#      H, iso = snf(G)
#      V = [preimage(iso, x) for x in gens(G)]
#      isoinv = hom(G, H, V)
#      W = R.d
#      W = [isoinv(W[i]) for i = 1:length(W)]
#      S, _ = graded_polynomial_ring(coefficient_ring(R), symbols(R), W)
#      change = hom(R, S, gens(S))
#      I = change(A.I)
#      R = S
#   end
#   m = ngens(grading_group(R))  
#   n = ngens(R)
#   W = R.d
#   MI = Matrix{Int}(undef, n, m)
#   for i=1:n
#       for j=1:m
#           MI[i, j] = Int(W[i][j])
#       end
#   end
#   minMI = minimum(MI)
#   if minMI<0
#      MI, T = transform_to_positive_orthant(MI)     
#   else
#      T = identity_matrix(ZZ, m)
#   end  
#   if m == 1
#      VAR = ["t"]
#   else
#      VAR = [_make_variable("t", i) for i = 1:m]
#   end
#   S, _ = polynomial_ring(ZZ, VAR) 
#   q = one(S)
#   for i = 1:n
#      e = [Int(MI[i, :][j]) for j = 1:m]
#      B = MPolyBuildCtx(S)
#      push_term!(B, 1, e)
#      q = q*(1-finish(B))
#   end
#   if iszero(I)
#      p = one(S)
#   else
#      LI = leading_ideal(I, ordering=degrevlex(gens(R)))
#     if minMI<0
#         RNEW, _ = graded_polynomial_ring(coefficient_ring(R), symbols(R), Matrix(transpose(MI)))
#         LI = ideal(RNEW, [RNEW(LI[i]) for i = 1:ngens(LI)])
#      end
#      p = _numerator_monomial_multi_hilbert_series(LI, S)
#   end
#   return  (p, q), T
#end

@doc raw"""
    multi_hilbert_series_reduced(A::MPolyQuoRing; algorithm::Symbol=:BayerStillmanA)

Return the reduced Hilbert series of the positively graded affine algebra `A`.

!!! note 
    The advanced user can select a `algorithm` for the computation; 
    see the code for details.

# Examples
```jldoctest
julia> W = [1 1 1; 0 0 -1];

julia> R, x = graded_polynomial_ring(QQ, ["x[1]", "x[2]", "x[3]"], W)
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], x[3]])

julia> I = ideal(R, [x[1]^3*x[2], x[2]*x[3]^2, x[2]^2*x[3], x[3]^4]);

julia> A, _ = quo(R, I);

julia> H = multi_hilbert_series_reduced(A);


julia> H[1][1]
-t[1]^5*t[2]^-1 + t[1]^3 + t[1]^3*t[2]^-3 + t[1]^2 + t[1]^2*t[2]^-1 + t[1]^2*t[2]^-2 + t[1] + t[1]*t[2]^-1 + 1

julia> H[1][2]
-t[1] + 1

julia> H[2][1]
GrpAb: Z^2

julia> H[2][2]
Identity map
  of GrpAb: Z^2

julia> G = abelian_group(ZZMatrix([1 -1]));

julia> g = gen(G, 1)
Element of G with components [0 1]

julia> W = [g, g, g, g];

julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, ["w", "x", "y", "z"], W);

julia> A, _ = quo(R, ideal(R, [w*y-x^2, w*z-x*y, x*z-y^2]));

julia> H = multi_hilbert_series_reduced(A);

julia> H[1][1]
2*t + 1

julia> H[1][2]
t^2 - 2*t + 1

julia> H[2][1]
GrpAb: Z

julia> H[2][2]
Map
  from GrpAb: Z
  to (General) abelian group with relation matrix
  [1 -1]
  with structure of GrpAb: Z
```
"""
function multi_hilbert_series_reduced(A::MPolyQuoRing; algorithm::Symbol=:BayerStillmanA)
   @req is_positively_graded(A) "ring must be positively graded"
   (p, q), (H, iso) = multi_hilbert_series(A, algorithm=algorithm)
   f = p//evaluate(q)
   p = numerator(f)
   q = denominator(f)
   sig = coeff(q.mpoly, -1 .* q.mindegs)
   return (sig*p, sig*q), (H, iso)
end

function _monomial_ideal_membership(m::MPolyRingElem, I::MPolyIdeal)
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

@doc raw"""
    multi_hilbert_function(A::MPolyQuoRing, g::GrpAbFinGenElem)

Given a positively graded affine algebra $A$ over a field $K$ with grading group $G$,
say, and given an element $g$ of $G$, return the value $H(A, g)$ of the Hilbert function

$H(A, \underline{\phantom{d}}): G \to \N, \; g\mapsto \dim_K(A_g).$

    multi_hilbert_function(A::MPolyQuoRing, g::Vector{<:IntegerUnion})

Given a positively $\mathbb  Z^m$-graded affine algebra $A$ over a field $K$,
and given a vector $g$ of $m$ integers, convert $g$ into an element 
of the grading group of $A$, and return the value $H(A, g)$
as above.

    multi_hilbert_function(A::MPolyQuoRing, g::IntegerUnion)

Given a positively $\mathbb  Z$-graded affine algebra $A$ over a field $K$,
and given an integer $g$, convert $g$ into an element of the grading group 
of $A$, and return the value $H(A, g)$ as above.

# Examples
```jldoctest
julia> W = [1 1 1; 0 0 -1];

julia> R, x = graded_polynomial_ring(QQ, ["x[1]", "x[2]", "x[3]"], W)
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], x[3]])

julia> I = ideal(R, [x[1]^3*x[2], x[2]*x[3]^2, x[2]^2*x[3], x[3]^4]);

julia> A, _ = quo(R, I);

julia> multi_hilbert_function(A::MPolyQuoRing, [1, 0])
2
```

```jldoctest
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, ["w", "x", "y", "z"], [-1, -1, -1, -1]);

julia> A, _ = quo(R, ideal(R, [w*y-x^2, w*z-x*y, x*z-y^2]));

julia> multi_hilbert_function(A, -7)
22
```

```jldoctest
julia> G = abelian_group(ZZMatrix([1 -1]));

julia> g = gen(G, 1);

julia> W = [g, g, g, g];

julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, ["w", "x", "y", "z"], W);

julia> A, _ = quo(R, ideal(R, [w*y-x^2, w*z-x*y, x*z-y^2]));

julia> multi_hilbert_function(A, 7*g)
22
```
"""
function multi_hilbert_function(A::MPolyQuoRing, g::GrpAbFinGenElem)
    R = base_ring(A)
    @req coefficient_ring(R) isa AbstractAlgebra.Field "The coefficient ring must be a field"
    
    L = monomial_basis(R, g)
    
    if size(L) == 0
       return 0
    end

    LI = leading_ideal(A.I, ordering=degrevlex(gens(R)))
    ### TODO: Decide whether we should check whether a GB with respect
    ### to another degree-compatible ordering is already available
    
    cc = 0
    for i in 1:length(L)
        if !(L[i] in LI)
            cc = cc+1
        end
    end
    return cc
end

function multi_hilbert_function(A::MPolyQuoRing, g::Vector{<:IntegerUnion})
  @assert is_zm_graded(A)
  return multi_hilbert_function(A, grading_group(A)(g))
end

function multi_hilbert_function(A::MPolyQuoRing, g::IntegerUnion)
  @assert is_z_graded(A)
  return multi_hilbert_function(A, grading_group(A)([g]))
end

##############################################################################
#
# Properties of affine algebras
#
##############################################################################

@doc raw"""
    is_reduced(A::MPolyQuoRing)

Given an affine algebra `A`, return `true` if `A` is reduced, `false` otherwise.

!!! warning
    The function computes the radical of the modulus of `A`. This may take some time.

# Examples
```jldoctest
julia> R, (x,) = polynomial_ring(QQ, ["x"]);

julia> A, _ = quo(R, ideal(R, [x^4]));

julia> is_reduced(A)
false
```
"""
function is_reduced(A::MPolyQuoRing) 
  I = A.I
  return I == radical(I)
end

@doc raw"""
    is_normal(A::MPolyQuoRing)

Given an affine algebra `A` over a perfect field,
return `true` if `A` is normal, `false` otherwise.

!!! note 
    This function performs the first step of the normalization algorithm of Greuel, Laplagne, and Seelisch [GLS10](@cite) and may, thus, be more efficient than computing the full normalization of `A`.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

julia> A, _ = quo(R, ideal(R, [z^2-x*y]));

julia> is_normal(A)
true
```
"""
function is_normal(A::MPolyQuoRing)
  @req coefficient_ring(A) isa AbstractAlgebra.Field "The coefficient ring must be a field"
  @req !(base_ring(A) isa MPolyDecRing) "Not implemented for quotients of decorated rings"

  I = A.I
  # TODO remove old1 & old2 once new Singular jll is out
  old1 = Singular.libSingular.set_option("OPT_REDSB", false)
  old2 = Singular.libSingular.set_option("OPT_RETURN_SB", false)
  f = Singular.LibNormal.isNormal(singular_generators(I))::Int
  Singular.libSingular.set_option("OPT_REDSB", old1)
  Singular.libSingular.set_option("OPT_RETURN_SB", old2)
  return Bool(f)
end

@doc raw"""
     is_cohen_macaulay(A::MPolyQuoRing) 

Given a $\mathbb Z$-graded affine algebra `A = R/I` over a field, say, `K`, where the grading 
is inherited from the standard $\mathbb Z$-grading on the polynomial ring `R`,
return `true` if `A` is a Cohen-Macaulay ring, `false` otherwise.

# Examples
```jldoctest
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, ["w", "x", "y", "z"]);

julia> I = ideal(R, [x*z-y^2, w*z-x*y, w*y-x^2]);

julia> A, _ = quo(R, I);

julia> is_cohen_macaulay(A)
true
```

```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> I = ideal(R, [x*z, y*z]);

julia> A, _ = quo(R, I);

julia> is_cohen_macaulay(A)
false
```
"""
function is_cohen_macaulay(A::MPolyQuoRing)
 I = A.I
 R = base_ring(I)
 @req coefficient_ring(R) isa AbstractAlgebra.Field "The coefficient ring must be a field"
 @req is_standard_graded(R) "The base ring must be standard ZZ-graded"

 sI = singular_generators(I.gens, negdegrevlex(gens(R)))
 res = Singular.LibHomolog.isCM(sI)
 if res == 1 return true end
 return false
end

##############################################################################
#
# Algebra Containment
#
##############################################################################

@doc raw"""
    subalgebra_membership(f::T, V::Vector{T}) where T <: Union{MPolyRingElem, MPolyQuoRingElem}

Given an element `f` of a multivariate polynomial ring over a field, or of a
quotient of such a ring, and given a vector `V` of further elements of that ring,
consider the subalgebra generated by the entries of `V` in the given ring. If `f`
is contained in the subalgebra, return `(true, h)`, where `h` is giving the
polynomial relation. Return, `(false, 0)`, otherwise.

# Examples
```jldoctest
julia> R, x = polynomial_ring(QQ, "x" => 1:3);

julia> f = x[1]^6*x[2]^6-x[1]^6*x[3]^6;

julia> V = [x[1]^3*x[2]^3-x[1]^3*x[3]^3, x[1]^3*x[2]^3+x[1]^3*x[3]^3]
2-element Vector{QQMPolyRingElem}:
 x[1]^3*x[2]^3 - x[1]^3*x[3]^3
 x[1]^3*x[2]^3 + x[1]^3*x[3]^3

julia> subalgebra_membership(f, V)
(true, t1*t2)
```
"""
function subalgebra_membership(f::T, v::Vector{T}) where T <: Union{MPolyRingElem, MPolyQuoRingElem}
  R = parent(f)
  @req coefficient_ring(R) isa Field "The coefficient ring must be a field"
  @req !isempty(v) "Input vector must not be empty"
  @req all(x -> parent(x) === R, v) "The polynomials must have the same parent"

  S, _ = polynomial_ring(coefficient_ring(R), length(v), "t")
  phi = hom(S, R, v)
  return has_preimage(phi, f)
end

@doc raw"""
    subalgebra_membership_homogeneous(f::T, V::Vector{T}; check::Bool = true)
      where T <: Union{MPolyDecRingElem, MPolyQuoRingElem{TT} where TT <: MPolyDecRingElem}

Given a homogeneous element `f` of a positively `Z`-graded multivariate polynomial
ring over a field or a quotient of such a ring, and given a vector `V` of homogeneous
elements in the same ring, consider the subalgebra generated by the entries of `V`
in that ring. If `f` is contained in the subalgebra, return `(true, h)`, where `h`
is giving the polynomial relation. Return, `(false, 0)`, otherwise.

If `check` is `true` (default), the homogeneity of all given polynomials is
checked.
"""
function subalgebra_membership_homogeneous(f::PolyRingElemT, v::Vector{PolyRingElemT}; check::Bool = true) where PolyRingElemT <: MPolyDecRingElem
  return _subalgebra_membership_homogeneous(f, v, ideal(parent(f), [ zero(parent(f)) ]), check = check)
end

function subalgebra_membership_homogeneous(f::PolyRingElemT, v::Vector{PolyRingElemT}; check::Bool = true) where PolyRingElemT <: MPolyQuoRingElem{T} where T <: MPolyDecRingElem
  return _subalgebra_membership_homogeneous(lift(f), [ lift(g) for g in v ], modulus(parent(f)), check = check)
end

function _subalgebra_membership_homogeneous(f::PolyRingElemT, v::Vector{PolyRingElemT}, I::MPolyIdeal{PolyRingElemT}; check::Bool = true) where PolyRingElemT <: MPolyDecRingElem
  R = parent(f)
  @req !isempty(v) "Input vector must not be empty"
  @req all(g -> parent(g) === R, v) "The polynomials must have the same parent"
  if check
    @req is_z_graded(R) "The base ring must be Z-graded"
    @req is_positively_graded(R) "The base ring must be positively graded"
    @req is_homogeneous(f) "The input must be homogeneous"
    @req all(is_homogeneous, v) "The input must be homogeneous"
  end

  # This is basically [GP09, p. 86, Solution 2], but we only compute a degree
  # truncated Gröbner basis

  T, _ = polynomial_ring(base_ring(R), ngens(R) + length(v))

  RtoT = hom(R, T, gens(T)[1:ngens(R)])

  J = RtoT(I) + ideal(T, [ RtoT(v[i]) - gen(T, ngens(R) + i) for i = 1:length(v) ])

  o = degrevlex(gens(T)[1:ngens(R)])*wdegrevlex(gens(T)[ngens(R) + 1:ngens(T)], [ Int(degree(g)[1]) for g in v ])

  # Everything is homogeneous, so a truncated Gröbner basis up to the degree
  # of f suffices to check containment of f in J
  GJ = _groebner_basis(J, Int(degree(f)[1]), ordering = o)

  ###
  # This computes the normal form of f w.r.t. the truncated Gröbner basis GJ.
  # Since we have a product ordering, we cannot use divrem, and since GJ is
  # "not really" a Gröbner basis, we cannot use normal_form...
  SR = singular_polynomial_ring(GJ)
  I = Singular.Ideal(SR, SR(RtoT(f)))
  K = ideal(T, reduce(I, singular_generators(GJ, GJ.ord)))
  @assert is_one(ngens(K.gens.S))
  nf = GJ.Ox(K.gens.S[1])
  ###

  S, _ = polynomial_ring(base_ring(R), [ "t$i" for i in 1:length(v) ])
  TtoS = hom(T, S, append!(zeros(S, ngens(R)), gens(S)))

  # f is in the subalgebra iff nf does not involve the variables
  # gen(T, 1), ..., gen(T, ngens(R)), that is, iff LM(nf) is strictly smaller
  # than gen(T, ngens(R)) w.r.t. the block ordering o.
  if isone(cmp(o, gen(T, ngens(R)), leading_monomial(nf, ordering = o)))
    return true, TtoS(nf)
  else
    return false, zero(S)
  end
end

################################################################################
#
#  Algebraic Independence
#
################################################################################

@doc raw"""
    are_algebraically_independent(V::Vector{T}) where T <: Union{MPolyRingElem, MPolyQuoRingElem}

Given a vector `V` of elements of a multivariate polynomial ring over a field `K`, say, or of a quotient of such a ring, 
return `(true, ideal(0))` if the elements of `V` are algebraically independent over `K`. Return, `false`
together with the ideal of `K`-algebra relations, otherwise.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> V = [x, y, x^2+y^3]
3-element Vector{QQMPolyRingElem}:
 x
 y
 x^2 + y^3

julia> are_algebraically_independent(V)
(false, ideal(t1^2 + t2^3 - t3))

julia> A, p = quo(R, [x*y]);

julia> are_algebraically_independent([p(x), p(y)])
(false, ideal(t1*t2))

```
"""
function are_algebraically_independent(V::Vector{T}) where T <: Union{MPolyRingElem, MPolyQuoRingElem}
  @req !isempty(V) "Input vector must not be empty"
  R = parent(V[1])
  @req coefficient_ring(R) isa Field "The coefficient ring must be a field"
  @req all(x -> parent(x) === R, V) "The elements must have the same parent"
  S, _ = polynomial_ring(coefficient_ring(R), length(V), "t"; cached = false)
  phi = hom(S, R, V)
  I = kernel(phi)
  return iszero(I), I
end

################################################################################
#
#  Minimalizing a set of subalgebra generators in graded case
#
################################################################################

@doc raw"""
    minimal_subalgebra_generators(V::Vector{T}; check::Bool = true) where T <: Union{MPolyRingElem, MPolyQuoRingElem}

Given a vector `V` of homogeneous elements of a positively graded multivariate
polynomial ring, or of a quotient of such a ring, return a minimal subset of the
elements in `V` which, in the given ring, generate
the same subalgebra as all elements in `V`.

If `check` is `true` (default), the conditions on `V` and the given ring are
checked.

# Examples
```jldoctest
julia> R, (x, y) = graded_polynomial_ring(QQ, ["x", "y"]);

julia> V = [x, y, x^2+y^2]
3-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x
 y
 x^2 + y^2

julia> minimal_subalgebra_generators(V)
2-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x
 y
```
"""
function minimal_subalgebra_generators(V::Vector{T}; check::Bool = true) where {T <: Union{MPolyDecRingElem, MPolyQuoRingElem{<: MPolyDecRingElem}}}
  @req !isempty(V) "Input vector must not be empty"
  I = ideal(parent(V[1]), [ zero(parent(V[1])) ])
  return _minimal_subalgebra_generators_with_relations(V, I, check = check)[1]
end

function minimal_subalgebra_generators_with_relations(V::Vector{<: MPolyDecRingElem}; check::Bool = true)
  @req !isempty(V) "Input vector must not be empty"
  I = ideal(parent(V[1]), [ zero(parent(V[1])) ])
  return _minimal_subalgebra_generators_with_relations(V, I, check = check)
end

function minimal_subalgebra_generators_with_relations(V::Vector{<: MPolyQuoRingElem{T}}; check::Bool = true) where T <: MPolyDecRingElem
  @req !isempty(V) "Input vector must not be empty"
  return _minimal_subalgebra_generators_with_relations([ lift(f) for f in V ], modulus(parent(V[1])), check = check)
end

function _minimal_subalgebra_generators_with_relations(V::Vector{PolyRingElemT}, I::MPolyIdeal{PolyRingElemT}; check::Bool = true, start::Int = 0) where PolyRingElemT <: MPolyDecRingElem
  @req !isempty(V) "Input vector must not be empty"
  @assert start >= 0
  R = parent(V[1])
  @req all(g -> parent(g) === R, V) "The polynomials must have the same parent"
  if check
    @req is_z_graded(R) "The base ring must be Z-graded"
    @req is_positively_graded(R) "The base ring must be positively graded"
    @req all(is_homogeneous, V) "The input must be homogeneous"
  end

  K = coefficient_ring(R)

  # Use S to keep track of relations
  S, _ = polynomial_ring(K, length(V), "t")
  rels = Vector{elem_type(S)}(undef, length(V))

  # If start > 0, we assume that the polynomials V[1:start] are known to be part
  # of a minimal generating set
  if start > 0
    W = V[start + 1:end]
    # Sort by increasing degree
    sp = sortperm(W, lt = (x, y) -> degree(x)[1] < degree(y)[1])
    W = W[sp]

    res = V[1:start]
    for i in 1:start
      rels[i] = gen(S, i)
    end
  else
    sp = sortperm(V, lt = (x, y) -> degree(x)[1] < degree(y)[1])
    W = V[sp]

    res = elem_type(R)[ W[1] ]
    rels[sp[1]] = gen(S, 1)
  end

  for i in 1:length(W)
    f = W[i]
    fl, t = _subalgebra_membership_homogeneous(f, res, I, check = false)
    if fl
      # f is in the span of the generators so far
      rels[start + sp[i]] = t(gens(S)[1:length(res)]...)
    else
      # f is a new generator
      rels[start + sp[i]] = gen(S, start + i)
      push!(res, f)
    end
  end

  # S might have too many variables
  if length(V) > length(res)
    S2, _ = polynomial_ring(K, length(res), "t")
    t = append!(gens(S2), [ zero(S2) for i in 1:ngens(S) - ngens(S2) ])
    rels = [ f(t...) for f in rels ]
  end

  return res, rels
end

##############################################################################
#
# Normalization
#
##############################################################################

function _conv_normalize_alg(algorithm::Symbol)
  if algorithm == :primeDec
    return "prim"
  elseif algorithm == :equidimDec
    return "equidim"
  else
    error("algorithm invalid")
  end
end

function _conv_normalize_data(A::MPolyQuoRing, l, br)
  return [
    begin
      newSR = l[1][i][1]::Singular.PolyRing
      newOR, _ = polynomial_ring(br, [string(x) for x in gens(newSR)])
      newA, newAmap = quo(newOR, ideal(newOR, l[1][i][2][:norid]))
      newgens = newOR.(gens(l[1][i][2][:normap]))
      _hom = hom(A, newA, newA.(newgens))
      idgens = base_ring(A).(gens(l[2][i]))
      (newA, _hom, (A(idgens[end]), ideal(A, idgens)))
    end
    for i in 1:length(l[1])]
end

@doc raw"""
    normalization(A::MPolyQuoRing; algorithm = :equidimDec)

Find the normalization of a reduced affine algebra over a perfect field $K$.
That is, given the quotient $A=R/I$ of a multivariate polynomial ring $R$ over $K$
modulo a radical ideal $I$, compute the integral closure $\overline{A}$ 
of $A$ in its total ring of fractions $Q(A)$, together with the embedding 
$f: A \to \overline{A}$. 

# Implemented Algorithms and how to Read the Output

The function relies on the algorithm 
of Greuel, Laplagne, and Seelisch which proceeds by finding a suitable decomposition 
$I=I_1\cap\dots\cap I_r$ into radical ideals $I_k$, together with
maps $A = R/I \to A_k=\overline{R/I_k}$ which give rise to the normalization map of $A$:

$A\hookrightarrow A_1\times \dots\times A_r=\overline{A}$

For each $k$, the function specifies two representations
of $A_k$: It returns an array of triples $(A_k, f_k, \mathfrak a_k)$,
where $A_k$ is represented as an affine $K$-algebra, and $f_k$ as a map of affine $K$-algebras.
The third entry $\mathfrak a_k$ is a tuple $(d_k, J_k)$, consisting of an element
$d_k\in A$ and an ideal $J_k\subset A$, such that $\frac{1}{d_k}J_k = A_k$ 
as $A$-submodules of the total ring of fractions of $A$.

By default (`algorithm = :equidimDec`), as a first step on its way to find the decomposition $I=I_1\cap\dots\cap I_r$, 
the algorithm computes an equidimensional decomposition of the radical ideal $I$.
Alternatively, if specified by `algorithm = :primeDec`, the algorithm computes $I=I_1\cap\dots\cap I_r$
as the prime decomposition of the radical ideal $I$.

See [GLS10](@cite).

!!! warning
    The function does not check whether $A$ is reduced. Use `is_reduced(A)` in case you are unsure (this may take some time).

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> A, _ = quo(R, ideal(R, [(x^2-y^3)*(x^2+y^2)*x]));

julia> L = normalization(A);

julia> size(L)
(2,)

julia> LL = normalization(A, algorithm = :primeDec);

julia> size(LL)
(3,)

julia> LL[1][1]
Quotient
  of multivariate polynomial ring in 3 variables T(1), x, y
    over rational field
  by ideal(-T(1)*y + x, -T(1)*x + y^2, T(1)^2 - y, -x^2 + y^3)

julia> LL[1][2]
Ring homomorphism
  from quotient of multivariate polynomial ring by ideal(x^5 - x^3*y^3 + x^3*y^2 - x*y^5)
  to quotient of multivariate polynomial ring by ideal(-T(1)*y + x, -T(1)*x + y^2, T(1)^2 - y, -x^2 + y^3)
defined by
  x -> x
  y -> y

julia> LL[1][3]
(y, ideal(x, y))
```
"""
function normalization(A::MPolyQuoRing; algorithm=:equidimDec)
  @req coefficient_ring(A) isa AbstractAlgebra.Field "The coefficient ring must be a field"
  @req !(base_ring(A) isa MPolyDecRing) "Not implemented for quotients of decorated rings"

  I = A.I
  br = base_ring(base_ring(A))
  l = Singular.LibNormal.normal(singular_generators(I), _conv_normalize_alg(algorithm))
  return _conv_normalize_data(A, l, br)
end

@doc raw"""
    normalization_with_delta(A::MPolyQuoRing; algorithm::Symbol = :equidimDec)

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
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> A, _ = quo(R, ideal(R, [(x^2-y^3)*(x^2+y^2)*x]));

julia> L = normalization_with_delta(A);

julia> L[2]
3-element Vector{Int64}:
 1
 1
 0

julia> L[3]
13
```
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

julia> A, _ = quo(R, ideal(R, [z^3-x*y^4]));

julia> L = normalization_with_delta(A);

julia> L[3]
-1
```
"""
function normalization_with_delta(A::MPolyQuoRing; algorithm::Symbol=:equidimDec)
  @req coefficient_ring(A) isa AbstractAlgebra.Field "The coefficient ring must be a field"
  @req !(base_ring(A) isa MPolyDecRing) "Not implemented for quotients of decorated rings"

  I = A.I
  br = base_ring(base_ring(A))
  l = Singular.LibNormal.normal(singular_generators(I), _conv_normalize_alg(algorithm), "withDelta")
  return (_conv_normalize_data(A, l, br), l[3][1]::Vector{Int}, l[3][2]::Int)
end


##############################################################################
#
# Noether Normalization
#
##############################################################################

@doc raw"""
    noether_normalization(A::MPolyQuoRing)

Given an affine algebra $A=R/I$ over a field $K$, return a triple $(V,F,G)$ such that:
- ``V`` is a vector of $d=\dim A$ elements of $A$, represented by linear forms $l_i\in R$, and such that $K[V]\hookrightarrow A$ is a Noether normalization for $A$; 
- ``F: A=R/I \to B = R/\phi(I)`` is an isomorphism, induced by a linear change $ \phi $ of coordinates of $R$ which maps the $l_i$ to the the last $d$ variables of $R$; 
- ``G = F^{-1}.``

!!! warning
    The algorithm may not terminate over a small finite field. If it terminates, the result is correct.

"""
function noether_normalization(A::MPolyQuoRing)
  @req coefficient_ring(A) isa AbstractAlgebra.Field "The coefficient ring must be a field"
  if base_ring(A) isa MPolyDecRing
    @req is_standard_graded(A) "If the base ring is decorated, it must be standard graded"
  end

 I = A.I
 R = base_ring(I)
 l = Singular.LibAlgebra.noetherNormal(singular_generators(I))
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

@doc raw"""
    integral_basis(f::MPolyRingElem, i::Int; algorithm::Symbol = :normal_local)

Given a polynomial $f$ in two variables with coefficients in a perfect field $K$, and
given an integer $i\in\{1,2\}$ specifying one of the variables, $f$ must be irreducible
and monic in the specified variable: Say, $f\in\mathbb K[x,y]$ is monic in $y$.
Then the normalization of $A = K[x,y]/\langle f \rangle$, that is, the
integral closure $\overline{A}$ of $A$ in its quotient field, is a free
module over $K[x]$ of finite rank, and any set of free generators for
$\overline{A}$ over $K[x]$ is called an *integral basis* for $\overline{A}$
over $K[x]$. The function returns a pair $(d, V)$, where $d$ is an element of $A$,
and $V$ is a vector of elements in $A$, such that the fractions $v/d, v\in V$,
form an integral basis for $\overline{A}$ over $K[x]$.

By default (`algorithm = :normal_local`), the function relies on the
local-to-global approach to normalization presented in [BDLPSS13](@cite).
Alternatively, if specified by `algorithm = :normal_global`, the global normalization
algorithm in [GLS10](@cite) is used. If $K = \mathbb Q$, it is recommended to
apply the algorithm in [BDLP19](@cite), which makes use of Puiseux expansions
and Hensel lifting (`algorithm = :hensel`).

!!! note
    The conditions on $f$ are automatically checked.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> f = (y^2-2)^2 + x^5
x^5 + y^4 - 4*y^2 + 4

julia> integral_basis(f, 2)
(x^2, MPolyQuoRingElem{QQMPolyRingElem}[x^2, x^2*y, y^2 - 2, y^3 - 2*y])
```
"""
function integral_basis(f::MPolyRingElem, i::Int; algorithm::Symbol = :normal_local)
  R = parent(f)

  if algorithm == :hensel
    options = ("hensel",)
  elseif algorithm == :normal_local
    options = ("normal", "local")
  elseif algorithm == :normal_global
    options = ("normal", "global")
  else
    throw(ArgumentError("unsupported algorithm $algorithm"))
  end

  @req !(R isa MPolyDecRing) "Not implemented for decorated rings"
  
  @req nvars(R) == 2 "The parent ring must be a polynomial ring in two variables"

  @req i == 1 || i == 2 "The index $i must be either 1 or 2, indicating the integral variable"

  @req isone(coeff(f, [i], [degree(f, i)])) "The input polynomial must be monic as a polynomial in $(gen(R,i))"

  SR = singular_poly_ring(R)

  if !(base_ring(SR) isa Singular.Rationals ||
       base_ring(SR) isa Singular.N_ZpField ||
       base_ring(SR) isa Singular.N_GField ||
       base_ring(SR) isa Singular.N_AlgExtField)
    throw(NotImplementedError(:integral_basis, f))
  end

  @req is_irreducible(f) "The input polynomial must be irreducible"

  l = Singular.LibIntegralbasis.integralBasis(SR(f), i, "isIrred", options...)
  A, p = quo(R, ideal(R, [f]))
  ###return (R(l[2]), R.(gens(l[1])))
  return (p(R(l[2])), [p(R(x)) for x in gens(l[1])])
end

