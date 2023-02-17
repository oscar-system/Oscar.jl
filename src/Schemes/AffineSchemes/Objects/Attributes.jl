export OO, coordinate_ring, base_ring, dim, codim, name

export ambient_space, ambient_coordinate_ring, ambient_coordinates, ambient_embedding

export ring_type, base_ring_type, base_ring_elem_type, poly_type, ring_type

export singular_locus, singular_locus_reduced

export reduced_scheme

########################################################################
# (1) Attributes of AbsSpec
#     coordinate ring and ambient space related methods
########################################################################

# Here is the interface for AbsSpec

@Markdown.doc """
    coordinate_ring(X::AbsSpec)

On an affine scheme ``X = Spec(R)`` this returns the ring ``R``.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> coordinate_ring(X)
Multivariate Polynomial Ring in x1, x2, x3 over Rational Field
```
We allow the shortcut `OO`
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> OO(X)
Multivariate Polynomial Ring in x1, x2, x3 over Rational Field
"""
coordinate_ring(X::AbsSpec) = OO(X)

@Markdown.doc """
    OO(X::AbsSpec)

On an affine scheme ``X = Spec(R)`` this returns the ring ``R``.
"""
function OO(X::AbsSpec{BRT, RT}) where {BRT, RT}
  OO(underlying_scheme(X))::RT
end


@Markdown.doc """
    ambient_space(X::AbsSpec)

Return the ambient affine space of ``X``. 

Use [`ambient_embedding(::AbsSpec)`](@ref) to obtain the embedding of ``X`` in
its ambient affine space.

# Examples
```jldoctest
julia> X = affine_space(QQ, [:x,:y])
Spec of Multivariate Polynomial Ring in x, y over Rational Field

julia> ambient_space(X) == X
true

julia> (x, y) = coordinates(X);

julia> Y = subscheme(X, [x])
Spec of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x)

julia> X == ambient_space(Y)
true

julia> Z = subscheme(Y, y)
Spec of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x, y)

julia> ambient_space(Z) == X
true

julia> V = hypersurface_complement(Y, y)
Spec of Localization of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x) at the multiplicative set powers of fmpq_mpoly[y]

julia> ambient_space(V) == X
true
```

We can create ``X``, ``Y`` and ``Z`` also by first constructing the corresponding
coordinate rings. The subset relations are inferred from the coordinate rings.
More precisely, for a polynomial ring ``P`` an ideal ``I ⊆ P `` and a multiplicatively closed subset
``U`` of ``P`` let ``R`` be one of ``P``, ``U^{-1}P``, ``P/I`` or ``U^{-1}(P/I)``.
In each case the ambient affine space is given by `Spec(P)`.

# Examples
```jldoctest ambient_via_spec
julia> P, (x, y) = PolynomialRing(QQ, [:x, :y])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> X = Spec(P)
Spec of Multivariate Polynomial Ring in x, y over Rational Field

julia> I = ideal(P, x)
ideal(x)

julia> RmodI, quotient_map = quo(P, I);

julia> Y = Spec(RmodI)
Spec of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x)

julia> ambient_space(Y) == X
true

julia> J = ideal(RmodI, y);

julia> RmodJ, quotient_map2 = quo(RmodI, J);

julia> Z = Spec(RmodJ)
Spec of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x, y)

julia> ambient_space(Z) == X
true

julia> U = powers_of_element(y)
powers of fmpq_mpoly[y]

julia> URmodI, _ = localization(RmodI, U);

julia> V = Spec(URmodI)
Spec of Localization of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x) at the multiplicative set powers of fmpq_mpoly[y]

julia> ambient_space(V) == X
true
```

Note: compare with `==`, as the same affine space could be represented
internally by different objects for technical reasons.
# Examples
```jldoctest ambient_via_spec
julia> AX = ambient_space(X);

julia> AY = ambient_space(Y);

julia> AX == AY
true

julia> AX === AY
false
```
"""
function ambient_space(X::AbsSpec)
  error("$X does not have an ambient affine space")
end

function ambient_space(X::AbsSpec{BRT, RT}) where {BRT, RT<:MPolyRing}
  return X
end

@attr function ambient_space(X::AbsSpec{BRT,RT}) where {BRT, RT <: Union{MPolyQuo,MPolyLocalizedRing,MPolyQuoLocalizedRing}}
  return Spec(ambient_coordinate_ring(X))
end

@Markdown.doc """
    ambient_embedding(X::AbsSpec)

Return the embedding of ``X`` in its ambient affine space.

# Examples
```jldoctest
julia> X = affine_space(QQ, [:x,:y])
Spec of Multivariate Polynomial Ring in x, y over Rational Field

julia> (x, y) = coordinates(X);

julia> Y = subscheme(X, [x]);

julia> inc = ambient_embedding(Y);

julia> inc == inclusion_morphism(Y, X)
true
```
"""
function ambient_embedding(X::AbsSpec)
  return inclusion_morphism(X, ambient_space(X), check=false)
end

@Markdown.doc """
    ambient_coordinate_ring(X::AbsSpec)

Return the coordinate ring of the ambient affine space of ``X``.

See also [`ambient_space(::AbsSpec)`](@ref).

# Examples
```jldoctest
julia> X = affine_space(QQ, [:x,:y])
Spec of Multivariate Polynomial Ring in x, y over Rational Field

julia> (x,y) = coordinates(X);

julia> Y = subscheme(X, [x])
Spec of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x)

julia> ambient_coordinate_ring(Y)
Multivariate Polynomial Ring in x, y over Rational Field
```
"""
function ambient_coordinate_ring(X::AbsSpec)
  return ambient_coordinate_ring(underlying_scheme(X))::MPolyRing
end

@Markdown.doc """
    ambient_coordinates(X::AbsSpec)

Return the coordinate functions of the ambient affine space of ``X``.

See also [`ambient_space(::AbsSpec)`](@ref).

# Examples
```jldoctest
julia> X = affine_space(QQ, [:x,:y])
Spec of Multivariate Polynomial Ring in x, y over Rational Field

julia> (x,y) = coordinates(X);

julia> Y = subscheme(X, [x])
Spec of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x)

julia> coordinates(X) == ambient_coordinates(Y)
true

julia> [x,y] == ambient_coordinates(Y)
true

```
"""
ambient_coordinates(X::AbsSpec) = gens(ambient_coordinate_ring(X))

@Markdown.doc """
    coordinates(X::AbsSpec)

Return the coordinate functions of ``X`` as elements of its coordinate ring.

For ``X`` a subscheme of an ambient affine space, the coordinate functions are induced
by the ambient affine space.

# Examples
```jldoctest
julia> X = affine_space(QQ, [:x,:y])
Spec of Multivariate Polynomial Ring in x, y over Rational Field

julia> (x, y) = coordinates(X)
2-element Vector{fmpq_mpoly}:
 x
 y

julia> Y = subscheme(X, [x])
Spec of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x)

julia> (xY, yY) = coordinates(Y)
2-element Vector{MPolyQuoElem{fmpq_mpoly}}:
 x
 y

julia> parent(xY) == coordinate_ring(Y)
true
```
"""
coordinates(X::AbsSpec) = gens(OO(X))

@Markdown.doc """
    base_ring(X::AbsSpec)

On an affine scheme ``X/𝕜`` over ``𝕜`` this returns the ring ``𝕜``.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> base_ring(X)
Rational Field
```
"""
function base_ring(X::AbsSpec{BRT, RT}) where {BRT, RT}
  return base_ring(underlying_scheme(X))::BRT
end

##############################################################################
# (2) Attributes of AbsSpec
#     dimension, codimension, name
##############################################################################

@Markdown.doc """
    dim(X::AbsSpec)

Return the dimension the affine scheme ``X = Spec(R)``.

By definition, this is the Krull dimension of ``R``.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> dim(X)
3

julia> Y = affine_space(ZZ, 2)
Spec of Multivariate Polynomial Ring in x1, x2 over Integer Ring

julia> dim(Y) # one dimension comes from ZZ and two from x1 and x2
3
```
"""
@attr function dim(X::AbsSpec{<:Ring, <:MPolyQuoLocalizedRing})
  return dim(saturated_ideal(modulus(OO(X))))
end

@attr function dim(X::AbsSpec{<:Ring, <:MPolyLocalizedRing})
  # the following line is supposed to refer the problem to the
  # algebra side of the problem
  return dim(ideal(ambient_coordinate_ring(X), [zero(ambient_coordinate_ring(X))]))
end

@attr function dim(X::AbsSpec{<:Ring, <:MPolyRing})
  return dim(ideal(ambient_coordinate_ring(X), [zero(ambient_coordinate_ring(X))]))
end

@attr function dim(X::AbsSpec{<:Ring, <:MPolyQuo})
  return dim(modulus(OO(X)))
end


@Markdown.doc """
    codim(X::AbsSpec)

Return the codimension of ``X`` in its ambient affine space.

Throws and error if ``X`` does not have an ambient affine space.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> codim(X)
0

julia> R = OO(X)
Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> (x1,x2,x3) = gens(R)
3-element Vector{fmpq_mpoly}:
 x1
 x2
 x3

julia> Y = subscheme(X, x1)
Spec of Quotient of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field by ideal(x1)

julia> codim(Y)
1
```
"""
@attr function codim(X::AbsSpec)
  return dim(ideal(ambient_coordinate_ring(X), [zero(ambient_coordinate_ring(X))])) - dim(X)
end


@doc Markdown.doc"""
    name(X::AbsSpec)

Return the current name of an affine scheme.

This name can be specified via `set_name!`.

# Examples
```jldoctest
julia> X = affine_space(QQ, 3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> name(X)
"unnamed affine variety"

julia> set_name!(X, "affine 3-dimensional space")

julia> name(X)
"affine 3-dimensional space"
```
"""
@attr String function name(X::AbsSpec)
  return "unnamed affine variety"
end


function set_name!(X::AbsSpec, name::String)
  return set_attribute!(X, :name, name)
end

#############################################################################
# (3) Attributes of AbsSpec
#     reduced scheme and singular locus
#############################################################################
# TODO: projective schemes, covered schemes

@Markdown.doc """
   reduced_scheme(X::AbsSpec{<:Field, <:MPolyAnyRing})

Return the induced reduced scheme of `X`.

Currently, this command is available for affine schemes and space germs.
 
This command relies on [`radical`](@ref).

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> J = ideal(R,[(x-y)^2])
ideal(x^2 - 2*x*y + y^2)

julia> X = Spec(R,J)
Spec of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x^2 - 2*x*y + y^2)

julia> U = MPolyComplementOfKPointIdeal(R,[0,0])
complement of maximal ideal corresponding to point with coordinates fmpq[0, 0]

julia> Y = Spec(R,J,U)
Spec of Localization of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x^2 - 2*x*y + y^2) at the multiplicative set complement of maximal ideal corresponding to point with coordinates fmpq[0, 0]

julia> reduced_scheme(X)
Spec of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x - y)

julia> reduced_scheme(Y)
Spec of Localization of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x - y) at the multiplicative set complement of maximal ideal corresponding to point with coordinates fmpq[0, 0]

```
"""

@attr function reduced_scheme(X::AbsSpec{<:Field, <:MPolyQuoLocalizedRing})
  I = modulus(OO(X))
  J = radical(pre_saturated_ideal(I))
  inc = ClosedEmbedding(X, ideal(OO(X), OO(X).(gens(J))))
  return domain(inc), inc
  return Spec(base_ring(J), J, inverted_set(OO(X)))
end

@attr function reduced_scheme(X::AbsSpec{<:Field, <:MPolyQuo})
  J = radical(modulus(OO(X)))
  inc = ClosedEmbedding(X, ideal(OO(X), OO(X).(gens(J))))
  return domain(inc), inc
  return Spec(base_ring(J), J)
end

## to make reduced_scheme agnostic for quotient ring
@attr function reduced_scheme(X::AbsSpec{<:Field, <:MPAnyNonQuoRing})
  return X, ClosedEmbedding(X, ideal(OO(X), one(OO(X))))
end

function reduced_scheme(X::AbsSpec)
  error("method 'reduced_scheme(X)' currently only implemented for affine Schemes and Space Germs over a field")
end

### TODO: The following two functions (singular_locus, 
###       singular_locus_reduced) need to be made type-sensitive 
###       and reduced=true needs to be set automatically for varieties 
###       as soon as not only schemes, but also varieties as separte
###       type have been introduced in OSCAR
### TODO: Make singular locus also available for projective schemes and
###       for covered schemes (using the workhorse here...).
 
@doc Markdown.doc"""
    singular_locus(X::Scheme{<:Field}) -> (Scheme, SchemeMor)

Return the singular locus of `X`.

For computing the singular locus of the reduced scheme induced by `X`,
please use [`singular_locus_reduced`](@ref).

Currently this command is available for affine schemes and for space germs.

Over non-perfect fields, this command returns the non-smooth locus and `X`
may still be regular at some points of the returned subscheme.

See also [`is_smooth`](@ref).

# Examples
``` jldoctest
julia> R, (x,y,z) = QQ["x", "y", "z"]
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> I = ideal(R, [x^2 - y^2 + z^2])
ideal(x^2 - y^2 + z^2)

julia> A3 = Spec(R)
Spec of Multivariate Polynomial Ring in x, y, z over Rational Field

julia> X = Spec(R,I)
Spec of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(x^2 - y^2 + z^2)

julia> singular_locus(A3)
(Spec of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(1), morphism from

	Spec of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(1)

to

	Spec of Multivariate Polynomial Ring in x, y, z over Rational Field

with coordinates

	0, 0, 0)

julia> singular_locus(X)
(Spec of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(x^2 - y^2 + z^2, z, y, x), morphism from

	Spec of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(x^2 - y^2 + z^2, z, y, x)

to

	Spec of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(x^2 - y^2 + z^2)

with coordinates

	0, 0, 0)

julia> U = MPolyComplementOfKPointIdeal(R,[0,0,0])
complement of maximal ideal corresponding to point with coordinates fmpq[0, 0, 0]

julia> Y = Spec(R,I,U)
Spec of Localization of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(x^2 - y^2 + z^2) at the multiplicative set complement of maximal ideal corresponding to point with coordinates fmpq[0, 0, 0]

julia> singular_locus(Y)
(Spec of Localization of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(x^2 - y^2 + z^2, z, y, x, x^2 - y^2 + z^2) at the multiplicative set complement of maximal ideal corresponding to point with coordinates fmpq[0, 0, 0], morphism from

	Spec of Localization of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(x^2 - y^2 + z^2, z, y, x, x^2 - y^2 + z^2) at the multiplicative set complement of maximal ideal corresponding to point with coordinates fmpq[0, 0, 0]

to

	Spec of Localization of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(x^2 - y^2 + z^2) at the multiplicative set complement of maximal ideal corresponding to point with coordinates fmpq[0, 0, 0]

with coordinates

	x, y, z)

```
"""
function singular_locus(X::AbsSpec{<:Field, <:MPAnyQuoRing})
  comp = _singular_locus_with_decomposition(X,false)
  if length(comp) == 0
    set_attribute!(X, :is_smooth, true)
    inc = ClosedEmbedding(X, ideal(OO(X), one(OO(X))))
    return domain(inc), inc
    return subscheme(X, ideal(OO(X),one(OO(X))))
  end
  R = base_ring(OO(X))
  I = prod([modulus(underlying_quotient(OO(Y))) for Y in comp])
  I = radical(I)
  set_attribute!(X, :is_smooth, false)
  inc = ClosedEmbedding(X, ideal(OO(X), OO(X).(gens(I))))
  return domain(inc), inc
end

# make singular_locus agnostic to quotient
function singular_locus(X::AbsSpec{<:Field, <:MPAnyNonQuoRing})
  set_attribute!(X, :is_smooth,true)
  inc = ClosedEmbedding(X, ideal(OO(X), one(OO(X))))
  return domain(inc), inc
end

# TODO: Covered schemes, projective schemes

@doc Markdown.doc"""
    singular_locus_reduced(X::Scheme{<:Field}) -> (Scheme, SchemeMor)

Return the singular locus of the reduced scheme ``X_{red}`` induced by `X`.

For computing the singular locus of `X` itself, please use 
['singular_locus](@ref)'.

Currently this command is available for affine schemes and for space germs.

Over non-perfect fields, this command returns the non-smooth locus and
`X_{red}` may still be regular at some points of the returned subscheme.

See also [`is_smooth`](@ref).

# Examples
``` jldoctest
julia> R, (x,y,z) = QQ["x", "y", "z"]
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> I = ideal(R, [(x^2 - y^2 + z^2)^2])
ideal(x^4 - 2*x^2*y^2 + 2*x^2*z^2 + y^4 - 2*y^2*z^2 + z^4)

julia> X = Spec(R,I)
Spec of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(x^4 - 2*x^2*y^2 + 2*x^2*z^2 + y^4 - 2*y^2*z^2 + z^4)

julia> singular_locus_reduced(X)
(Spec of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(x^4 - 2*x^2*y^2 + 2*x^2*z^2 + y^4 - 2*y^2*z^2 + z^4, z, y, x), morphism from

	Spec of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(x^4 - 2*x^2*y^2 + 2*x^2*z^2 + y^4 - 2*y^2*z^2 + z^4, z, y, x)

to

	Spec of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(x^4 - 2*x^2*y^2 + 2*x^2*z^2 + y^4 - 2*y^2*z^2 + z^4)

with coordinates

	0, 0, 0)

julia> singular_locus(X)
(Spec of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(x^4 - 2*x^2*y^2 + 2*x^2*z^2 + y^4 - 2*y^2*z^2 + z^4, x^2 - y^2 + z^2), morphism from

	Spec of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(x^4 - 2*x^2*y^2 + 2*x^2*z^2 + y^4 - 2*y^2*z^2 + z^4, x^2 - y^2 + z^2)

to

	Spec of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(x^4 - 2*x^2*y^2 + 2*x^2*z^2 + y^4 - 2*y^2*z^2 + z^4)

with coordinates

	x, y, z)

```
"""
function singular_locus_reduced(X::AbsSpec{<:Field, <:MPAnyQuoRing})
  comp =  _singular_locus_with_decomposition(X, true)
  I= ideal(ambient_coordinate_ring(X),one(ambient_coordinate_ring(X)))
  for Z in comp
    I = intersect(I, modulus(underlying_quotient(OO(Z))))
    # TODO: Already compute intermediate radicals?
  end
  I = radical(I)
  inc = ClosedEmbedding(X, ideal(OO(X), OO(X).(gens(I))))
  return domain(inc), inc
end

# make singular_locus_reduced agnostic to quotient
function singular_locus_reduced(X::AbsSpec{<:Field, <:MPAnyNonQuoRing})
  inc = ClosedEmbedding(X, ideal(OO(X), one(OO(X))))
  return domain(inc), inc
end

# internal workhorse, not user-facing
function _singular_locus_with_decomposition(X::AbsSpec{<:Field, <:MPAnyQuoRing}, reduced::Bool=true)
  I = saturated_ideal(modulus(OO(X)))
  empty = typeof(X)[]
  result = empty

# equidimensional decompositon to allow Jacobi criterion on each component
  P = Ideal[]

  if has_attribute(X, :is_equidimensional) && is_equidimensional(X) && !reduced 
    push!(P, I)
  else 
    if reduced
      P = equidimensional_decomposition_radical(I)
    else
      P = equidimensional_decomposition_weak(I)
    end
  end

# if irreducible, just do Jacobi criterion
  if length(P)==1 && !reduced
    d = dim(X)
    R = base_ring(I)
    n = nvars(R) 
    M = _jacobi_matrix_modulus(X)
    minvec = minors(M, n-d)
    J = ideal(R, minvec)
    JX = ideal(OO(X),minvec)
    one(OO(X)) in JX && return empty
    return [subscheme(X, J)]
  else
# if reducible, determine pairwise intersection loci
    components = [subscheme(X, J) for J in P]
    for i in 1:length(components)
      for j in (i+1):length(components)
        W = intersect(components[i],components[j])
        if !isempty(W) 
          push!(result, W)
        end
      end
    end
# and singular loci of components
    for Y in components
      result = vcat(result, singular_locus(Y)[1])
    end
  one(OO(X)) in result && return empty
  end
  return result
end

## cheaper version of jacobi_matrix specifically for Jacobi matrix of modulus
## compute *some* representative of the jacobian matrix of gens(modulus),
## forgetting about the denominators (contribution killed by modulus anyway)

function _jacobi_matrix_modulus(X::AbsSpec{<:Ring, <:MPAnyQuoRing})
  g = gens(modulus(underlying_quotient(OO(X))))
  L = base_ring(underlying_quotient(OO(X)))
  n = nvars(L)
  M = matrix(L, n, length(g),[derivative(f,i) for i=1:n for f in g])
  return M
end

########################################################################
# (X) Attributes for AbsSpec   to be deleted
#     
########################################################################

# TODO: ambient_closure_ideal should be deleted

@Markdown.doc """
    ambient_closure_ideal(X::AbsSpec{<:Any, <:MPolyRing})

Return the defining ideal of the closure of ``X`` in its ambient affine space.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> R = OO(X)
Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> (x1,x2,x3) = gens(R)
3-element Vector{fmpq_mpoly}:
 x1
 x2
 x3

julia> Y = subscheme(X, ideal(R, [x1*x2]))
Spec of Quotient of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field by ideal(x1*x2)

julia> I = Oscar.ambient_closure_ideal(Y)
ideal(x1*x2)

julia> base_ring(I) == OO(Y)
false

julia> base_ring(I) == R
true
```
"""
@attr ambient_closure_ideal(X::AbsSpec{<:Any, <:MPolyRing}) = ideal(OO(X), [zero(OO(X))])
ambient_closure_ideal(X::AbsSpec{<:Any, <:MPolyQuo}) = modulus(OO(X))
@attr ambient_closure_ideal(X::AbsSpec{<:Any, <:MPolyLocalizedRing}) = ideal(ambient_coordinate_ring(X), [zero(ambient_coordinate_ring(X))])
ambient_closure_ideal(X::AbsSpec{<:Any, <:MPolyQuoLocalizedRing}) = saturated_ideal(modulus(OO(X)))


########################################################################
# (4) Implementation of the AbsSpec interface for the basic Spec
########################################################################

OO(X::Spec) = X.OO
base_ring(X::Spec) = X.kk
ambient_coordinate_ring(X::Spec{<:Any, <:MPolyRing}) = OO(X)
ambient_coordinate_ring(X::Spec{<:Any, <:MPolyQuo}) = base_ring(OO(X))
ambient_coordinate_ring(X::Spec{<:Any, <:MPolyLocalizedRing}) = base_ring(OO(X))
ambient_coordinate_ring(X::Spec{<:Any, <:MPolyQuoLocalizedRing}) = base_ring(OO(X))
ambient_coordinate_ring(X::Spec{T, T}) where {T<:Field} = base_ring(X)



########################################################################
# (5) Type getters
########################################################################

# TODO: Needed?

ring_type(::Type{SpecType}) where {BRT, RT, SpecType<:AbsSpec{BRT, RT}} = RT
ring_type(X::AbsSpec) = ring_type(typeof(X))

base_ring_type(::Type{SpecType}) where {BRT, RT, SpecType<:AbsSpec{BRT, RT}} = BRT
base_ring_type(X::AbsSpec) = base_ring_type(typeof(X))
base_ring_elem_type(::Type{SpecType}) where {BRT, RT, SpecType<:AbsSpec{BRT, RT}} = elem_type(BRT)
base_ring_elem_type(X::AbsSpec) = base_ring_elem_type(typeof(X))

poly_type(::Type{SpecType}) where {BRT, RT<:MPolyRing, SpecType<:AbsSpec{BRT, RT}} = elem_type(RT)
poly_type(::Type{SpecType}) where {BRT, T, RT<:MPolyQuo{T}, SpecType<:AbsSpec{BRT, RT}} = T
poly_type(::Type{SpecType}) where {BRT, T, RT<:MPolyLocalizedRing{<:Any, <:Any, <:Any, T}, SpecType<:AbsSpec{BRT, RT}} = T
poly_type(::Type{SpecType}) where {BRT, T, RT<:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, T}, SpecType<:AbsSpec{BRT, RT}} = T
poly_type(X::AbsSpec) = poly_type(typeof(X))

ring_type(::Type{Spec{BRT, RT}}) where {BRT, RT} = RT
ring_type(X::Spec) = ring_type(typeof(X))
base_ring_type(::Type{Spec{BRT, RT}}) where {BRT, RT} = BRT
base_ring_type(X::Spec) = base_ring_type(typeof(X))

