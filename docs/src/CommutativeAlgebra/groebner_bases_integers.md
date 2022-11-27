```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["groebner_bases_integers.md"]
```

# Gröbner/Standard Bases Over $\mathbb Z$

Over the integers the coefficients of the polynomials 
are not invertible, thus their handling when computing
Gröbner bases and normal forms plays an important role. This is done when 
computing strong Gröbner bases which ensure the following property: 
For any element of an ideal its leading term is divisible by a leading term of an 
element of a corresponding strong Gröbner basis.

The textbook [AL94](@cite) provides details on theory and algorithms as well as references.

```jldoctest
julia> R, (x,y) = PolynomialRing(ZZ, ["x","y"])
(Multivariate Polynomial Ring in x, y over Integer Ring, fmpz_mpoly[x, y])

julia> I = ideal(R, [2x,3x,4y])
ideal(2*x, 3*x, 4*y)

julia> H = groebner_basis(I)
Gröbner basis with elements
1 -> 4*y
2 -> x
with respect to the ordering
degrevlex([x, y])

```
