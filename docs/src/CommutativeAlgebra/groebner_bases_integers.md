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

# Gröbner and Standard Bases Over $\mathbb Z$

Over the integers the coefficients of the polynomials 
are not invertible, thus their handling when computing
Gröbner bases and normal forms plays an important role. This is done when 
computing strong Gröbner bases which ensure the following property: 
For any element of an ideal its leading term is divisible by a leading term of an 
element of a corresponding strong Gröbner basis.

The textbook [AL94](@cite) provides details on theory and algorithms as well as references.

```jldoctest
julia> R, (x,y) = PolynomialRing(ZZ, ["x","y"])

julia> I = ideal(R, [2x,3x,4y])

julia> H = groebner_basis(I)

```
