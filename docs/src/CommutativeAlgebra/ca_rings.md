```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["ca_rings.md"]
```

# Multivariate Polynomial Rings


## Constructors

```julia
PolynomialRing(P::Ring, v::Vector{AbstractString}; ordering=:lex)
```

!!! note
    Consider fields and ZZ as coefficient ring `P`. Which fields are allowed


### Examples
```@repl oscar
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
f = x^2+y*z
typeof(f)
S, (a,b) = PolynomialRing(ZZ, ["a", "b"])
g = a*b
typeof(g)
```

!!! note
    Add docu for convenient constructor using indexed variables.





Return a tuple, say `R, vars`, consisting of a polynomial ring `R` and an array
`vars` of generators (variables) which print according to the strings in the supplied
vector `v` . The ordering can at present be `:lex`, `:deglex` or `:degrevlex`.
If no ordering is specified, `ordering=:lex` is chosen by default.


```@docs
    grade(R::MPolyRing)
```
