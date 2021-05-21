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
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> f = x^2+y*z
x^2 + y*z

julia> typeof(f)
fmpq_mpoly

julia> S, (a,b) = PolynomialRing(ZZ, ["a", "b"])
(Multivariate Polynomial Ring in a, b over Integer Ring, fmpz_mpoly[a, b])

julia> g = a*b
a*b

julia> typeof(g)
fmpz_mpoly
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
