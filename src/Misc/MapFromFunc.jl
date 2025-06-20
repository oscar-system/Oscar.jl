
"""
    map_from_func(D, C, f, [g])

Creates the map `D -> C, x -> f(x)` given the callable
object `f`. If `g` is provided, it is assumed to satisfy
`f(g(x)) = x` and will be used as the preimage function.

# Example

```jldoctest
julia> F = GF(2);

julia> f = map_from_func(QQ, F, x -> F(numerator(x)) * inv(F(denominator(x))))
Map defined by a julia-function
  from rational field
  to prime field of characteristic 2

julia> f(QQ(1//3))
1

julia> println(f)
Map: QQ -> F

julia> f = map_from_func(QQ, F, x -> F(numerator(x)) * inv(F(denominator(x))), y -> QQ(lift(ZZ, y)),)
Map defined by a julia-function with inverse
  from rational field
  to prime field of characteristic 2

julia> preimage(f, F(1))
1

julia> println(f)
Map: QQ -> F

```
"""
map_from_func(D, C, f) = MapFromFunc(D, C, f)
map_from_func(D, C, f, g) = MapFromFunc(D, C, f, g)
