# Unit groups of finite rings

## Aims

Implementation of the algorithms from [Hof26](@cite) for computing unit groups and the first K-group of finite rings and finite associative algebras.

## Examples

```jldoctest
julia> R, = finite_ring(GF(2)[symmetric_group(4)]); # this the modular group ring F_2[S_4]

julia> U, mU = unit_group(R)
(Finitely presented group of order 3145728, Map: U -> finite ring)

julia> R, = finite_ring(GF(2)[small_group(4, 2)]);

julia> Oscar.k1(R)
((Z/2)^3, Map: (Z/2)^3 -> finite ring)
```
