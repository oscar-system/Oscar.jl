"""
A *group action* of a group G on a set Ω (from the right) is defined by
a map μ: Ω × G → Ω that satisfies the compatibility conditions
μ(μ(x, g), h) = μ(x, g*h) and μ(x, one(G)) == x for all x ∈ Ω.

The maps μ are implemented as functions that take two arguments, an element
x of Ω and a group element g, and return the image of x under g.

In many cases, a natural action is given by the types of the elements in Ω
and in G.
For example permutation groups act on positive integers by just applying
the permutations.
In such situations, the function `^` can be used as action function,
and `^` is taken as the default whenever no other function is prescribed.

However, the action is not always determined by the types of the involved
objects.
For example, permutations can act on arrays of positive integers by
applying the permutations pointwise, or by permuting the entries;
matrices can act on vectors by multiplying the vector with the matrix,
or by multiplying the inverse of the matrix with the vector;
and of course one can construct new custom actions in situations where
default actions are already available.

Thus it is in general necessary to specify the action function explicitly.
The following ones are commonly used.
"""

#############################################################################
##
##  common actions of group elements
##
##  The idea is to delegate the action of `GAPGroupElem` objects
##  on `GAP.GapObj` objects to the corresponding GAP action,
##  and to implement the action on native Julia objects case by case.

import Base.:^
import Base.:*

export on_tuples, on_sets, on_indeterminates, permuted

"""
We try to avoid introducing `on_points` and `on_right`.
Note that the GAP functions `OnPoints` and `OnRight` just delegate
to powering `^` and right multiplication `*`, respectively.
Thus we have to make sure that `^` and `*` are installed in all
relevant situations.
One such case is the action of `GAPGroupElem` objects on `GapObj`
objects, for example wrapped GAP matrices on GAP vectors:

```
julia> g = GL(2,3);

julia> m = g[1]
[ [ Z(3), 0*Z(3) ], [ 0*Z(3), Z(3)^0 ] ]

julia> v = m.X[1]
GAP: [ Z(3), 0*Z(3) ]

julia> v^m
GAP: [ Z(3)^0, 0*Z(3) ]

julia> v*m
GAP: [ Z(3)^0, 0*Z(3) ]

```
"""

^(pnt::GAP.GapObj, x::GAPGroupElem) = GAP.Globals.:^(pnt, x.X)

*(pnt::GAP.GapObj, x::GAPGroupElem) = GAP.Globals.:*(pnt, x.X)


"""
    on_tuples(tuple::GAP.GapObj, x::GAPGroupElem)
    on_tuples(tuple::Vector{T}, x::GAPGroupElem) where T
    on_tuples(tuple::T, x::GAPGroupElem) where T <: Tuple

Return the image of `tuple` under `x`,
where the action is given by applying `^` to the entries of `tuple`.

For `Vector` and `Tuple` objects,
one can also call `^` instead of `on_tuples`.

# Examples
```jldoctest
julia> g = symmetric_group(3);  g[1]
(1,2,3)

julia> l = GAP.julia_to_gap([1, 2, 4])
GAP: [ 1, 2, 4 ]

julia> on_tuples(l, g[1])
GAP: [ 2, 3, 4 ]

julia> on_tuples([1, 2, 4], g[1])
3-element Array{Int64,1}:
 2
 3
 4

julia> on_tuples((1, 2, 4), g[1])
(2, 3, 4)

```
"""
on_tuples(tuple::GAP.GapObj, x::GAPGroupElem) = GAP.Globals.OnTuples(tuple, x.X)

on_tuples(tuple::Vector{T}, x::GAPGroupElem) where T = T[pnt^x for pnt in tuple]
^(tuple::Vector{T}, x::GAPGroupElem) where T = on_tuples(tuple, x)

on_tuples(tuple::T, x::GAPGroupElem) where T <: Tuple = T([pnt^x for pnt in tuple])
^(tuple::T, x::GAPGroupElem) where T <: Tuple = on_tuples(tuple, x)


"""
    on_sets(set::GAP.GapObj, x::GAPGroupElem)
    on_sets(set::Vector{T}, x::GAPGroupElem) where T
    on_sets(set::T, x::GAPGroupElem) where T <: Union{Tuple, Set}

Return the image of `set` under `x`,
where the action is given by applying `^` to the entries
of `set`, and then turning the result into a sorted array/tuple or a set,
respectively.

For `Set` objects, one can also call `^` instead of `on_sets`.

# Examples
```jldoctest
julia> g = symmetric_group(3);  g[1]
(1,2,3)

julia> l = GAP.julia_to_gap([1,3])
GAP: [ 1, 3 ]

julia> on_sets(l, g[1])
GAP: [ 1, 2 ]

julia> on_sets([1, 3], g[1])
2-element Array{Int64,1}:
 1
 2

julia> on_sets((1, 3), g[1])
(1, 2)

julia> on_sets(Set([1, 3]), g[1])
Set{Int64} with 2 elements:
  2
  1

```
"""
on_sets(set::GAP.GapObj, x::GAPGroupElem) = GAP.Globals.OnSets(set, x.X)

function on_sets(set::Vector{T}, x::GAPGroupElem) where T
    res = T[pnt^x for pnt in set]
    sort!(res)
    return res
end

function on_sets(set::T, x::GAPGroupElem) where T <: Union{Tuple, Set}
    res = [pnt^x for pnt in set]
    sort!(res)
    return T(res)
end

^(set::T, x::GAPGroupElem) where T <: Set = on_sets(set, x)


"""
    permuted(pnt::GAP.GapObj, x::PermGroupElem)
    permuted(pnt::Vector{T}, x::PermGroupElem) where T
    permuted(pnt::T, x::PermGroupElem) where T <: Tuple

Return the image of `pnt` under `x`,
where the action is given by permuting the entries of `pnt` with `x`.

# Examples
```jldoctest
julia> g = symmetric_group(3);  g[1]
(1,2,3)

julia> a = ["a", "b", "c"]
3-element Array{String,1}:
 "a"
 "b"
 "c"

julia> permuted(a, g[1])
3-element Array{String,1}:
 "c"
 "a"
 "b"

julia> permuted(("a", "b", "c"), g[1])
("c", "a", "b")

julia> l = GAP.julia_to_gap(a, recursive = true)
GAP: [ "a", "b", "c" ]

julia> permuted(l, g[1])
GAP: [ "c", "a", "b" ]

```
"""
permuted(pnt::GAP.GapObj, x::PermGroupElem) = GAP.Globals.Permuted(pnt, x.X)

function permuted(pnt::Vector{T}, x::PermGroupElem) where T
   invx = inv(x)
   return pnt[[i^invx for i in 1:length(pnt)]]
end

function permuted(pnt::T, x::PermGroupElem) where T <: Tuple
   invx = inv(x)
   return T(pnt[[i^invx for i in 1:length(pnt)]])
end


"""
    on_indeterminates(f::GAP.GapObj, p::PermGroupElem)
    on_indeterminates(f::Nemo.MPolyElem{T}, p::PermGroupElem) where T

Returns the image of `f` under `p`, w.r.t. permuting the indeterminates
with `p`.

For `Nemo.MPolyElem{T}` objects, one can also call `^` instead of
`on_indeterminates`.

# Examples
```
julia> g = symmetric_group(3);  p = g[1]
(1,2,3)

julia> R, x = PolynomialRing(QQ, ["x1", "x2", "x3"]);

julia> f = x[1]*x[2] + x[2]*x[3]
x1*x2 + x2*x3

julia> f^p
x1*x3 + x2*x3

julia> x = [GAP.Globals.X( GAP.Globals.Rationals, i ) for i in 1:3];

julia> f = x[1]*x[2] + x[2]*x[3]
GAP: x_1*x_2+x_2*x_3

julia> on_indeterminates(f, p)
GAP: x_1*x_3+x_2*x_3

```
"""
on_indeterminates(f::GAP.GapObj, p::PermGroupElem) = GAP.Globals.OnIndeterminates(f, p.X)

function on_indeterminates(f::Nemo.MPolyElem{T}, p::PermGroupElem) where T
    pnt = gens(parent(f))
    return evaluate(f, pnt[[i^p for i in 1:length(pnt)]])
end

^(f::Nemo.MPolyElem{T}, p::PermGroupElem) where T = on_indeterminates(f, p)
