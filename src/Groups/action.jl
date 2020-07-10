#############################################################################
##
##  common actions of group elements
##
##  The idea is to delegate the action of `GAPGroupElem` objects
##  on `GAP.GapObj` objects to the corresponding GAP action,
##  and to implement the action on native Julia objects case by case.

export on_points, on_right, on_tuples, on_sets, permuted


"""
    on_points(pnt::GAP.GapObj, x::GAPGroupElem)
    on_points(pnt::GAPGroupElem, x::GAPGroupElem)
    on_points(pnt::Int, x::PermGroupElem)

Return the image of `pnt` under `x`.
If `pnt` is a `GAP.GapObj` then the action is given by `GAP.Globals.OnPoints`.
If `pnt` is a `GAPGroupElem` then the action is given by conjugation,
that is, the result is `x`^-1 * `pnt` * `x`.
If `pnt` is a positive integer and `x` is a permutation then the action
is given by the natural action that maps `pnt` to `x(pnt)`.

# Examples
```jldoctest
julia> g = symmetric_group(3);  g[1]
(1,2,3)

julia> g[2]
(1,2)

julia> on_points(g[2], g[1])
(2,3)

julia> on_points(g[2].X, g[1])
GAP: (2,3)

julia> on_points(1, g[1])
2

```
"""
on_points(pnt::GAP.GapObj, x::GAPGroupElem) = GAP.Globals.OnPoints(pnt, x.X)

on_points(pnt::GAPGroupElem, x::GAPGroupElem) = inv(x) * pnt * x

on_points(pnt::Int, x::PermGroupElem) = GAP.Globals.OnPoints(pnt, x.X)


"""
    on_right(pnt::GAP.GapObj, x::GAPGroupElem)
    on_right(pnt::GAPGroupElem, x::GAPGroupElem)
    on_right(v::Vector{T}, x::MatrixGroupElem) where T

Return the image of the point `pnt` under `x`,
where the action is given by right multiplication with `x`.

# Examples
```jldoctest
julia> g = general_linear_group(2,3);  mat = g[1]
[ [ Z(3), 0*Z(3) ], [ 0*Z(3), Z(3)^0 ] ]

julia> on_right(mat, mat)
[ [ Z(3)^0, 0*Z(3) ], [ 0*Z(3), Z(3)^0 ] ]

julia> r = mat.X[1]
GAP: [ Z(3), 0*Z(3) ]

julia> on_right( r, mat )
GAP: [ Z(3)^0, 0*Z(3) ]

julia> v = Vector{GAP.FFE}(r)
2-element Array{FFE,1}:
 GAP: Z(3)
 GAP: 0*Z(3)

julia> on_right(v, mat)
2-element Array{GAP.FFE,1}:
 GAP: Z(3)^0
 GAP: 0*Z(3)

```
"""
on_right(pnt::GAP.GapObj, x::GAPGroupElem) = GAP.Globals.OnRight(pnt, x.X)

on_right(pnt::GAPGroupElem, x::GAPGroupElem) = pnt * x

function on_right(v::Vector{T}, x::MatrixGroupElem) where T
    gapv = GAP.julia_to_gap(v)
    img = GAP.Globals.OnRight(gapv, x.X)
    return Vector{T}(img)
end


"""
    on_tuples(tuple::GAP.GapObj, x::GAPGroupElem)
    on_tuples(tuple::Vector{T}, x::GAPGroupElem) where T
    on_tuples(tuple::T, x::GAPGroupElem) where T <: Tuple

Return the image of `tuple` under `x`,
where the action is given by applying [`on_points`](@ref) to the entries
of `tuple`.

# Examples
```jldoctest
julia> g = symmetric_group(3);  g[1]
(1,2,3)

julia> l = GAP.julia_to_gap([1, 2, 3, 4, 5])
GAP: [ 1, 2, 3, 4, 5 ]

julia> on_tuples(l, g[1])
GAP: [ 2, 3, 1, 4, 5 ]

julia> on_tuples([1, 2, 3, 4, 5], g[1])
5-element Array{Int64,1}:
 2
 3
 1
 4
 5

julia> on_tuples((1, 2, 3, 4, 5), g[1])
(2, 3, 1, 4, 5)

```
"""
on_tuples(tuple::GAP.GapObj, x::GAPGroupElem) = GAP.Globals.OnTuples(tuple, x.X)

on_tuples(tuple::Vector{T}, x::GAPGroupElem) where T = T[on_points(pnt, x) for pnt in tuple]

on_tuples(tuple::T, x::GAPGroupElem) where T <: Tuple = T([on_points(pnt, x) for pnt in tuple])


"""
    on_sets(set::GAP.GapObj, x::GAPGroupElem)
    on_sets(set::Vector{T}, x::GAPGroupElem) where T
    on_sets(set::T, x::GAPGroupElem) where T <: Union{Tuple, Set}

Return the image of `set` under `x`,
where the action is given by applying [`on_points`](@ref) to the entries
of `set`, and then turning the result into a sorted array/tuple or a set,
respectively.

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
    res = T[on_points(pnt, x) for pnt in set]
    sort!(res)
    return res
end

function on_sets(set::T, x::GAPGroupElem) where T <: Union{Tuple, Set}
    res = [on_points(pnt, x) for pnt in set]
    sort!(res)
    return T(res)
end


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

julia> a = ["a", "b", "c", "d", "e"]
5-element Array{String,1}:
 "a"
 "b"
 "c"
 "d"
 "e"

julia> l = GAP.julia_to_gap(a, recursive = true)
GAP: [ "a", "b", "c", "d", "e" ]

julia> permuted(l, g[1])
GAP: [ "c", "a", "b", "d", "e" ]

julia> permuted(a, g[1])
5-element Array{String,1}:
 "c"
 "a"
 "b"
 "d"
 "e"

julia> permuted(("a", "b", "c", "d", "e"), g[1])
("c", "a", "b", "d", "e")

```
"""
permuted(pnt::GAP.GapObj, x::PermGroupElem) = GAP.Globals.Permuted(pnt, x.X)

function permuted(pnt::Vector{T}, x::PermGroupElem) where T
   invx = inv(x)
   return Vector{T}(pnt[[invx(i) for i in 1:length(pnt)]])
end

function permuted(pnt::T, x::PermGroupElem) where T <: Tuple
   invx = inv(x)
   return T(pnt[[invx(i) for i in 1:length(pnt)]])
end

