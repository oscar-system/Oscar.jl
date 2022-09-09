export on_indeterminates,
       on_sets,
       on_sets_sets,
       on_tuples,
       permuted,
       right_coset_action

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
For example, permutations can act on vectors of positive integers by
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

^(pnt::GAP.Obj, x::GAPGroupElem) = GAP.Globals.:^(pnt, x.X)

*(pnt::GAP.Obj, x::GAPGroupElem) = GAP.Globals.:*(pnt, x.X)


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

julia> l = GAP.GapObj([1, 2, 4])
GAP: [ 1, 2, 4 ]

julia> on_tuples(l, g[1])
GAP: [ 2, 3, 4 ]

julia> on_tuples([1, 2, 4], g[1])
3-element Vector{Int64}:
 2
 3
 4

julia> on_tuples((1, 2, 4), g[1])
(2, 3, 4)

julia> (1, 2, 4)^g[1]
(2, 3, 4)

```
"""
on_tuples(tuple::GAP.GapObj, x::GAPGroupElem) = GAP.Globals.OnTuples(tuple, x.X)

on_tuples(tuple::Vector{T}, x::GAPGroupElem) where T = T[pnt^x for pnt in tuple]
^(tuple::Vector{T}, x::GAPGroupElem) where T = on_tuples(tuple, x)

on_tuples(tuple::T, x::GAPGroupElem) where T <: Tuple = T(pnt^x for pnt in tuple)
^(tuple::T, x::GAPGroupElem) where T <: Tuple = on_tuples(tuple, x)


"""
    on_sets(set::GAP.GapObj, x::GAPGroupElem)
    on_sets(set::Vector, x::GAPGroupElem)
    on_sets(set::Tuple, x::GAPGroupElem)
    on_sets(set::AbstractSet, x::GAPGroupElem)

Return the image of `set` under `x`,
where the action is given by applying `^` to the entries
of `set`, and then turning the result into a sorted vector/tuple or a set,
respectively.

For `Set` objects, one can also call `^` instead of `on_sets`.

# Examples
```jldoctest
julia> g = symmetric_group(3);  g[1]
(1,2,3)

julia> l = GAP.GapObj([1,3])
GAP: [ 1, 3 ]

julia> on_sets(l, g[1])
GAP: [ 1, 2 ]

julia> on_sets([1, 3], g[1])
2-element Vector{Int64}:
 1
 2

julia> on_sets((1, 3), g[1])
(1, 2)

julia> on_sets(Set([1, 3]), g[1])
Set{Int64} with 2 elements:
  2
  1

julia> BitSet([1, 3])^g[1]
BitSet with 2 elements:
  1
  2

```
"""
on_sets(set::GAP.GapObj, x::GAPGroupElem) = GAP.Globals.OnSets(set, x.X)

function on_sets(set::Vector{T}, x::GAPGroupElem) where T
    res = T[pnt^x for pnt in set]
    sort!(res)
    return res
end

on_sets(set::T, x::GAPGroupElem) where T <: AbstractSet = T(pnt^x for pnt in set)

function on_sets(set::T, x::GAPGroupElem) where T <: Tuple
    res = [pnt^x for pnt in set]
    sort!(res)
    return T(res)
end

^(set::AbstractSet, x::GAPGroupElem) = on_sets(set, x)

"""
    on_sets_sets(set::GAP.GapObj, x::GAPGroupElem)
    on_sets_sets(set::Vector, x::GAPGroupElem)
    on_sets_sets(set::Tuple, x::GAPGroupElem)
    on_sets_sets(set::AbstractSet, x::GAPGroupElem)

Return the image of `set` under `x`,
where the action is given by applying `on_sets` to the entries
of `set`, and then turning the result into a sorted vector/tuple or a set,
respectively.

# Examples
```jldoctest
julia> g = symmetric_group(3);  g[1]
(1,2,3)

julia> l = GAP.GapObj([[1, 2], [3, 4]], recursive = true)
GAP: [ [ 1, 2 ], [ 3, 4 ] ]

julia> on_sets_sets(l, g[1])
GAP: [ [ 1, 4 ], [ 2, 3 ] ]

julia> on_sets_sets([[1, 2], [3, 4]], g[1])
2-element Vector{Vector{Int64}}:
 [1, 4]
 [2, 3]

julia> on_sets_sets(((1,2), (3,4)), g[1])
((1, 4), (2, 3))

julia> on_sets_sets(Set([[1, 2], [3, 4]]), g[1])
Set{Vector{Int64}} with 2 elements:
  [1, 4]
  [2, 3]

julia> setset = Set([BitSet([1, 2]), BitSet([3, 4])]);

julia> on_sets_sets(setset, g[1])
Set{BitSet} with 2 elements:
  BitSet([1, 4])
  BitSet([2, 3])

julia> ans == setset^g[1]
true

```
"""
on_sets_sets(set::GAP.GapObj, x::GAPGroupElem) = GAP.Globals.OnSetsSets(set, x.X)

function on_sets_sets(set::Vector{T}, x::GAPGroupElem) where T
    res = T[on_sets(pnt, x) for pnt in set]
    sort!(res)
    return res
end

on_sets_sets(set::T, x::GAPGroupElem) where T <: AbstractSet = T(on_sets(pnt, x) for pnt in set)

function on_sets_sets(set::T, x::GAPGroupElem) where T <: Tuple
    res = [on_sets(pnt, x) for pnt in set]
    sort!(res)
    return T(res)
end


"""
    permuted(pnt::GAP.GapObj, x::PermGroupElem)
    permuted(pnt::Vector, x::PermGroupElem)
    permuted(pnt::Tuple, x::PermGroupElem)

Return the image of `pnt` under `x`,
where the action is given by permuting the entries of `pnt` with `x`.

# Examples
```jldoctest
julia> g = symmetric_group(3);  g[1]
(1,2,3)

julia> a = ["a", "b", "c"]
3-element Vector{String}:
 "a"
 "b"
 "c"

julia> permuted(a, g[1])
3-element Vector{String}:
 "c"
 "a"
 "b"

julia> permuted(("a", "b", "c"), g[1])
("c", "a", "b")

julia> l = GAP.GapObj(a, recursive = true)
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
    on_indeterminates(f::Nemo.MPolyElem, p::PermGroupElem)
    on_indeterminates(f::MPolyIdeal, p::PermGroupElem)
    on_indeterminates(f::GAP.GapObj, p::MatrixGroupElem)
    on_indeterminates(f::Nemo.MPolyElem{T}, p::MatrixGroupElem{T, S}) where T where S
    on_indeterminates(f::MPolyIdeal, p::MatrixGroupElem)

Return the image of `f` under `p`.
If `p` is a `PermGroupElem` then it acts via permuting the indeterminates,
if `p` is a `MatrixGroupElem` then it acts via evaluating `f` at the
vector obtained by multiplying `p` with the (column) vector of indeterminates.

For `Nemo.MPolyElem` and `MPolyIdeal` objects,
one can also call `^` instead of `on_indeterminates`.

# Examples
```jldoctest
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

julia> g = general_linear_group(2, 5);  m = g[2]
[4   1]
[4   0]

julia> R, x = PolynomialRing(base_ring(g), degree(g));

julia> f = x[1]*x[2] + x[1]
x1*x2 + x1

julia> f^m
x1^2 + 4*x1*x2 + 4*x1 + x2

```
"""
on_indeterminates(f::GAP.GapObj, p::PermGroupElem) = GAP.Globals.OnIndeterminates(f, p.X)

function on_indeterminates(f::Nemo.MPolyElem, s::PermGroupElem)
  G = parent(s)
  @assert ngens(parent(f)) == degree(G)

  g = Generic.MPolyBuildCtx(parent(f))
  for (c, e) = Base.Iterators.zip(Generic.MPolyCoeffs(f), Generic.MPolyExponentVectors(f))
    s_e = zeros(Int, degree(G))
    for i=1:degree(G)
      s_e[s(i)] = e[i]
    end
    push_term!(g, c, s_e)
  end
  return finish(g)
end

function on_indeterminates(f::GAP.GapObj, p::MatrixGroupElem)
  # We assume that we act on the indeterminates with numbers 1, ..., nrows(p).
  # (Note that `f` does not know about a polynomial ring to which it belongs.)
  n = nrows(p)
  fam = GAP.Globals.CoefficientsFamily(GAP.Globals.FamilyObj(f))
  indets = GAP.GapObj([GAP.Globals.Indeterminate(fam, i) for i in 1:n])
  return GAP.Globals.Value(f, indets, p.X * indets)
end

function on_indeterminates(f::Nemo.MPolyElem{T}, p::MatrixGroupElem{T, S}) where T where S
  @assert base_ring(f) == base_ring(p)
  @assert ngens(parent(f)) == degree(parent(p))
  act = Oscar.right_action(parent(f), p)
  return act(f)
end

function on_indeterminates(I::MPolyIdeal, p::PermGroupElem)
  @assert ngens(parent(gen(I, 1))) == degree(parent(p))
  imggens = [on_indeterminates(x, p) for x in gens(I)]
  return ideal(parent(imggens[1]), imggens)
end

function on_indeterminates(I::MPolyIdeal, p::MatrixGroupElem)
  @assert base_ring(gen(I, 1)) == base_ring(p)
  @assert ngens(parent(gen(I, 1))) == degree(parent(p))
  imggens = [on_indeterminates(x, p) for x in gens(I)]
  return ideal(parent(imggens[1]), imggens)
end

^(f::Nemo.MPolyElem, p::PermGroupElem) = on_indeterminates(f, p)

^(f::Nemo.MPolyElem{T}, p::MatrixGroupElem{T, S}) where T where S = on_indeterminates(f, p)

^(I::MPolyIdeal, p::PermGroupElem) = on_indeterminates(I, p)

^(I::MPolyIdeal, p::MatrixGroupElem) = on_indeterminates(I, p)


@doc Markdown.doc"""
    stabilizer(G::Oscar.GAPGroup, pnt::Any[, actfun::Function])

Return the subgroup of `G` that consists of all those elements `g`
that fix `pnt` under the action given by `actfun`,
that is, `actfun(pnt, g) == pnt` holds.

The default for `actfun` depends on the types of `G` and `pnt`:
If `G` is a `PermGroup` then the default actions on integers,
`Vector`s of  integers, and `Set`s of integers are given by
`^`, `on_tuples`, and `on_sets`, respectively.
If `G` is a `MatrixGroup` then the default actions on `FreeModuleElem`s,
`Vector`s of them, and `Set`s of them are given by
`*`, `on_tuples`, and `on_sets`, respectively.

# Examples
```jldoctest
julia> G = symmetric_group(5);

julia> S = stabilizer(G, 1);  order(S[1])
24

julia> S = stabilizer(G, [1, 2]);  order(S[1])
6

julia> S = stabilizer(G, Set([1, 2]));  order(S[1])
12

julia> S = stabilizer(G, [1,1,2,2,3], permuted);  order(S[1])
4

```
"""
function stabilizer(G::Oscar.GAPGroup, pnt::Any, actfun::Function)
    return Oscar._as_subgroup(G, GAP.Globals.Stabilizer(G.X, pnt,
        GAP.GapObj([x.X for x in gens(G)]), GAP.GapObj(gens(G)),
        GAP.WrapJuliaFunc(actfun)))
end

# natural stabilizers in permutation groups
stabilizer(G::PermGroup, pnt::T) where T <: Oscar.IntegerUnion = stabilizer(G, pnt, ^)

stabilizer(G::PermGroup, pnt::Vector{T}) where T <: Oscar.IntegerUnion = stabilizer(G, pnt, on_tuples)

stabilizer(G::PermGroup, pnt::AbstractSet{T}) where T <: Oscar.IntegerUnion = stabilizer(G, pnt, on_sets)

# natural stabilizers in matrix groups
stabilizer(G::MatrixGroup{ET,MT}, pnt::AbstractAlgebra.Generic.FreeModuleElem{ET}) where {ET,MT} = stabilizer(G, pnt, *)

stabilizer(G::MatrixGroup{ET,MT}, pnt::Vector{AbstractAlgebra.Generic.FreeModuleElem{ET}}) where {ET,MT} = stabilizer(G, pnt, on_tuples)

stabilizer(G::MatrixGroup{ET,MT}, pnt::AbstractSet{AbstractAlgebra.Generic.FreeModuleElem{ET}}) where {ET,MT} = stabilizer(G, pnt, on_sets)


"""
    right_coset_action(G::T, U::T) where T <: GAPGroup

Compute the action of `G` on the right cosets of its subgroup `U`.

# Examples
```jldoctest
julia> G = symmetric_group(6);

julia> H = sylow_subgroup(G, 2)[1]
Group([ (1,2), (3,4), (1,3)(2,4), (5,6) ])

julia> index(G, H)
45

julia> act = right_coset_action(G, H);

julia> degree(codomain(act)) == index(G, H)
true

```
"""
function right_coset_action(G::T, U::T) where T <: GAPGroup
  mp = GAP.Globals.FactorCosetAction(G.X, U.X)
  mp == GAP.Globals.fail && throw(ArgumentError("Invalid input"))
  H = PermGroup(GAP.Globals.Range(mp))
  return GAPGroupHomomorphism(G, H, mp)
end
