#############################################################################
##
##  common actions of group elements
##
##  The idea is to delegate the action of `GAPGroupElem` objects
##  on `GapObj` objects to the corresponding GAP action,
##  and to implement the action on native Julia objects case by case.

"""
We try to avoid introducing `on_points` and `on_right`.
Note that the GAP functions `OnPoints` and `OnRight` just delegate
to powering `^` and right multiplication `*`, respectively.
Thus we have to make sure that `^` and `*` are installed in all
relevant situations.
One such case is the action of `GAPGroupElem` objects on `GapObj`
objects, for example wrapped GAP matrices on GAP vectors:

```jldoctest
julia> g = GL(2,3);

julia> m = g[1]
[2   0]
[0   1]

julia> v = GapObj(m)[1]
GAP: [ Z(3), 0*Z(3) ]

julia> v^m
GAP: [ Z(3)^0, 0*Z(3) ]

julia> v*m
GAP: [ Z(3)^0, 0*Z(3) ]
```
"""

^(pnt::GAP.Obj, x::GAPGroupElem) = GAP.Globals.:^(pnt, GapObj(x))

*(pnt::GAP.Obj, x::GAPGroupElem) = GAP.Globals.:*(pnt, GapObj(x))


"""
    on_tuples(tuple::GapObj, x::GAPGroupElem)
    on_tuples(tuple::Vector, x::GAPGroupElem)
    on_tuples(tuple::T, x::GAPGroupElem) where T <: Tuple

Return the image of `tuple` under `x`,
where the action is given by applying `^` to the entries of `tuple`.

For `Vector` and `Tuple` objects,
one can also call `^` instead of `on_tuples`.

# Examples
```jldoctest
julia> g = symmetric_group(3);  g[1]
(1,2,3)

julia> l = GapObj([1, 2, 4])
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
on_tuples(tuple::GapObj, x::GAPGroupElem) = GAPWrap.OnTuples(tuple, GapObj(x))

on_tuples(tuple::Vector{T}, x::GroupElem) where {T} = T[pnt^x for pnt in tuple]
^(tuple::Vector{T}, x::GroupElem) where {T} = on_tuples(tuple, x)

on_tuples(tuple::T, x::GroupElem) where {T <: Tuple} = T(pnt^x for pnt in tuple)
^(tuple::T, x::GroupElem) where {T <: Tuple} = on_tuples(tuple, x)


"""
    on_sets(set::GapObj, x::GAPGroupElem)
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

julia> l = GapObj([1, 3])
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
on_sets(set::GapObj, x::GAPGroupElem) = GAPWrap.OnSets(set, GapObj(x))

function on_sets(set::Vector{T}, x::GroupElem) where {T}
    res = T[pnt^x for pnt in set]
    sort!(res)
    return res
end

on_sets(set::T, x::GroupElem) where {T<:AbstractSet} = T(pnt^x for pnt in set)

function on_sets(set::T, x::GroupElem) where {T<:Tuple}
    res = [pnt^x for pnt in set]
    sort!(res)
    return T(res)
end

^(set::AbstractSet, x::GroupElem) = on_sets(set, x)

"""
    on_sets_sets(set::GapObj, x::GAPGroupElem)
    on_sets_sets(set::Vector, x::GAPGroupElem)
    on_sets_sets(set::Tuple, x::GAPGroupElem)
    on_sets_sets(set::AbstractSet, x::GAPGroupElem)

Return the image of `set` under `x`,
where the action is given by applying `on_sets` to the entries
of `set`, and then turning the result into a sorted vector/tuple or a set,
respectively.

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> g = symmetric_group(3);  g[1]
(1,2,3)

julia> l = GapObj([[1, 2], [3, 4]]; recursive = true)
GAP: [ [ 1, 2 ], [ 3, 4 ] ]

julia> on_sets_sets(l, g[1])
GAP: [ [ 1, 4 ], [ 2, 3 ] ]

julia> on_sets_sets([[1, 2], [3, 4]], g[1])
2-element Vector{Vector{Int64}}:
 [1, 4]
 [2, 3]

julia> on_sets_sets(((1, 2), (3, 4)), g[1])
((1, 4), (2, 3))

julia> on_sets_sets(Set([[1, 2], [3, 4]]), g[1])
Set{Vector{Int64}} with 2 elements:
  [2, 3]
  [1, 4]

julia> setset = Set([BitSet([1, 2]), BitSet([3, 4])]);

julia> on_sets_sets(setset, g[1])
Set{BitSet} with 2 elements:
  BitSet([1, 4])
  BitSet([2, 3])

julia> ans == setset^g[1]
true
```
"""
on_sets_sets(set::GapObj, x::GAPGroupElem) = GAPWrap.OnSetsSets(set, GapObj(x))

function on_sets_sets(set::Vector{T}, x::GroupElem) where {T}
    res = T[on_sets(pnt, x) for pnt in set]
    sort!(res)
    return res
end

on_sets_sets(set::T, x::GroupElem) where {T<:AbstractSet} =
  T(on_sets(pnt, x) for pnt in set)

function on_sets_sets(set::T, x::GroupElem) where {T<:Tuple}
    res = [on_sets(pnt, x) for pnt in set]
    sort!(res)
    return T(res)
end


"""
    permuted(pnt::GapObj, x::PermGroupElem)
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

julia> l = GapObj(a; recursive = true)
GAP: [ "a", "b", "c" ]

julia> permuted(l, g[1])
GAP: [ "c", "a", "b" ]
```
"""
permuted(pnt::GapObj, x::PermGroupElem) = GAPWrap.Permuted(pnt, GapObj(x))

function permuted(pnt::Vector{T}, x::PermGroupElem) where T
   invx = inv(x)
   return pnt[[i^invx for i in 1:length(pnt)]]
end

function permuted(pnt::T, x::PermGroupElem) where T <: Tuple
   invx = inv(x)
   return T(pnt[[i^invx for i in 1:length(pnt)]])
end


@doc raw"""
    on_indeterminates(f::GapObj, p::PermGroupElem)
    on_indeterminates(f::MPolyRingElem, p::PermGroupElem)
    on_indeterminates(f::FreeAssociativeAlgebraElem, p::PermGroupElem)
    on_indeterminates(f::MPolyIdeal, p::PermGroupElem)

Return the image of `f` under `p` where `p` acts via permuting the indeterminates.

For `MPolyRingElem`, `FreeAssociativeAlgebraElem`, and `MPolyIdeal` objects,
one can also call `^` instead of `on_indeterminates`.

# Examples
```jldoctest
julia> g = symmetric_group(3);  p = g[1]
(1,2,3)

julia> R, x = polynomial_ring(QQ, [:x1, :x2, :x3]);

julia> f = x[1]*x[2] + x[2]*x[3]
x1*x2 + x2*x3

julia> f^p
x1*x3 + x2*x3

julia> x = [GAP.Globals.X(GAP.Globals.Rationals, i) for i in 1:3];

julia> f = x[1]*x[2] + x[2]*x[3]
GAP: x_1*x_2+x_2*x_3

julia> on_indeterminates(f, p)
GAP: x_1*x_3+x_2*x_3
```
"""
on_indeterminates(f::GapObj, p::PermGroupElem) = GAPWrap.OnIndeterminates(f, GapObj(p))

function on_indeterminates(f::MPolyRingElem, s::PermGroupElem)
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

function on_indeterminates(f::FreeAssociativeAlgebraElem{T}, s::PermGroupElem) where T
  G = parent(s)
  S = parent(f)
  @assert ngens(S) == degree(G)

  e = exponent_words(f)
  c = Vector{T}(collect(coefficients(f)))
  es = [on_tuples(x, s) for x in e]
  return S(c, es)
end

@doc raw"""
    on_indeterminates(f::GapObj, p::MatrixGroupElem)
    on_indeterminates(f::MPolyRingElem{T}, p::MatrixGroupElem{T}) where T
    on_indeterminates(f::MPolyIdeal, p::MatrixGroupElem)

Return the image of `f` under `p` where `p` acts via evaluating `f` at the
vector obtained by multiplying `p` with the (column) vector of indeterminates.
This corresponds to considering the variables of the polynomial ring containing
`f` as the basis of a vector space on which `p` acts by multiplication from the
right.

For `MPolyRingElem` and `MPolyIdeal` objects, one can also call `^` instead of
`on_indeterminates`.

# Examples
```jldoctest
julia> g = general_linear_group(2, 5);  m = g[2]
[4   1]
[4   0]

julia> R, x = polynomial_ring(base_ring(g), degree(g));

julia> f = x[1]*x[2] + x[1]
x1*x2 + x1

julia> f^m
x1^2 + 4*x1*x2 + 4*x1 + x2
```
"""
function on_indeterminates(f::GapObj, p::MatrixGroupElem)
  # We assume that we act on the indeterminates with numbers 1, ..., nrows(p).
  # (Note that `f` does not know about a polynomial ring to which it belongs.)
  n = nrows(p)
  fam = GAPWrap.CoefficientsFamily(GAPWrap.FamilyObj(f))
  indets = GapObj([GAPWrap.Indeterminate(fam, i) for i in 1:n])
  return GAPWrap.Value(f, indets, GapObj(p) * indets)
end

function on_indeterminates(f::MPolyRingElem{T}, p::MatrixGroupElem{T}) where T
  @assert base_ring(f) == base_ring(p)
  @assert ngens(parent(f)) == degree(parent(p))
  act = right_action(parent(f), p)
  return act(f)
end

function on_indeterminates(I::MPolyIdeal, p::PermGroupElem)
  @assert ngens(base_ring(I)) == degree(parent(p))
  imggens = [on_indeterminates(x, p) for x in gens(I)]
  return ideal(parent(imggens[1]), imggens)
end

function on_indeterminates(I::MPolyIdeal, p::MatrixGroupElem)
  @assert base_ring(gen(I, 1)) == base_ring(p)
  @assert ngens(base_ring(I)) == degree(parent(p))
  imggens = [on_indeterminates(x, p) for x in gens(I)]
  return ideal(parent(imggens[1]), imggens)
end

^(f::MPolyRingElem, p::PermGroupElem) = on_indeterminates(f, p)

^(f::FreeAssociativeAlgebraElem, p::PermGroupElem) = on_indeterminates(f, p)

^(f::MPolyRingElem{T}, p::MatrixGroupElem{T, S}) where T where S = on_indeterminates(f, p)

^(I::MPolyIdeal, p::PermGroupElem) = on_indeterminates(I, p)

^(I::MPolyIdeal, p::MatrixGroupElem) = on_indeterminates(I, p)


# We do not support `on_lines(line::Vector, x::GAPGroupElem)`
# because for example `*(::fpFieldElem, ::Vector{fpFieldElem})`
# is not supported.
"""
    on_lines(line::GapObj, x::GAPGroupElem)
    on_lines(line::AbstractAlgebra.Generic.FreeModuleElem, x::GAPGroupElem)

Return the image of the nonzero vector `line` under `x`,
where the action is given by first computing `line * x`
and then normalizing the result by scalar multiplication from the left
such that the first nonzero entry is the `one` of the `base_ring` of `line`.

# Examples
```jldoctest
julia> n = 2;  F = GF(5);  g = general_linear_group(n, F);

julia> v = gen(free_module(F, n), 1)
(1, 0)

julia> m = gen(g, 2)
[4   1]
[4   0]

julia> v * m
(4, 1)

julia> on_lines(v, m)
(1, 4)
```
"""
on_lines(line::GapObj, x::GAPGroupElem) = GAPWrap.OnLines(line, GapObj(x))

function on_lines(line::AbstractAlgebra.Generic.FreeModuleElem, x::GAPGroupElem)
    res = line * x
    @assert ! iszero(res) "vector must be nonzero"
    i = 1
    while iszero(res[i])
      i = i+1
    end
    return inv(res[i]) * res
end

@doc raw"""
    on_subgroups(x::GapObj, g::GAPGroupElem) -> GapObj
    on_subgroups(x::T, g::GAPGroupElem) where T <: GAPGroup -> T

Return the image of the group `x` under `g`. Note that `x` must
be a subgroup of the domain of `g`.

# Examples
```jldoctest
julia> C = cyclic_group(20)
Pc group of order 20

julia> S = automorphism_group(C)
Automorphism group of
  pc group of order 20

julia> H, _ = sub(C, [gens(C)[1]^4])
(Sub-pc group of order 5, Hom: H -> C)

julia> all(g -> on_subgroups(H, g) == H, S)
true
```
"""
function on_subgroups(x::GapObj, g::GAPGroupElem)
  return GAPWrap.Image(GapObj(g), x)
end

on_subgroups(x::T, g::GAPGroupElem) where T <: GAPGroup = T(on_subgroups(GapObj(x), g))


@doc raw"""
    on_echelon_form_mats(m::MatElem{T}, x::MatrixGroupElem) where T <: FinFieldElem

Return the image of `m` under `x`,
where the action is given by first computing the product `m * x`
and then normalizing the result by computing its reduced row echelon form
with `echelon_form`.

Identifying `m` with the subspace of the natural module for the group of `x`
that is generated by the rows of `m`,
`on_echelon_form_mats` describes the action on subspaces of this natural module.
Note that for computing the orbit and stabilizer of `m` w.r.t.
`on_echelon_form_mats`, `m` must be in reduced row echelon form.

# Examples
```jldoctest
julia> n = 3;  q = 2;  F = GF(q);  V = free_module(F, n);

julia> G = GL(n, F);

julia> W, embW = sub(V, [gen(V,1), gen(V,3)])
(Subspace over F with 2 generators and no relations, Hom: W -> V)

julia> m = matrix(embW)
[1   0   0]
[0   0   1]

julia> S, _ = stabilizer(G, m, on_echelon_form_mats);  order(S)
24

julia> orb = orbit(G, on_echelon_form_mats, m);  length(orb)
7
```
"""
function on_echelon_form_mats(m::MatElem{T}, x::MatrixGroupElem) where T <: FinFieldElem
  return echelon_form(m * x)
end

@doc raw"""
    induced_action(actfun::Function, phi::GAPGroupHomomorphism)

Return the action function that is obtained by inducing `actfun` along `phi`.

That means, given a groups ``G`` and ``H``, a set ``\Omega`` with action function ``f: \Omega \times G \to \Omega``
and a homomorphism ``\phi: H \to G``, construct the action function
$\Omega \times H \to \Omega, (\omega, h) \mapsto f(\omega, \phi(h))$.
"""
function induced_action(actfun::Function, phi::GAPGroupHomomorphism)
  return _induced_action(actfun, phi)
end

# This method is not documented as we need `phi` to be a group homomorphism, but in many cases
# there is no dedicated type for this (WeylGroup, FinGenAbGroup, etc.).
# This should be restricted to group homomorphisms once we have a type for them.
function induced_action(actfun::Function, phi::Map{<:Union{Group,FinGenAbGroup}, <:Union{Group,FinGenAbGroup}})
  return _induced_action(actfun, phi)
end

function _induced_action(actfun::Function, phi::Map{<:Union{Group,FinGenAbGroup}, <:Union{Group,FinGenAbGroup}})
  return function (omega, g)
    return actfun(omega, phi(g))
  end
end

@doc raw"""
    stabilizer(G::GAPGroup, pnt::Any[, actfun::Function])

Return `S, emb` where `S` is the subgroup of `G` that consists of
all those elements `g` that fix `pnt` under the action given by `actfun`,
that is, `actfun(pnt, g) == pnt` holds,
and `emb` is the embedding of `S` into `G`.

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

julia> S = stabilizer(G, [1, 1, 2, 2, 3], permuted);  order(S[1])
4
```
"""
stabilizer(G::GAPGroup, pnt::Any, actfun::Function) = _stabilizer_generic(G, pnt, actfun)

function _stabilizer_generic(G::GAPGroup, pnt::Any, actfun::Function)
    return Oscar._as_subgroup(G, GAPWrap.Stabilizer(GapObj(G), pnt,
        GapObj(gens(G), recursive = true), GapObj(gens(G)),
        GapObj(actfun)))
end

# natural stabilizers in permutation groups
# Construct the arguments on the GAP side such that GAP's method selection
# can choose the special method.
# - stabilizer in a perm. group of an integer via `^`
# - stabilizer in a perm. group of a vector of integers via `on_tuples`
# - stabilizer in a perm. group of a set of integers via `on_sets`
function stabilizer(G::PermGroup, pnt::T) where T <: IntegerUnion
    return Oscar._as_subgroup(G, GAPWrap.Stabilizer(GapObj(G),
        GapObj(pnt),
        GAP.Globals.OnPoints))  # Do not use GAPWrap.OnPoints!
end

function stabilizer(G::PermGroup, pnt::Union{Vector{T}, Tuple{T, Vararg{T}}}) where T <: Oscar.IntegerUnion
    return Oscar._as_subgroup(G, GAPWrap.Stabilizer(GapObj(G),
        GapObj(pnt, recursive = true),
        GAP.Globals.OnTuples))  # Do not use GAPWrap.OnTuples!
end

function stabilizer(G::PermGroup, pnt::AbstractSet{T}) where T <: Oscar.IntegerUnion
    return Oscar._as_subgroup(G, GAPWrap.Stabilizer(GapObj(G),
        GapObj(pnt, recursive = true),
        GAP.Globals.OnSets))  # Do not use GAPWrap.OnSets!
end

# now the same with given action function,
# these calls may come from delegations from G-sets
function stabilizer(G::PermGroup, pnt::T, actfun::Function) where T <: IntegerUnion
    return (actfun == ^) ? stabilizer(G, pnt) : _stabilizer_generic(G, pnt, actfun)
end

function stabilizer(G::PermGroup, pnt::Union{Vector{T},Tuple{T,Vararg{T}}}, actfun::Function) where T <: IntegerUnion
    return actfun == on_tuples ? stabilizer(G, pnt) : _stabilizer_generic(G, pnt, actfun)
end

function stabilizer(G::PermGroup, pnt::AbstractSet{T}, actfun::Function) where T <: IntegerUnion
    return actfun == on_sets ? stabilizer(G, pnt) : _stabilizer_generic(G, pnt, actfun)
end

# natural stabilizers in matrix groups
# Construct the arguments on the GAP side such that GAP's method selection
# can choose the special method.
# - stabilizer in a matrix group (over a finite field)
#   of a `FreeModuleElem` via `*` (or `^`)
# - stabilizer in a matrix group (over a finite field)
#   of a vector or tuple of `FreeModuleElem`s via `on_tuples`
# - stabilizer in a matrix group (over a finite field)
#   of a `Set` of `FreeModuleElem`s via `on_sets`
function stabilizer(G::MatrixGroup{ET,MT}, pnt::AbstractAlgebra.Generic.FreeModuleElem{ET}) where {ET,MT}
    iso = Oscar.iso_oscar_gap(base_ring(parent(pnt)))
    return Oscar._as_subgroup(G, GAPWrap.Stabilizer(GapObj(G),
        map_entries(iso, AbstractAlgebra.Generic._matrix(pnt)),
        GAP.Globals.OnRight))
end

function stabilizer(G::MatrixGroup{ET,MT}, pnt::Vector{AbstractAlgebra.Generic.FreeModuleElem{ET}}) where {ET,MT}
    length(pnt) == 0 && return G
    iso = Oscar.iso_oscar_gap(base_ring(parent(pnt[1])))
    return Oscar._as_subgroup(G, GAPWrap.Stabilizer(GapObj(G),
        GapObj([GapObj(map_entries(iso, AbstractAlgebra.Generic._matrix(v)))[1] for v in pnt]),
        GAP.Globals.OnTuples))
end

function stabilizer(G::MatrixGroup{ET,MT}, pnt::Tuple{AbstractAlgebra.Generic.FreeModuleElem{ET},Vararg{AbstractAlgebra.Generic.FreeModuleElem{ET}}}) where {ET,MT}
    length(pnt) == 0 && return G
    iso = Oscar.iso_oscar_gap(base_ring(parent(pnt[1])))
    return Oscar._as_subgroup(G, GAPWrap.Stabilizer(GapObj(G),
        GapObj([GapObj(map_entries(iso, AbstractAlgebra.Generic._matrix(v)))[1] for v in pnt]),
        GAP.Globals.OnTuples))
end

function stabilizer(G::MatrixGroup{ET,MT}, pnt::AbstractSet{AbstractAlgebra.Generic.FreeModuleElem{ET}}) where {ET,MT}
    length(pnt) == 0 && return G
    iso = Oscar.iso_oscar_gap(base_ring(parent(iterate(pnt)[1])))
    return Oscar._as_subgroup(G, GAPWrap.Stabilizer(GapObj(G),
        GAPWrap.Set(GapObj([GapObj(map_entries(iso, AbstractAlgebra.Generic._matrix(v)))[1] for v in pnt])),
        GAP.Globals.OnSets))
end

# now the same with given action function,
# these calls may come from delegations from G-sets
function stabilizer(G::MatrixGroup{ET,MT}, pnt::AbstractAlgebra.Generic.FreeModuleElem{ET}, actfun::Function) where {ET,MT}
    return (actfun == *) ? stabilizer(G, pnt) : _stabilizer_generic(G, pnt, actfun)
end

function stabilizer(G::MatrixGroup{ET,MT}, pnt::Vector{AbstractAlgebra.Generic.FreeModuleElem{ET}}, actfun::Function) where {ET,MT}
    return (actfun == on_tuples) ? stabilizer(G, pnt) : _stabilizer_generic(G, pnt, actfun)
end

function stabilizer(G::MatrixGroup{ET,MT}, pnt::Tuple{AbstractAlgebra.Generic.FreeModuleElem{ET},Vararg{AbstractAlgebra.Generic.FreeModuleElem{ET}}}, actfun::Function) where {ET,MT}
    return (actfun == on_tuples) ? stabilizer(G, pnt) : _stabilizer_generic(G, pnt, actfun)
end

function stabilizer(G::MatrixGroup{ET,MT}, pnt::AbstractSet{AbstractAlgebra.Generic.FreeModuleElem{ET}}, actfun::Function) where {ET,MT}
    return (actfun == on_sets) ? stabilizer(G, pnt) : _stabilizer_generic(G, pnt, actfun)
end

# stabilizer in a matrix group (over a finite field)
# of a row reduced matrix via `on_echelon_form_mats`
function stabilizer(G::MatrixGroup{ET,<:MT}, pnt::MatElem{<:MT}, actfun::Function) where {ET,MT}
    (actfun === on_echelon_form_mats) || return _stabilizer_generic(G, pnt, actfun)
    nrows(pnt) == 0 && return (G, identity_map(G))
    iso = Oscar.iso_oscar_gap(base_ring(pnt))
    return Oscar._as_subgroup(G, GAPWrap.Stabilizer(GapObj(G),
        map_entries(iso, pnt),
        GAP.Globals.OnSubspacesByCanonicalBasis))
end

"""
    right_coset_action(G::GAPGroup, U::GAPGroup)

Compute the action of `G` on the right cosets of its subgroup `U`.

# Examples
```jldoctest
julia> G = symmetric_group(6);

julia> H = sylow_subgroup(G, 2)[1]
Permutation group of degree 6 and order 16

julia> index(G, H)
45

julia> act = right_coset_action(G, H);

julia> degree(codomain(act)) == index(G, H)
true
```
"""
function right_coset_action(G::GAPGroup, U::GAPGroup)
  _check_compatible(G, U)
  mp = GAPWrap.FactorCosetAction(GapObj(G), GapObj(U))
  @req mp !== GAP.Globals.fail "Invalid input"
  H = PermGroup(GAPWrap.Range(mp))
  return GAPGroupHomomorphism(G, H, mp)
end
