#
#
#

Base.isfinite(G::PermGroup) = true

==(x::PermGroup, y::PermGroup) = x.deg == y.deg && x.X == y.X

==(x::PermGroupElem, y::PermGroupElem) = degree(x) == degree(y) && x.X == y.X

Base.:<(x::PermGroupElem, y::PermGroupElem) = x.X < y.X

Base.isless(x::PermGroupElem, y::PermGroupElem) = x<y


@doc raw"""
    degree(G::PermGroup) -> Int

Return the degree of `G` as a permutation group, that is,
an integer `n` that is stored in `G`, with the following meaning.

- `G` embeds into `symmetric_group(n)`.
- Two permutation groups of different degrees are regarded as not equal,
  even if they contain the same permutations.
- Subgroups constructed with `derived_subgroup`, `sylow_subgroup`, etc.,
  get the same degree as the given group.
- The range `1:degree(G)` is used as the default set of points on which
  `G` and its element acts.

!!! note
    The degree of a group of permutations is not necessarily equal to the
    largest moved point of the group `G`. For example, the trivial subgroup of
    `symmetric_group(n)` has degree `n` even though it fixes `n`.

# Examples
```jldoctest
julia> degree(symmetric_group(4))
4

julia> t4 = trivial_subgroup(symmetric_group(4))[1];

julia> degree(t4)
4

julia> t4 == trivial_subgroup(symmetric_group(5))[1]
false

julia> show(Vector(gen(symmetric_group(4), 2)))
[2, 1, 3, 4]
julia> show(Vector(gen(symmetric_group(5), 2)))
[2, 1, 3, 4, 5]
```
"""
degree(x::PermGroup) = x.deg

@doc raw"""
    degree(g::PermGroupElem) -> Int

Return the degree of the parent of `g`.
This value is always greater or equal `number_of_moved_points(g)`

"""
degree(g::PermGroupElem) = degree(parent(g))

@doc raw"""
    moved_points(x::PermGroupElem) -> Vector{Int}
    moved_points(G::PermGroup) -> Vector{Int}

Return the vector of those points in `1:degree(x)` or `1:degree(G)`,
respectively, that are not mapped to themselves under the action `^`.

# Examples
```jldoctest
julia> g = symmetric_group(4);  s = sylow_subgroup(g, 3)[1];

julia> length(moved_points(s))
3

julia> length(moved_points(gen(s, 1)))
3
```
"""
@gapattribute moved_points(x::Union{PermGroupElem,PermGroup}) = Vector{Int}(GAP.Globals.MovedPoints(x.X))

@doc raw"""
    number_of_moved_points(x::PermGroupElem) -> Int
    number_of_moved_points(G::PermGroup) -> Int

Return the number of those points in `1:degree(x)` or `1:degree(G)`,
respectively, that are moved (i.e., not fixed) under the action `^`.

# Examples
```jldoctest
julia> g = symmetric_group(4);  s = sylow_subgroup(g, 3)[1];

julia> number_of_moved_points(s)
3

julia> number_of_moved_points(gen(s, 1))
3
```
"""
@gapattribute number_of_moved_points(x::Union{PermGroupElem,PermGroup}) = GAP.Globals.NrMovedPoints(x.X)::Int

@doc raw"""
    perm(L::AbstractVector{<:IntegerUnion})

Return the permutation $x$ which maps every $i$ from `1` to $n$` = length(L)`
to `L`$[i]$.
The parent of $x$ is set to [`symmetric_group`](@ref)$(n)$.
An exception is thrown if `L` does not contain every integer from 1 to $n$
exactly once.

The parent group of $x$ is set to [`symmetric_group`](@ref)$(n)$.

# Examples
```jldoctest
julia> x = perm([2,4,6,1,3,5])
(1,2,4)(3,6,5)

julia> parent(x)
Sym(6)
```
"""
function perm(L::AbstractVector{<:IntegerUnion})
  return PermGroupElem(symmetric_group(length(L)), GAPWrap.PermList(GAP.GapObj(L;recursive=true)))
end


@doc raw"""
    perm(G::PermGroup, L::AbstractVector{<:IntegerUnion})
    (G::PermGroup)(L::AbstractVector{<:IntegerUnion})

Return the permutation $x$ which maps every `i` from 1 to $n$` = length(L)`
to `L`$[i]$. The parent of $x$ is `G`.
An exception is thrown if $x$ is not contained in `G`
or `L` does not contain every integer from 1 to $n$ exactly once.

# Examples
```jldoctest
julia> perm(symmetric_group(6),[2,4,6,1,3,5])
(1,2,4)(3,6,5)
```

Equivalent permutations can be created using [`cperm`](@ref) and [`@perm`](@ref)
```jldoctest
julia> x = perm(symmetric_group(8),[2,3,1,5,4,7,8,6])
(1,2,3)(4,5)(6,7,8)

julia> y = cperm([1,2,3],[4,5],[6,7,8])
(1,2,3)(4,5)(6,7,8)

julia> x == y
true

julia> z = @perm (1,2,3)(4,5)(6,7,8)
(1,2,3)(4,5)(6,7,8)

julia> x == z
true
```
"""
function perm(g::PermGroup, L::AbstractVector{<:IntegerUnion})
   x = GAPWrap.PermList(GAP.GapObj(L;recursive=true))
   @req x !== GAP.Globals.fail "the list does not describe a permutation"
   @req (length(L) <= degree(g) && x in g.X) "the element does not embed in the group"
   return PermGroupElem(g, x)
end

perm(g::PermGroup, L::AbstractVector{<:ZZRingElem}) = perm(g, [Int(y) for y in L])

function (g::PermGroup)(L::AbstractVector{<:IntegerUnion})
   x = GAPWrap.PermList(GAP.GapObj(L;recursive=true))
   @req (length(L) <= degree(g) && x in g.X) "the element does not embed in the group"
   return PermGroupElem(g, x)
end

(g::PermGroup)(L::AbstractVector{<:ZZRingElem}) = g([Int(y) for y in L])

# cperm stands for "cycle permutation", but we can change name if we want
# takes as input a list of vectors (not necessarily disjoint)
@doc raw"""
    cperm(L::AbstractVector{<:T}...) where T <: IntegerUnion
    cperm(G::PermGroup, L::AbstractVector{<:T}...)
    cperm(L::Vector{Vector{T}}) where T <: IntegerUnion
    cperm(g::PermGroup,L::Vector{Vector{T}}) where T <: IntegerUnion

For given lists $[a_1, a_2, \ldots, a_n], [b_1, b_2, \ldots , b_m], \ldots$
of positive integers, return the
permutation $x = (a_1, a_2, \ldots, a_n) * (b_1, b_2, \ldots, b_m) * \ldots$.
Arrays of the form `[n, n+1, ..., n+k]` can be replaced by `n:n+k`.

The parent of $x$ is `G`.
If `G` is not specified then the parent of $x$ is set to
[`symmetric_group`](@ref)$(n)$,
where $n$ is the largest integer that occurs in an entry of `L`.

An exception is thrown if $x$ is not contained in `G`
or one of the given vectors is empty or contains duplicates.

# Examples
```jldoctest
julia> cperm([1,2,3],4:7)
(1,2,3)(4,5,6,7)

julia> cperm([1,2],[2,3])
(1,3,2)

julia> cperm()
()

julia> p = cperm([1,2,3],[7])
(1,2,3)

julia> degree(p)
7
```

Two permutations coincide if, and only if, they move the same points and their parent groups have the same degree.
```jldoctest
julia> G=symmetric_group(5);

julia> A=alternating_group(5);

julia> x=cperm(G,[1,2,3]);

julia> y=cperm(A,[1,2,3]);

julia> z=cperm([1,2,3]); parent(z)
Sym(3)

julia> x==y
true

julia> x==z
false
```
In the example above, `x` and `y` are equal because both act on a set of cardinality `5`, while `x` and `z` are different because `x` belongs to `Sym(5)` and `z` belongs to `Sym(3)`.

cperm can also handle cycles passed in inside of a vector
```jldoctest
julia> x = cperm([[1,2],[3,4]])
(1,2)(3,4)

julia> y = cperm([1,2],[3,4])
(1,2)(3,4)

julia> x == y
true
```

```jldoctest
julia> G=symmetric_group(5)
Sym(5)

julia> x = cperm(G,[[1,2],[3,4]])
(1,2)(3,4)

julia> parent(x)
Sym(5)
```

Equivalent permutations can be created using [`perm`](@ref) and [`@perm`](@ref):
```jldoctest
julia> x = cperm([1,2,3],[4,5],[6,7,8])
(1,2,3)(4,5)(6,7,8)

julia> y = perm(symmetric_group(8),[2,3,1,5,4,7,8,6])
(1,2,3)(4,5)(6,7,8)

julia> x == y
true

julia> z = @perm (1,2,3)(4,5)(6,7,8)
(1,2,3)(4,5)(6,7,8)

julia> x == z
true
```

At the moment, the input vectors of the function `cperm` need not be disjoint.

"""
function cperm()
  return one(symmetric_group(1))
end

function cperm(L1::AbstractVector{T}, L::AbstractVector{T}...) where T <: IntegerUnion
  return prod([PermGroupElem(symmetric_group(maximum(y)), GAPWrap.CycleFromList(GAP.Obj([Int(k) for k in y]))) for y in [L1, L...]])
  #TODO: better create the product of GAP permutations?
end

# cperm stays for "cycle permutation", but we can change name if we want
# takes as input a list of vectors (not necessarily disjoint)
# WARNING: we allow e.g. PermList([2,3,1,4,5,6]) in Sym(3)
function cperm(g::PermGroup)
  return one(g)
end

function cperm(g::PermGroup, L1::AbstractVector{T}, L::AbstractVector{T}...) where T <: IntegerUnion
  x = prod(y -> GAPWrap.CycleFromList(GAP.Obj([Int(k) for k in y])), [L1, L...])
  @req x in g.X "the element does not embed in the group"
  return PermGroupElem(g, x)
end

function cperm(L::Vector{Vector{T}}) where T <: IntegerUnion
    return cperm(L...)
end

function cperm(g::PermGroup,L::Vector{Vector{T}}) where T <: IntegerUnion
    return cperm(g,L...)
end

@doc raw"""
    Vector{T}(x::PermGroupElem, n::Int = x.parent.deg) where T <: IntegerUnion
    Vector(x::PermGroupElem, n::Int = x.parent.deg)

Return the list of length `n` that contains `x(i)` at position `i`. If not specified, `T` is set as `Int`.

# Examples
```jldoctest
julia> pi = cperm(1:3)
(1,2,3)
julia> Vector(pi)
3-element Vector{Int64}:
 2
 3
 1
julia> Vector(pi, 2)
2-element Vector{Int64}:
 2
 3
julia> Vector(pi, 4)
4-element Vector{Int64}:
 2
 3
 1
 4
julia> Vector{ZZRingElem}(pi, 2)
2-element Vector{ZZRingElem}:
 2
 3
```
"""
Base.Vector{T}(x::PermGroupElem, n::Int = x.parent.deg) where T <: IntegerUnion = T[x(i) for i in 1:n]
Base.Vector(x::PermGroupElem, n::Int = x.parent.deg) = Vector{Int}(x,n)

#evaluation function
(x::PermGroupElem)(n::IntegerUnion) = n^x

^(n::T, x::PermGroupElem) where T <: IntegerUnion = T(GAP.Obj(n)^x.X)

^(n::Int, x::PermGroupElem) = (n^x.X)::Int


@doc raw"""
    sign(g::PermGroupElem) -> Int

Return the sign of the permutation `g`.

The sign of a permutation ``g`` is defined as ``(-1)^k`` where ``k`` is the number of
cycles of ``g`` of even length.

# Examples
```jldoctest
julia> sign(cperm(1:2))
-1

julia> sign(cperm(1:3))
1
```
"""
Base.sign(g::PermGroupElem) = GAPWrap.SignPerm(g.X)

# TODO: document the following?
Base.sign(G::PermGroup) = GAPWrap.SignPermGroup(G.X)


@doc raw"""
    isodd(g::PermGroupElem)

Return `true` if the permutation `g` is odd, `false` otherwise.

A permutation is odd if it has an odd number of cycles of even length.
Equivalently, a permutation is odd if it has sign ``-1``.

# Examples
```jldoctest
julia> isodd(cperm(1:2))
true

julia> isodd(cperm(1:3))
false

julia> isodd(cperm(1:2,3:4))
false
```
"""
Base.isodd(g::PermGroupElem) = sign(g) == -1

@doc raw"""
    iseven(g::PermGroupElem)

Return `true` if the permutation `g` is even, `false` otherwise.

A permutation is even if it has an even number of cycles of even length.
Equivalently, a permutation is even if it has sign ``+1``.

# Examples
```jldoctest
julia> iseven(cperm(1:2))
false

julia> iseven(cperm(1:3))
true

julia> iseven(cperm(1:2,3:4))
true
```
"""
Base.iseven(n::PermGroupElem) = !isodd(n)


# TODO: document the following?
Base.isodd(G::PermGroup) = sign(G) == -1
Base.iseven(n::PermGroup) = !isodd(n)

##
# cycle types and support
##
struct CycleType <: AbstractVector{Pair{Int64, Int64}}
  # pairs 'cycle length => number of times it occurs'
  # so 'n => 1' is a single n-cycle and  '1 => n' is the identity on n points
  s::Vector{Pair{Int, Int}}

  # take a vector of cycle lengths
  function CycleType(c::Vector{Int})
    s = Vector{Pair{Int, Int}}()
    for i = c
      _push_cycle!(s, i)
    end
    sort!(s, by = x -> x[1])
    return new(s)
  end
  function CycleType(v::Vector{Pair{Int, Int}}; sorted::Bool = false)
    sorted && return new(v)
    return new(sort(v, by = x -> x[1]))
#TODO: check that each cycle length is specified at most once?
  end
end

Base.iterate(C::CycleType) = iterate(C.s)
Base.iterate(C::CycleType, x) = iterate(C.s, x)
Base.length(C::CycleType) = length(C.s)
Base.eltype(C::CycleType) = Pair{Int, Int}
Base.getindex(C::CycleType, i::Int) = C.s[i]
Base.size(C::CycleType) = size(C.s)

function Base.hash(c::CycleType, u::UInt = UInt(121324))
  return hash(c.s, u)
end

function Base.show(io::IO, C::CycleType)
  print(io, C.s)
end

function _push_cycle!(s::Vector{Pair{Int, Int}}, i::Int, j::Int = 1)
  # TODO: rewrite this to use searchsortedfirst instead of findfirst,
  # then avoid the sort! below
  f = findfirst(x->x[1] == i, s)
  if f === nothing
    push!(s, i=>j)
    sort!(s, by = x -> x[1])
  else
    s[f] = s[f][1]=>s[f][2] + j
  end
end

function ^(c::CycleType, e::Int)
  t = Vector{Pair{Int, Int}}()
  for (i,j) in c.s
    g = gcd(i, e)
    _push_cycle!(t, divexact(i, g), g*j)
  end
  return CycleType(t; sorted=true)
end


@doc raw"""
    order(::Type{T} = ZZRingElem, c::CycleType) where T <: IntegerUnion

Return the order of the permutations with cycle structure `c`.

# Examples
```jldoctest
julia> g = symmetric_group(3);

julia> all(x -> order(cycle_structure(x)) == order(x), gens(g))
true
```
"""
order(::Type{T}, c::CycleType) where T = mapreduce(x->T(x[1]), lcm, c.s, init = T(1))
order(c::CycleType) = order(ZZRingElem, c)


@doc raw"""
    degree(c::CycleType) -> Int

Return the degree of the permutations with cycle structure `c`.

# Examples
```jldoctest
julia> g = symmetric_group(3);

julia> all(x -> degree(cycle_structure(x)) == degree(g), gens(g))
true
```
"""
degree(c::CycleType) = mapreduce(x->x[1]*x[2], +, c.s, init = 0)


@doc raw"""
    sign(c::CycleType) -> Int

Return the sign of the permutations with cycle structure `c`.

# Examples
```jldoctest
julia> g = symmetric_group(3);

julia> all(x -> sign(cycle_structure(x)) == sign(x), gens(g))
true
```
"""
function Base.sign(c::CycleType)
    res = 1
    for (a, b) in c.s
      if iseven(a) && isodd(b)
        res = - res
      end
    end
    return res
end

@doc raw"""
    isodd(c::CycleType) -> Bool

Return whether the permutations with cycle structure `c` are odd.

# Examples
```jldoctest
julia> g = symmetric_group(3);

julia> all(x -> isodd(cycle_structure(x)) == isodd(x), gens(g))
true
```
"""
Base.isodd(c::CycleType) = sign(c) == -1


@doc raw"""
    iseven(c::CycleType) -> Bool

Return whether the permutations with cycle structure `c` are even.

# Examples
```jldoctest
julia> g = symmetric_group(3);

julia> all(x -> iseven(cycle_structure(x)) == iseven(x), gens(g))
true
```
"""
Base.iseven(c::CycleType) = !isodd(c)


@doc raw"""
    cycle_structure(g::PermGroupElem) -> CycleType

Return the cycle structure of the permutation `g` as a cycle type.
A cycle type behaves similar to a vector of pairs `k => n`
indicating that there are `n` cycles of length `k`.

# Examples
```jldoctest
julia> g = cperm(1:3, 4:5, 6:7, 8:10, 11:15)
(1,2,3)(4,5)(6,7)(8,9,10)(11,12,13,14,15)

julia> cycle_structure(g)
3-element Oscar.CycleType:
 2 => 2
 3 => 2
 5 => 1

julia> g = cperm()
()

julia> cycle_structure(g)
1-element Oscar.CycleType:
 1 => 1
```
"""
function cycle_structure(g::PermGroupElem)
    c = GAPWrap.CycleStructurePerm(g.X)
    # TODO: use SortedDict from DataStructures.jl ?
    ct = Pair{Int, Int}[ i+1 => c[i] for i in 1:length(c) if GAP.Globals.ISB_LIST(c, i) ]
    s = degree(CycleType(ct, sorted = true))
    if s < degree(g)
      @assert length(c) == 0 || ct[1][1] > 1
      insert!(ct, 1, 1=>degree(g)-s)
    end
    return CycleType(ct, sorted = true)
end

function cycle_structure(x::GroupConjClass{PermGroup, PermGroupElem})
  return cycle_structure(representative(x))
end


@doc raw"""
    cycle_structures(G::PermGroup) -> Set{CycleType}

Return the set of cycle structures of elements in `G`,
see [`cycle_structure`](@ref).

# Examples
```jldoctest
julia> g = symmetric_group(3);

julia> sort!(collect(cycle_structures(g)))
3-element Vector{Oscar.CycleType}:
 [1 => 1, 2 => 1]
 [1 => 3]
 [3 => 1]
```
"""
function cycle_structures(G::PermGroup)
  r = conjugacy_classes(G)
  return Set(cycle_structure(x) for x in r)
end


################################################################################
#
#   _perm_helper
#
# The following code implements a new way to input permutations in Julia. For example
# it is possible to create a permutation as follow
# pi = Oscar.Permutations.@perm (1,2,3)(4,5)(6,7,8)
# > (1,2,3)(4,5)(6,7,8)
# For this we use macros to modify the syntax tree of (1,2,3)(4,5)(6,7,8) such that
# Julia can deal with the expression.

function _perm_helper(ex::Expr)

    ex == :( () ) && return []
    ex isa Expr || error("Input is not a permutation expression")

    res = []
    while ex isa Expr && ex.head == :call
        push!(res, Expr(:vect, ex.args[2:end]...))
        ex = ex.args[1]
    end

    if !(ex isa Expr) || ex.head != :tuple
        error("Input is not a permutation.")
    end

    push!(res, Expr(:vect,ex.args...))

    # reverse `res` to match the original order; this ensures
    # the evaluation order is as the user expects
    reverse!(res)

    return res
end


################################################################################
#
#   perm
#
@doc raw"""
    @perm ex
    
Input a permutation in cycle notation. Supports arbitrary expressions for
generating the integer entries of the cycles. The parent group is inferred 
to be the symmetric group with a degree of the highest integer referenced 
in the permutation.

The actual work is done by [`cperm`](@ref). Thus, for the time being,
cycles which are *not* disjoint actually are supported.

# Examples
```jldoctest
julia> x = @perm (1,2,3)(4,5)(factorial(3),7,8)
(1,2,3)(4,5)(6,7,8)

julia> parent(x)
Sym(8)

julia> y = cperm([1,2,3],[4,5],[6,7,8])
(1,2,3)(4,5)(6,7,8)

julia> x == y
true

julia> z = perm(symmetric_group(8),[2,3,1,5,4,7,8,6])
(1,2,3)(4,5)(6,7,8)

julia> x == z
true
```
"""
macro perm(ex)
    res = _perm_helper(ex)
    return esc(:(Oscar.cperm($(res...))))
end


################################################################################
#
#   perm(n,gens)
#
@doc raw"""
    @perm n gens
    
Input a list of permutations in cycle notation, created as elements of the
symmetric group of degree `n`, i.e., `symmetric_group(n)`, by invoking
[`cperm`](@ref) suitably.

# Examples
```jldoctest
julia> gens = @perm 14 [
              (1,10)
              (2,11)
              (3,12)
              (4,13)
              (5,14)
              (6,8)
              (7,9)
              (1,2,3,4,5,6,7)(8,9,10,11,12,13,14)
              (1,2)(10,11)
             ]
9-element Vector{PermGroupElem}:
 (1,10)
 (2,11)
 (3,12)
 (4,13)
 (5,14)
 (6,8)
 (7,9)
 (1,2,3,4,5,6,7)(8,9,10,11,12,13,14)
 (1,2)(10,11)
 
julia> parent(gens[1])
Sym(14)
```
"""
macro perm(n,gens)

    ores = Expr[]
    for ex in gens.args
        res = _perm_helper(ex)
        push!(ores, esc(:(  [$(res...)]  )))
    end

    return quote
       let g = symmetric_group($n)
           [ cperm(g, pi...) for pi in [$(ores...)] ]
       end
    end
end


@doc raw"""
    permutation_group(n::IntegerUnion, perms::Vector{PermGroupElem})

Return the permutation group of degree `n` that is generated by the
elements in `perms`.

# Examples
```jldoctest
julia> x = cperm([1,2,3], [4,5]);  y = cperm([1,4]);

julia> permutation_group(5, [x, y])
Permutation group of degree 5
```
"""
function permutation_group(n::IntegerUnion, perms::Vector{PermGroupElem})
  return sub(symmetric_group(n), perms)[1]
end

@doc raw"""
    @permutation_group(n, gens...)

Input the permutation group of degree `n` with generators `gens...`,
given by permutations in cycle notation.

# Examples
```jldoctest
julia> g = @permutation_group(7, (1,2), (1,2,3)(4,5))
Permutation group of degree 7

julia> degree(g)
7
```
"""
macro permutation_group(n, gens...)
    ores = Expr[]
    for ex in gens
        res = _perm_helper(ex)
        push!(ores, esc(:([$(res...)])))
    end

    return quote
       let g = symmetric_group($n)
           sub(g, [cperm(g, pi...) for pi in [$(ores...)]], check = false)[1]
       end
    end
end
