#
#
#

Base.isfinite(G::PermGroup) = true

==(x::PermGroup, y::PermGroup) = x.deg == y.deg && x.X == y.X

==(x::PermGroupElem, y::PermGroupElem) = degree(parent(x)) == degree(parent(y)) && x.X == y.X

Base.:<(x::PermGroupElem, y::PermGroupElem) = x.X < y.X

Base.isless(x::PermGroupElem, y::PermGroupElem) = x<y


"""
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

"""
    moved_points(x::PermGroupElem)
    moved_points(G::PermGroup)

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

"""
    number_moved_points(::Type{T} = fmpz, x::PermGroupElem) where T <: IntegerUnion
    number_moved_points(::Type{T} = fmpz, G::PermGroup) where T <: IntegerUnion

Return the number of those points in `1:degree(x)` or `1:degree(G)`,
respectively, that are not mapped to themselves under the action `^`,
as an instance of `T`.

# Examples
```jldoctest
julia> g = symmetric_group(4);  s = sylow_subgroup(g, 3)[1];

julia> number_moved_points(s)
3

julia> number_moved_points(Int, s)
3

julia> number_moved_points(gen(s, 1))
3

```
"""
@gapattribute number_moved_points(x::Union{PermGroupElem,PermGroup}) = fmpz(GAP.Globals.NrMovedPoints(x.X))::fmpz

number_moved_points(::Type{T}, x::Union{PermGroupElem,PermGroup}) where T <: IntegerUnion = T(GAP.Globals.NrMovedPoints(x.X))::T

# FIXME: clashes with AbstractAlgebra.perm method
#function perm(L::AbstractVector{<:IntegerUnion})
#   return PermGroupElem(symmetric_group(length(L)), GAP.Globals.PermList(GAP.GapObj(L;recursive=true)))
#end
# FIXME: use name gap_perm for now
@doc Markdown.doc"""
    gap_perm(L::AbstractVector{<:IntegerUnion})

Return the permutation $x$ which maps every $i$ from `1` to $n$` = length(L)`
to `L`$[i]$.
The parent of $x$ is set to [`symmetric_group`](@ref)$(n)$.
An exception is thrown if `L` does not contain every integer from 1 to $n$
exactly once.

# Examples
```jldoctest
julia> gap_perm([2,4,6,1,3,5])
(1,2,4)(3,6,5)
```
"""
function gap_perm(L::AbstractVector{<:IntegerUnion})
  return PermGroupElem(symmetric_group(length(L)), GAP.Globals.PermList(GAP.GapObj(L;recursive=true)))
end

@doc Markdown.doc"""
    perm(G::PermGroup, L::AbstractVector{<:IntegerUnion})
    (G::PermGroup)(L::AbstractVector{<:IntegerUnion})

Return the permutation $x$ which maps every `i` from 1 to $n$` = length(L)`
to `L`$[i]$.
The parent of $x$ is `G`.
An exception is thrown if $x$ is not contained in `G`
or `L` does not contain every integer from 1 to $n$ exactly once.

For [`gap_perm`](@ref),
the parent group of $x$ is set to [`symmetric_group`](@ref)$(n)$.

# Examples
```jldoctest
julia> perm(symmetric_group(6),[2,4,6,1,3,5])
(1,2,4)(3,6,5)
```
"""
function perm(g::PermGroup, L::AbstractVector{<:IntegerUnion})
   x = GAP.Globals.PermList(GAP.GapObj(L;recursive=true))
   x == GAP.Globals.fail && throw(ArgumentError("the list does not describe a permutation"))
   if length(L) <= degree(g) && x in g.X
     return PermGroupElem(g, x)
   end
   throw(ArgumentError("the element does not embed in the group"))
end

perm(g::PermGroup, L::AbstractVector{<:fmpz}) = perm(g, [Int(y) for y in L])

function (g::PermGroup)(L::AbstractVector{<:IntegerUnion})
   x = GAP.Globals.PermList(GAP.GapObj(L;recursive=true))
   if length(L) <= degree(g) && x in g.X
     return PermGroupElem(g, x)
   end
   throw(ArgumentError("the element does not embed in the group"))
end

(g::PermGroup)(L::AbstractVector{<:fmpz}) = g([Int(y) for y in L])

# cperm stays for "cycle permutation", but we can change name if we want
# takes as input a list of vectors (not necessarly disjoint)
@doc Markdown.doc"""
    cperm(L::AbstractVector{<:T}...) where T <: IntegerUnion
    cperm(G::PermGroup, L::AbstractVector{<:T}...)

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

julia> p = cperm([1,2,3],[7])
(1,2,3)

julia> degree(parent(p))
7

```

At the moment, the input vectors of the function `cperm` need not be disjoint.

!!! warning
    If the function `perm` is evaluated in a vector of integers
    without specifying the group `G`,
    then the returned value is an element of the AbstractAlgebra.jl type
    `Perm{Int}`.
    For this reason, if one wants a permutation of type
    `GAPGroupElem{PermGroup}` without specifying a parent,
    one has to use the function `gap_perm`.
"""
function cperm(L::AbstractVector{T}...) where T <: IntegerUnion
   if length(L)==0
      return one(symmetric_group(1))
   else
      return prod([PermGroupElem(symmetric_group(maximum(y)), GAP.Globals.CycleFromList(GAP.julia_to_gap([Int(k) for k in y]))) for y in L])
#TODO: better create the product of GAP permutations?
   end
end

# cperm stays for "cycle permutation", but we can change name if we want
# takes as input a list of vectors (not necessarly disjoint)
# WARNING: we allow e.g. PermList([2,3,1,4,5,6]) in Sym(3)
function cperm(g::PermGroup,L::AbstractVector{T}...) where T <: IntegerUnion
   if length(L)==0
      return one(g)
   else
      x=prod(y -> GAP.Globals.CycleFromList(GAP.julia_to_gap([Int(k) for k in y])), L)
      if length(L) <= degree(g) && x in g.X
         return PermGroupElem(g, x)
      else
         throw(ArgumentError("the element does not embed in the group"))
      end
   end
end

"""
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
julia> Vector{fmpz}(pi, 2)
2-element Vector{fmpz}:
 2
 3

```
"""
Base.Vector{T}(x::PermGroupElem, n::Int = x.parent.deg) where T <: IntegerUnion = T[x(i) for i in 1:n]
Base.Vector(x::PermGroupElem, n::Int = x.parent.deg) = Vector{Int}(x,n)

#embedding of a permutation in permutation group
function (G::PermGroup)(x::PermGroupElem)
   x.X in G.X && return group_element(G, x.X)
   throw(ArgumentError("the element does not embed in the group"))
end

#evaluation function
(x::PermGroupElem)(n::IntegerUnion) = n^x

^(n::T, x::PermGroupElem) where T <: IntegerUnion = T(GAP.Obj(n)^x.X)

^(n::Int, x::PermGroupElem) = (n^x.X)::Int


"""
    sign(g::PermGroupElem)

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
Base.sign(g::PermGroupElem) = GAP.Globals.SignPerm(g.X)::Int

# TODO: document the following?
Base.sign(G::PermGroup) = GAP.Globals.SignPermGroup(G.X)


"""
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

"""
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

"""
    cycle_structure(g::PermGroupElem)

Return the cycle structure of the permutation `g` as a sorted vector of
pairs. A pair `k => n` in this vector indicates that the number
of ``k``-cycles in the cycle decomposition of `g` is equal to `n`.

# Examples
```jldoctest
julia> g = cperm(1:3, 4:5, 6:7, 8:10, 11:15)
(1,2,3)(4,5)(6,7)(8,9,10)(11,12,13,14,15)

julia> cycle_structure(g)
3-element Vector{Pair{Int64, Int64}}:
 2 => 2
 3 => 2
 5 => 1

julia> cperm()
()

julia> cycle_structure(ans)
Pair{Int64, Int64}[]
```
"""
function cycle_structure(g::PermGroupElem)
    c = GAP.Globals.CycleStructurePerm(GAP.GapObj(g))
    # TODO: use SortedDict from DataStructures.jl ?
    return Pair{Int,Int}[ i+1 => c[i] for i in 1:length(c) if GAP.Globals.ISB_LIST(c, i) ]
end


# The following code implements a new way to input permutations in Julia. For example
# it is possible to create a permutation as follow
# pi = Oscar.Permutations.@perm (1,2,3)(4,5)(6,7,8)
# > (1,2,3)(4,5)(6,7,8)
# For this we use macros to modify the syntax tree of (1,2,3)(4,5)(6,7,8) such that
# Julia can deal with the expression.


################################################################################
#
#   perm
#
@doc Markdown.doc"""
    @perm(ex)
    
Macro to input a permutation as
`pi = @perm (1,2,3)(4,5)(6,7,8)` to obtain
the permutation `(1,2,3)(4,5)(6,7,8)`, that is, the output of
`cperm([1,2,3],[4,5],[6,7,8])`.
# Examples
```jldoctest
julia> @perm (1,2,3)(4,5)(6,7,8)
(1,2,3)(4,5)(6,7,8)
```
"""
macro perm(ex)
    res = []
    
    if typeof(ex) != Expr
        error("Input is not a permutation.")
    end
    
    if ex.head != :call && ex.head != :tuple
        error("Input is not a permutation.")
    end
   
    while ex.head == :call
        pushfirst!(res, Expr(:vect, ex.args[2:end]...))
        ex = ex.args[1]
        if typeof(ex) != Expr
            error("Input is not a permutation.")
        end
        if ex.head != :call && ex.head != :tuple
            error("Input is not a permutation.")
        end
    end
    
    if ex.head != :tuple
        error("Input is not a permutation.")
    end
    
    pushfirst!(res, Expr(:vect,ex.args...))

    return esc(:(Oscar.cperm($(res...))))
end


################################################################################
#
#   perm(n,gens)
#
@doc Markdown.doc"""
    @perm(n,gens)
    
Macro to input a list of permutations which are generated as elements of
the `symmetric_group(n)` with the function `cperm`.
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
Sym( [ 1 .. 14 ] )
```
"""
macro perm(n,gens)

    ores = Vector{Expr}(undef,length(gens.args))
    i = 1
    for ex in gens.args
        res = []
        
        if typeof(ex) != Expr
            throw(ArgumentError("Input is not a permutation."))
            error("Input is not a permutation.")
        end
        
        if ex.head != :call && ex.head != :tuple
            error("Input is not a permutation.")
        end
       
        while ex.head == :call
            pushfirst!(res, Expr(:vect, ex.args[2:end]...))
            ex = ex.args[1]
            if typeof(ex) != Expr
                error("Input is not a permutation.")
            end
            if ex.head != :call && ex.head != :tuple
                error("Input is not a permutation.")
            end
        end
        
        if ex.head != :tuple
            error("Input is not a permutation.")
        end
        
        pushfirst!(res, Expr(:vect,ex.args...))

        ores[i] = esc(:(Oscar.cperm(symmetric_group($n),$(res...))))
        i = i + 1
    end

    return Expr(:vect,ores...)
end


################################################################################
#
#   permgroup(n::Int64,gens::Vector{PermGroupElem})
#
@doc Markdown.doc"""
    permgroup(n::Int64,gens::Vector{PermGroupElem})
    
Generates a `PermGroup` with generators gens as a subgroup of `symmetric_group(n)`.
# Examples
```jldoctest
julia> gens = Oscar.Permutations.@perm 14 [
              (1,10)
              (2,11)
              (3,12)
              (4,13)
              (5,14)
              (6,8)
              (7,9)
              (1,2,3,4,5,6,7)(8,9,10,11,12,13,14)
              (1,2)(10,11)
             ];
julia> G = Oscar.Permutations.permgroup(14,gens)
<permutation group with 9 generators>
```
"""
function permgroup(n::Int64,gens::Vector{PermGroupElem})

    return PermGroup(GAP.Globals.Subgroup(GAP.Globals.SymmetricGroup(GAP.Obj(n)),GAP.Obj([GAP.Obj(x) for x in gens ])))
end

export @perm,permgroup
