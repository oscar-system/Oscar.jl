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

@doc Markdown.doc"""
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
Sym( [ 1 .. 6 ] )
```
"""
function perm(L::AbstractVector{<:IntegerUnion})
  return PermGroupElem(symmetric_group(length(L)), GAP.Globals.PermList(GAP.GapObj(L;recursive=true)))
end


@doc Markdown.doc"""
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

# cperm stands for "cycle permutation", but we can change name if we want
# takes as input a list of vectors (not necessarily disjoint)
@doc Markdown.doc"""
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

julia> p = cperm([1,2,3],[7])
(1,2,3)

julia> degree(parent(p))
7
```

Two permutations coincide if, and only if, they move the same points and their parent groups have the same degree.
```jldoctest
julia> G=symmetric_group(5);

julia> A=alternating_group(5);

julia> x=cperm(G,[1,2,3]);

julia> y=cperm(A,[1,2,3]);

julia> z=cperm([1,2,3]); parent(z)
Sym( [ 1 .. 3 ] )

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
Sym( [ 1 .. 5 ] )

julia> x = cperm(G,[[1,2],[3,4]])
(1,2)(3,4)

julia> parent(x)
Sym( [ 1 .. 5 ] )
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
function cperm(L::AbstractVector{T}...) where T <: IntegerUnion
   if length(L)==0
      return one(symmetric_group(1))
   else
      return prod([PermGroupElem(symmetric_group(maximum(y)), GAP.Globals.CycleFromList(GAP.Obj([Int(k) for k in y]))) for y in L])
#TODO: better create the product of GAP permutations?
   end
end

# cperm stays for "cycle permutation", but we can change name if we want
# takes as input a list of vectors (not necessarily disjoint)
# WARNING: we allow e.g. PermList([2,3,1,4,5,6]) in Sym(3)
function cperm(g::PermGroup,L::AbstractVector{T}...) where T <: IntegerUnion
   if length(L)==0
      return one(g)
   else
      x=prod(y -> GAP.Globals.CycleFromList(GAP.Obj([Int(k) for k in y])), L)
      if length(L) <= degree(g) && x in g.X
         return PermGroupElem(g, x)
      else
         throw(ArgumentError("the element does not embed in the group"))
      end
   end
end

function cperm(L::Vector{Vector{T}}) where T <: IntegerUnion
    return cperm(L...)
end

function cperm(g::PermGroup,L::Vector{Vector{T}}) where T <: IntegerUnion
    return cperm(g,L...)
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

##
# cycle-types and support
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
  function CycleType()
    return new(Pair{Int, Int}[])
  end
  function CycleType(v::Vector{Pair{Int, Int}}; sorted::Bool = false)
    sorted && return new(v)
    return new(sort(v, by = x -> x[1]))
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

function order(c::CycleType)
  return reduce(lcm, map(x->fmpz(x[1]), c.s), init = fmpz(1))
end

"""
    cycle_structure(g::PermGroupElem)

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

julia> cperm()
()

julia> cycle_structure(ans)
1-element Oscar.CycleType:
 1 => 1
```
"""
function cycle_structure(g::PermGroupElem)
    c = GAP.Globals.CycleStructurePerm(GAP.GapObj(g))
    # TODO: use SortedDict from DataStructures.jl ?
    ct = Pair{Int, Int}[ i+1 => c[i] for i in 1:length(c) if GAP.Globals.ISB_LIST(c, i) ]
    s = mapreduce(x->x[1]*x[2], +, ct, init = Int(0))
    if s < degree(parent(g))
      @assert length(c) == 0 || ct[1][1] > 1
      insert!(ct, 1, 1=>degree(parent(g))-s)
    end
    return CycleType(ct, sorted = true)
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
@doc Markdown.doc"""
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
Sym( [ 1 .. 8 ] )

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
@doc Markdown.doc"""
    @perm n gens
    
Input a list of permutations in cycle notation, created as elements of the
symmetric group of degree `n`, i.e., `symmetric_group(n)`, by invoking
[`cperm`](@ref) suitably..

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

export @perm
