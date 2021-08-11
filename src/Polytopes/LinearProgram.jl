@doc Markdown.doc"""
    LinearProgram(P, c; k = 0, convention = :max)

The linear program on the feasible set `P` (a Polyhedron) with
respect to the function x ↦ dot(c,x)+k.

"""
struct LinearProgram
   feasible_region::Polyhedron
   polymake_lp::Polymake.BigObject
   convention::Symbol
end

function LinearProgram(Q::Polyhedron, objective::AbstractVector; k = 0, convention = :max)
   if convention != :max && convention != :min
      throw(ArgumentError("convention must be set to :min or :max."))
   end
   P=Polyhedron(Polymake.polytope.Polytope(pm_polytope(Q)))
   ambDim = ambient_dim(P)
   size(objective, 1) == ambDim || error("objective has wrong dimension.")
   lp = Polymake.polytope.LinearProgram(LINEAR_OBJECTIVE=homogenize(objective, k))
   if convention == :max
      Polymake.attach(lp, "convention", "max")
   elseif convention == :min
      Polymake.attach(lp, "convention", "min")
   end
   pm_polytope(P).LP = lp
   LinearProgram(P, lp, convention)
end

LinearProgram(A::Union{Oscar.MatElem,AbstractMatrix}, b, c::AbstractVector; k = 0, convention = :max) =
   LinearProgram(Polyhedron(A, b), c;  k = k, convention = convention)


###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, LP::LinearProgram)
    c=dehomogenize(LP.polymake_lp.LINEAR_OBJECTIVE)
    k=LP.polymake_lp.LINEAR_OBJECTIVE[1]
    print(io, "The linear program\n")
    if LP.convention == :max
      print(io, "   max")
    elseif LP.convention == :min
      print(io, "   min")
    end
    print(io, "{c⋅x + k | x ∈ P}\n")
    print(io, "where P is a "*string(typeof(LP.feasible_region)))
    print(io, " and\n   c=")
    print(io, string(c'))
    print(io, "\n   k=")
    print(io, string(k))
end


###############################################################################
###############################################################################
### Access
###############################################################################
###############################################################################

"""
    objective_function(LP::LinearProgram; as = :pair)

Return the objective function x ↦ dot(c,x)+k of the linear program LP.
The allowed values for `as` are
* `pair`: Return the pair `(c,k)`
* `function`: Return the objective function as a function.


"""
function objective_function(lp::LinearProgram; as::Symbol = :pair)
   if as == :pair
      return dehomogenize(lp.polymake_lp.LINEAR_OBJECTIVE),lp.polymake_lp.LINEAR_OBJECTIVE[1]
   elseif as == :function
      (c,k) = objective_function(lp, as = :pair)
      return x -> sum(x.*c)+k
   else
       throw(ArgumentError("Unsupported `as` argument :" * string(as)))
   end
end

"""
    feasible_region(lp::LinearProgram)

Return the feasible region of the linear program `lp`, which is a `Polyhedron`.
"""
feasible_region(lp::LinearProgram) = lp.feasible_region


###############################################################################
###############################################################################
### Solving of linear programs
###############################################################################
###############################################################################

"""
    minimal_vertex(LP::LinearProgram)

Return either a point of the feasible region of `LP` which minimizes the objective
function of `LP`, or `nothing` if no such point exists.

# Examples
The following example constructs a linear program over the three dimensional cube.
Although the linear program is given using the `:max` convention, one may still call
`minimal_vertex`.
```jldoctest
julia> C=cube(3)
A polyhedron in ambient dimension 3

julia> LP=LinearProgram(C,[1,2,-3])
The linear program
   max{c⋅x + k | x ∈ P}
where P is a Polyhedron and
   c=Polymake.Rational[1 2 -3]
   k=0

julia> minimal_vertex(LP)
pm::Vector<pm::Rational>
-1 -1 1
```
"""
function minimal_vertex(lp::LinearProgram)
   mv = lp.polymake_lp.MINIMAL_VERTEX
   if mv != nothing
      return dehomogenize(mv)
   else
      return nothing
   end
end

"""
    maximal_vertex(LP::LinearProgram)

    Return, if it exists, a point of the feasible region of `LP` which maximizes the objective
    function of `LP`.

    Return `nothing` if the objective function does not obtain a maximal value over the feasible region.

# Examples
The following example constructs a linear program over the three dimensional cube, and
obtains the vertex of the cube which maximizes the function (x,y,z) ↦ x+2y-3z.
```jldoctest
julia> C=cube(3)
A polyhedron in ambient dimension 3

julia> LP=LinearProgram(C,[1,2,-3])
The linear program
   max{c⋅x + k | x ∈ P}
where P is a Polyhedron and
   c=Polymake.Rational[1 2 -3]
   k=0

julia> maximal_vertex(LP)
pm::Vector<pm::Rational>
1 1 -1
```
"""
function maximal_vertex(lp::LinearProgram)
   mv = lp.polymake_lp.MAXIMAL_VERTEX
   if mv != nothing
      return dehomogenize(mv)
   else
      return nothing
   end
end


"""
    minimal_value(LP::LinearProgram)

Return, if it exists, the minimal value of the objective function of `LP` over the feasible region
of `LP`. Otherwise, return `-inf`.

# Examples
The following example constructs a linear program over the three dimensional cube, and
obtains the minimal value of the function (x,y,z) ↦ x+2y-3z over that cube.
```jldoctest
julia> C=cube(3)
A polyhedron in ambient dimension 3

julia> LP=LinearProgram(C,[1,2,-3])
The linear program
   max{c⋅x + k | x ∈ P}
where P is a Polyhedron and
   c=Polymake.Rational[1 2 -3]
   k=0

julia> minimal_value(LP)
-6
```
"""
minimal_value(lp::LinearProgram) = lp.polymake_lp.MINIMAL_VALUE


"""
    maximal_value(LP::LinearProgram)

Return the maximal value of the objective function of `LP` over the feasible region
of `LP`. Otherwise, return `inf`

# Examples
The following example constructs a linear program over the three dimensional cube, and
obtains the maximal value of the function (x,y,z) ↦ x+2y-3z over that cube.
```jldoctest
julia> C=cube(3)
A polyhedron in ambient dimension 3

julia> LP=LinearProgram(C,[1,2,-3])
The linear program
   max{c⋅x + k | x ∈ P}
where P is a Polyhedron and
   c=Polymake.Rational[1 2 -3]
   k=0

julia> maximal_value(LP)
6
```
"""
maximal_value(lp::LinearProgram) = lp.polymake_lp.MAXIMAL_VALUE

"""
    solve_lp(LP::LinearProgram)

Return a pair `(m,v)` where the optimal value `m` of the objective
 function of `LP` is attained at `v` (if `m` exists). If the optimum
 is not attained, `m` may be `inf` or `-inf` in which case `v` is
 `nothing`.
"""
function solve_lp(lp::LinearProgram)
   if lp.convention == :max
      return maximal_value(lp),maximal_vertex(lp)
   elseif lp.convention == :min
      return minimal_value(lp),minimal_vertex(lp)
   end
end
