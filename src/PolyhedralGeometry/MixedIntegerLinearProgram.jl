@doc Markdown.doc"""
    MixedIntegerLinearProgram(P, c; integer_variables = [], k = 0, convention = :max)

The mixed integer linear program on the feasible set `P` (a Polyhedron) with
respect to the function x ↦ dot(c,x)+k, where $x_i\in\mathbb{Z}$ for all $i$ in
`integer_variables`. If `integer_variables` is empty, or not specified, all
entries of `x` are required to be integral. If this is not intended, consider
using a [`LinearProgram`](@ref) instead.
"""
struct MixedIntegerLinearProgram{T}
   feasible_region::Polyhedron{T}
   polymake_milp::Polymake.BigObject
   convention::Symbol
   
   MixedIntegerLinearProgram{T}(fr::Polyhedron{T}, milp::Polymake.BigObject, c::Symbol) where T<:scalar_types = new{T}(fr, milp, c)
end

# no default = `fmpq` here; scalar type can be derived from the feasible region
MixedIntegerLinearProgram(p::Polyhedron{T}, x...) where T<:scalar_types = MixedIntegerLinearProgram{T}(p, x...)

function MixedIntegerLinearProgram{T}(P::Polyhedron{T}, objective::AbstractVector; integer_variables=Int64[], k = 0, convention = :max) where T<:scalar_types
   if convention != :max && convention != :min
      throw(ArgumentError("convention must be set to :min or :max."))
   end
   ambDim = ambient_dim(P)
   if isnothing(integer_variables) || isempty(integer_variables)
      integer_variables = 1:ambDim
   end
   size(objective, 1) == ambDim || error("objective has wrong dimension.")
   milp = Polymake.polytope.MixedIntegerLinearProgram{scalar_type_to_polymake[T]}(LINEAR_OBJECTIVE=homogenize(objective, k), INTEGER_VARIABLES=Vector(integer_variables))
   if convention == :max
      Polymake.attach(milp, "convention", "max")
   elseif convention == :min
      Polymake.attach(milp, "convention", "min")
   end
   Polymake.add(pm_object(P), "MILP", milp)
   MixedIntegerLinearProgram{T}(P, milp, convention)
end

MixedIntegerLinearProgram(Q::Polyhedron{T}, objective::AbstractVector; integer_variables=Vector{Int64}([]), k = 0, convention = :max) where T<:scalar_types = MixedIntegerLinearProgram{T}(Q, objective; integer_variables = integer_variables, k = k, convention = convention)

MixedIntegerLinearProgram{T}(A::Union{Oscar.MatElem,AbstractMatrix}, b, c::AbstractVector; integer_variables=Vector{Int64}([]), k = 0, convention = :max)  where T<:scalar_types =
   MixedIntegerLinearProgram{T}(Polyhedron{T}(A, b), c; integer_variables = integer_variables, k = k, convention = convention)

MixedIntegerLinearProgram(x...) = MixedIntegerLinearProgram{fmpq}(x...)

pm_object(milp::MixedIntegerLinearProgram) = milp.polymake_milp

###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function describe(io::IO, MILP::MixedIntegerLinearProgram)
    c=dehomogenize(MILP.polymake_milp.LINEAR_OBJECTIVE)
    k=MILP.polymake_milp.LINEAR_OBJECTIVE[1]
    print(io, "The mixed integer linear program\n")
    if MILP.convention == :max
      print(io, "   max")
    elseif MILP.convention == :min
      print(io, "   min")
    end
    print(io, "{c⋅x + k | x ∈ P}\n")
    print(io, "where P is a "*string(typeof(MILP.feasible_region)))
    print(io, "\n   c=")
    print(io, string(c'))
    print(io, "\n   k=")
    print(io, string(k))
    print(io, "\n   ")
    ivar = _integer_variables(MILP)
    if length(ivar) == 0
        print(io, "all entries of x in ZZ")
    else
        print(io, join(["x"*string(i) for i in ivar],",")*" in ZZ")
    end
end

function Base.show(io::IO, MILP::MixedIntegerLinearProgram)
    print(io, "A mixed integer linear program")
end


###############################################################################
###############################################################################
### Access
###############################################################################
###############################################################################
function _integer_variables(milp::MixedIntegerLinearProgram)
    return milp.polymake_milp.INTEGER_VARIABLES
end

@doc Markdown.doc"""
    objective_function(MILP::MixedIntegerLinearProgram; as = :pair)

Return the objective function x ↦ dot(c,x)+k of the mixed integer linear program MILP.
The allowed values for `as` are
* `pair`: Return the pair `(c,k)`
* `function`: Return the objective function as a function.


"""
function objective_function(milp::MixedIntegerLinearProgram{T}; as::Symbol = :pair) where T<:scalar_types
   if as == :pair
      return Vector{T}(dehomogenize(milp.polymake_milp.LINEAR_OBJECTIVE)),convert(T, milp.polymake_milp.LINEAR_OBJECTIVE[1])
   elseif as == :function
      (c,k) = objective_function(milp, as = :pair)
      return x -> sum(x.*c)+k
   else
       throw(ArgumentError("Unsupported `as` argument: $as"))
   end
end

@doc Markdown.doc"""
    feasible_region(milp::MixedIntegerLinearProgram)

Return the feasible region of the mixed integer linear program `milp`, which is
a `Polyhedron`.
"""
feasible_region(milp::MixedIntegerLinearProgram) = milp.feasible_region


###############################################################################
###############################################################################
### Solving of mixed integer linear programs
###############################################################################
###############################################################################


@doc Markdown.doc"""
    optimal_solution(MILP::MixedIntegerLinearProgram)

Return either a point of the feasible region of `MILP` which optimizes the
objective function of `MILP`, or `nothing` if no such point exists.

# Examples
Take the square $[-1/2,3/2]^2$ and optimize $[1,1]$ in different settings.
```jldoctest
julia> c = cube(2, -1//2, 3//2)
A polyhedron in ambient dimension 2

julia> milp = MixedIntegerLinearProgram(c, [1,1], integer_variables=[1])
A mixed integer linear program

julia> optimal_solution(milp)
2-element PointVector{fmpq}:
 1
 3//2

julia> milp = MixedIntegerLinearProgram(c, [1,1])
A mixed integer linear program

julia> optimal_solution(milp)
2-element PointVector{fmpq}:
 1
 1
```
"""
function optimal_solution(milp::MixedIntegerLinearProgram{T}) where T<:scalar_types
   opt_vert = nothing
   if milp.convention == :max
      opt_vert = milp.polymake_milp.MAXIMAL_SOLUTION
   else
      opt_vert = milp.polymake_milp.MINIMAL_SOLUTION
   end
   if opt_vert != nothing
      return PointVector{T}(dehomogenize(opt_vert))
   else
      return nothing
   end
end


@doc Markdown.doc"""
    optimal_value(MILP::MixedIntegerLinearProgram)

Return, if it exists, the optimal value of the objective function of `MILP`
over the feasible region of `MILP`. Otherwise, return `-inf` or `inf` depending
on convention.

# Examples
Take the square $[-1/2,3/2]^2$ and optimize $[1,1]$ in different settings.
```jldoctest
julia> c = cube(2, -1//2, 3//2)
A polyhedron in ambient dimension 2

julia> milp = MixedIntegerLinearProgram(c, [1,1], integer_variables=[1])
A mixed integer linear program

julia> optimal_value(milp)
5/2

julia> milp = MixedIntegerLinearProgram(c, [1,1])
A mixed integer linear program

julia> optimal_value(milp)
2
```
"""
function optimal_value(milp::MixedIntegerLinearProgram{T}) where T<:scalar_types
   if milp.convention == :max
      # TODO: consider inf
      return milp.polymake_milp.MAXIMAL_VALUE
   else
      return milp.polymake_milp.MINIMAL_VALUE
   end
end


@doc Markdown.doc"""
    solve_milp(MILP::MixedIntegerLinearProgram)

Return a pair `(m,v)` where the optimal value `m` of the objective function of
`MILP` is attained at `v` (if `m` exists). If the optimum is not attained, `m`
may be `inf` or `-inf` in which case `v` is `nothing`.

# Examples
Take the square $[-1/2,3/2]^2$ and optimize $[1,1]$ in different settings.
```jldoctest
julia> c = cube(2, -1//2, 3//2)
A polyhedron in ambient dimension 2

julia> milp = MixedIntegerLinearProgram(c, [1,1], integer_variables=[1])
A mixed integer linear program

julia> solve_milp(milp)
(5/2, fmpq[1, 3//2])

julia> milp = MixedIntegerLinearProgram(c, [1,1])
A mixed integer linear program

julia> solve_milp(milp)
(2, fmpq[1, 1])
```
"""
function solve_milp(milp::MixedIntegerLinearProgram)
   return optimal_value(milp),optimal_solution(milp)
end
