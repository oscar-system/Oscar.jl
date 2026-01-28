struct LinearProgram{T} <: PolyhedralObject{T}
  feasible_region::Polyhedron{T}
  polymake_lp::Polymake.BigObject
  convention::Symbol
  parent_field::Field

  LinearProgram{T}(
    fr::Polyhedron{T}, lp::Polymake.BigObject, c::Symbol, p
  ) where {T<:scalar_types} = new{T}(fr, lp, c, p)
  LinearProgram{QQFieldElem}(
    fr::Polyhedron{QQFieldElem}, lp::Polymake.BigObject, c::Symbol
  ) = new{QQFieldElem}(fr, lp, c, QQ)
end

# no default = `QQFieldElem` here; scalar type can be derived from the feasible region
linear_program(p::Polyhedron{T}, lp, c) where {T<:scalar_types} =
  LinearProgram{T}(p, lp, c, coefficient_field(p))

@doc raw"""
    linear_program(P, c; k = 0, convention = :max)

The linear program on the feasible set `P` (a Polyhedron) with
respect to the function x ↦ dot(c,x)+k.

"""
function linear_program(
  P::Polyhedron{T}, objective::AbstractVector; k=0, convention=:max
) where {T<:scalar_types}
  if convention != :max && convention != :min
    throw(ArgumentError("convention must be set to :min or :max."))
  end
  ambDim = ambient_dim(P)
  cf = coefficient_field(P)
  objective = cf.(objective)
  size(objective, 1) == ambDim || error("objective has wrong dimension.")
  lp = Polymake.polytope.LinearProgram{_scalar_type_to_polymake(T)}(;
    LINEAR_OBJECTIVE=homogenize(cf, objective, k)
  )
  if convention == :max
    Polymake.attach(lp, "convention", "max")
  elseif convention == :min
    Polymake.attach(lp, "convention", "min")
  end
  Polymake.add(pm_object(P), "LP", lp)
  LinearProgram{T}(P, lp, convention, cf)
end

linear_program(
  f::scalar_type_or_field,
  A::AbstractCollection[AffineHalfspace],
  b,
  c::AbstractVector;
  k=0,
  convention=:max,
) = linear_program(polyhedron(f, A, b), c; k=k, convention=convention)

pm_object(lp::LinearProgram) = lp.polymake_lp

###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, LP::LinearProgram)
  (c, k) = objective_function(LP; as=:pair)
  print(io, "Linear program\n")
  if LP.convention == :max
    print(io, "   max")
  elseif LP.convention == :min
    print(io, "   min")
  end
  if is_unicode_allowed()
    print(io, "{c⋅x + k | x ∈ P}\n")
  else
    print(io, "{c*x + k | x in P}\n")
  end
  print(io, "where P is a " * string(typeof(LP.feasible_region)))
  print(io, " and\n   c=")
  print(io, string(c))
  print(io, "\n   k=")
  print(io, string(k))
end

###############################################################################
###############################################################################
### Access
###############################################################################
###############################################################################

@doc raw"""
    objective_function(LP::LinearProgram; as = :pair)

Return the objective function x ↦ dot(c,x)+k of the linear program LP.
The allowed values for `as` are
* `pair`: Return the pair `(c,k)`
* `function`: Return the objective function as a function.


"""
function objective_function(lp::LinearProgram{T}; as::Symbol=:pair) where {T<:scalar_types}
  if as == :pair
    cf = coefficient_field(lp)
    return T[cf(x) for x in lp.polymake_lp.LINEAR_OBJECTIVE[2:end]],
    cf.(lp.polymake_lp.LINEAR_OBJECTIVE[1])
  elseif as == :function
    (c, k) = objective_function(lp; as=:pair)
    return x -> sum(x .* c) + k
  else
    throw(ArgumentError("Unsupported `as` argument: $as"))
  end
end

@doc raw"""
    feasible_region(lp::LinearProgram)

Return the feasible region of the linear program `lp`, which is a `Polyhedron`.
"""
feasible_region(lp::LinearProgram) = lp.feasible_region

###############################################################################
###############################################################################
### Solving of linear programs
###############################################################################
###############################################################################

@doc raw"""
    optimal_vertex(LP::LinearProgram)

Return either a point of the feasible region of `LP` which optimizes the objective
function of `LP`, or `nothing` if no such point exists.

# Examples
The following example constructs a linear program over the three dimensional cube, and
obtains the vertex of the cube which maximizes the function (x,y,z) ↦ x+2y-3z.
```jldoctest
julia> C=cube(3)
Polytope in ambient dimension 3

julia> LP=linear_program(C,[1,2,-3])
Linear program
   max{c*x + k | x in P}
where P is a Polyhedron{QQFieldElem} and
   c=QQFieldElem[1, 2, -3]
   k=0

julia> optimal_vertex(LP)
3-element PointVector{QQFieldElem}:
 1
 1
 -1
```
"""
function optimal_vertex(lp::LinearProgram{T}) where {T<:scalar_types}
  opt_vert = nothing
  if lp.convention == :max
    opt_vert = lp.polymake_lp.MAXIMAL_VERTEX
  else
    opt_vert = lp.polymake_lp.MINIMAL_VERTEX
  end
  if !isnothing(opt_vert)
    return point_vector(
      coefficient_field(lp), view(dehomogenize(opt_vert), :)
    )::PointVector{T}
  else
    return nothing
  end
end

@doc raw"""
    optimal_value(LP::LinearProgram)

Return, if it exists, the optimal value of the objective function of `LP` over the feasible region
of `LP`. Otherwise, return `-infinity` or `infinity` depending on convention, or `nothing` if the feasible region is empty.

# Examples
The following example constructs a linear program over the three dimensional cube, and
obtains the minimal value of the function (x,y,z) ↦ x+2y-3z over that cube.
```jldoctest
julia> C=cube(3)
Polytope in ambient dimension 3

julia> LP=linear_program(C,[1,2,-3]; convention = :min)
Linear program
   min{c*x + k | x in P}
where P is a Polyhedron{QQFieldElem} and
   c=QQFieldElem[1, 2, -3]
   k=0

julia> optimal_value(LP)
-6
```

Optimizing in an unbounded direction yields infinity.
```jldoctest
julia> lp = linear_program(convex_hull([0 0], [1 0; 0 1]), [1, 1])
Linear program
   max{c*x + k | x in P}
where P is a Polyhedron{QQFieldElem} and
   c=QQFieldElem[1, 1]
   k=0

julia> optimal_value(lp)
infinity
```
"""
function optimal_value(lp::LinearProgram{T}) where {T<:scalar_types}
  cf = coefficient_field(lp)
  conv = lp.convention
  mv = conv === :max ? lp.polymake_lp.MAXIMAL_VALUE : lp.polymake_lp.MINIMAL_VALUE
  is_feasible(feasible_region(lp)) || return nothing
  isinf(mv) && return conv === :max ? PosInf() : NegInf()
  return cf(mv)
end

@doc raw"""
    solve_lp(LP::LinearProgram)

Return a pair `(m,v)` where the optimal value `m` of the objective
function of `LP` is attained at `v` (if `m` exists). If the optimum is not
attained or the feasible region is empty, `m` may be `infinity`, `-infinity`,
or `nothing` in which case `v` is `nothing`.
"""
solve_lp(lp::LinearProgram) = optimal_value(lp), optimal_vertex(lp)

@doc raw"""
    ambient_dim(LP::LinearProgram)

Return the ambient dimension of the feasible reagion of `LP`.
"""
ambient_dim(lp::LinearProgram) = ambient_dim(feasible_region(lp))
