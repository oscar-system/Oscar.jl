struct MixedIntegerLinearProgram{T} <: PolyhedralObject{T}
  feasible_region::Polyhedron{T}
  polymake_milp::Polymake.BigObject
  convention::Symbol
  parent_field::Field

  MixedIntegerLinearProgram{T}(
    fr::Polyhedron{T}, milp::Polymake.BigObject, c::Symbol, parent_field::Field
  ) where {T<:scalar_types} = new{T}(fr, milp, c, parent_field)
end

# no default = `QQFieldElem` here; scalar type can be derived from the feasible region
mixed_integer_linear_program(p::Polyhedron{T}, x...) where {T<:scalar_types} =
  MixedIntegerLinearProgram{T}(p, x..., coefficient_field(p))

@doc raw"""
    mixed_integer_linear_program(P, c; integer_variables = [], k = 0, convention = :max)

The mixed integer linear program on the feasible set `P` (a Polyhedron) with
respect to the function x ↦ dot(c,x)+k, where $x_i\in\mathbb{Z}$ for all $i$ in
`integer_variables`. If `integer_variables` is empty, or not specified, all
entries of `x` are required to be integral. If this is not intended, consider
using a [`linear_program`](@ref) instead.
"""
function mixed_integer_linear_program(
  P::Polyhedron{T},
  objective::AbstractVector;
  integer_variables=Int64[],
  k=0,
  convention=:max,
) where {T<:scalar_types}
  if convention != :max && convention != :min
    throw(ArgumentError("convention must be set to :min or :max."))
  end
  ambDim = ambient_dim(P)
  if isnothing(integer_variables) || isempty(integer_variables)
    integer_variables = 1:ambDim
  end
  size(objective, 1) == ambDim || error("objective has wrong dimension.")
  cf = coefficient_field(P)
  objective = cf.(objective)
  milp = Polymake.polytope.MixedIntegerLinearProgram{_scalar_type_to_polymake(T)}(;
    LINEAR_OBJECTIVE=homogenize(cf, objective, k), INTEGER_VARIABLES=Vector(integer_variables)
  )
  if convention == :max
    Polymake.attach(milp, "convention", "max")
  elseif convention == :min
    Polymake.attach(milp, "convention", "min")
  end
  Polymake.add(pm_object(P), "MILP", milp)
  MixedIntegerLinearProgram{T}(P, milp, convention, cf)
end

function mixed_integer_linear_program(
  ::Type{T},
  A::Union{Oscar.MatElem,AbstractMatrix},
  b,
  c::AbstractVector;
  integer_variables=Vector{Int64}([]),
  k=0,
  convention=:max,
) where {T<:scalar_types}
  P = polyhedron(T, A, b)
  return mixed_integer_linear_program(
    P,
    c,
    coefficient_field(P);
    integer_variables=integer_variables,
    k=k,
    convention=convention,
  )
end

pm_object(milp::MixedIntegerLinearProgram) = milp.polymake_milp

###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function describe(io::IO, MILP::MixedIntegerLinearProgram)
  c = MILP.polymake_milp.LINEAR_OBJECTIVE[2:end]
  k = MILP.polymake_milp.LINEAR_OBJECTIVE[1]
  print(io, "The mixed integer linear program\n")
  if MILP.convention == :max
    print(io, "   max")
  elseif MILP.convention == :min
    print(io, "   min")
  end
  if is_unicode_allowed()
    print(io, "{c⋅x + k | x ∈ P}\n")
  else
    print(io, "{c*x + k | x in P}\n")
  end
  print(io, "where P is a " * string(typeof(MILP.feasible_region)))
  print(io, "\n   c=")
  print(io, string(c'))
  print(io, "\n   k=")
  print(io, string(k))
  print(io, "\n   ")
  ivar = _integer_variables(MILP)
  if length(ivar) == 0
    print(io, "all entries of x in ZZ")
  else
    print(io, join(["x" * string(i) for i in ivar], ",") * " in ZZ")
  end
end

Base.show(io::IO, MILP::MixedIntegerLinearProgram) =
  print(io, "Mixed integer linear program")

###############################################################################
###############################################################################
### Access
###############################################################################
###############################################################################
_integer_variables(milp::MixedIntegerLinearProgram) = milp.polymake_milp.INTEGER_VARIABLES

@doc raw"""
    objective_function(MILP::MixedIntegerLinearProgram; as = :pair)

Return the objective function x ↦ dot(c,x)+k of the mixed integer linear program MILP.
The allowed values for `as` are
* `pair`: Return the pair `(c,k)`
* `function`: Return the objective function as a function.


"""
function objective_function(
  milp::MixedIntegerLinearProgram{T}; as::Symbol=:pair
) where {T<:scalar_types}
  if as == :pair
    cf = coefficient_field(milp)
    return T[cf(x) for x in milp.polymake_milp.LINEAR_OBJECTIVE[2:end]],
    cf.(milp.polymake_milp.LINEAR_OBJECTIVE[1])
  elseif as == :function
    (c, k) = objective_function(milp; as=:pair)
    return x -> sum(x .* c) + k
  else
    throw(ArgumentError("Unsupported `as` argument: $as"))
  end
end

@doc raw"""
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

@doc raw"""
    optimal_solution(MILP::MixedIntegerLinearProgram)

Return either a point of the feasible region of `MILP` which optimizes the
objective function of `MILP`, or `nothing` if no such point exists.

# Examples
Take the square $[-1/2,3/2]^2$ and optimize $[1,1]$ in different settings.
```jldoctest
julia> c = cube(2, -1//2, 3//2)
Polytope in ambient dimension 2

julia> milp = mixed_integer_linear_program(c, [1,1], integer_variables=[1])
Mixed integer linear program

julia> optimal_solution(milp)
2-element PointVector{QQFieldElem}:
 1
 3//2

julia> milp = mixed_integer_linear_program(c, [1,1])
Mixed integer linear program

julia> optimal_solution(milp)
2-element PointVector{QQFieldElem}:
 1
 1
```
"""
function optimal_solution(milp::MixedIntegerLinearProgram{T}) where {T<:scalar_types}
  opt_vert = nothing
  if milp.convention == :max
    opt_vert = milp.polymake_milp.MAXIMAL_SOLUTION
  else
    opt_vert = milp.polymake_milp.MINIMAL_SOLUTION
  end
  if !isnothing(opt_vert)
    return point_vector(coefficient_field(milp), dehomogenize(opt_vert))::PointVector{T}
  else
    return nothing
  end
end

@doc raw"""
    optimal_value(MILP::MixedIntegerLinearProgram)

Return, if it exists, the optimal value of the objective function of `MILP`
over the feasible region of `MILP`. Otherwise, return `-inf` or `inf` depending
on convention.

# Examples
Take the square $[-1/2,3/2]^2$ and optimize $[1,1]$ in different settings.
```jldoctest
julia> c = cube(2, -1//2, 3//2)
Polytope in ambient dimension 2

julia> milp = mixed_integer_linear_program(c, [1,1], integer_variables=[1])
Mixed integer linear program

julia> optimal_value(milp)
5/2

julia> milp = mixed_integer_linear_program(c, [1,1])
Mixed integer linear program

julia> optimal_value(milp)
2
```
"""
function optimal_value(milp::MixedIntegerLinearProgram{T}) where {T<:scalar_types}
  if milp.convention == :max
    # TODO: consider inf
    return milp.polymake_milp.MAXIMAL_VALUE
  else
    return milp.polymake_milp.MINIMAL_VALUE
  end
end

@doc raw"""
    solve_milp(MILP::MixedIntegerLinearProgram)

Return a pair `(m,v)` where the optimal value `m` of the objective function of
`MILP` is attained at `v` (if `m` exists). If the optimum is not attained, `m`
may be `inf` or `-inf` in which case `v` is `nothing`.

# Examples
Take the square $[-1/2,3/2]^2$ and optimize $[1,1]$ in different settings.
```jldoctest
julia> c = cube(2, -1//2, 3//2)
Polytope in ambient dimension 2

julia> milp = mixed_integer_linear_program(c, [1,1], integer_variables=[1])
Mixed integer linear program

julia> solve_milp(milp)
(5/2, QQFieldElem[1, 3//2])

julia> milp = mixed_integer_linear_program(c, [1,1])
Mixed integer linear program

julia> solve_milp(milp)
(2, QQFieldElem[1, 1])
```
"""
solve_milp(milp::MixedIntegerLinearProgram) = optimal_value(milp), optimal_solution(milp)

@doc raw"""
    ambient_dim(MILP::MixedIntegerLinearProgram)

Return the ambient dimension of the feasible reagion of `MILP`.
"""
ambient_dim(milp::MixedIntegerLinearProgram) = ambient_dim(feasible_region(milp))
