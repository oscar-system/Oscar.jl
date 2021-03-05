@doc Markdown.doc"""
    LinearProgram(P, c [, k])

The linear program on the feasible set P (a Polyhedron) with
 respect to the function x ↦ dot(c,x)+k where k is optional (default 0).

"""
struct LinearProgram
   feasible_region::Polyhedron
   polymake_lp::Polymake.BigObject
   convention::Symbol
   function LinearProgram(Q::Polyhedron, objective::AbstractVector, k = 0; convention = :max)
      P=Polyhedron(Polymake.polytope.Polytope(pm_polytope(Q)))
      ambDim = ambient_dim(P)
      size(objective, 1) == ambDim || error("objective has wrong dimension.")
      lp = Polymake.polytope.LinearProgram(LINEAR_OBJECTIVE=homogenize(objective, k))
      pm_polytope(P).LP = lp
      new(P, lp, convention)
   end
end


###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, lp::LinearProgram)
    c=dehomogenize(lp.polymake_lp.LINEAR_OBJECTIVE)
    k=lp.polymake_lp.LINEAR_OBJECTIVE[1]
    print(io, "The linear program\n")
    if lp.convention == :max
      print(io, "   max")
    elseif lp.convention == :min
      print(io, "   min")
    end
    print(io, "{dot(c,x)+k | x ∈ P}\n")
    print(io, "where P is a "*string(typeof(lp.feasible_region)))
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

Returns the objective function x ↦ dot(c,x)+k of the linear program LP.
The allowed values for `as` are
* `pair`: Returns the pair `(c,k)`
* `function`: Returns the objective function as a function.


"""
function objective_function(LP::LinearProgram; as::Symbol = :pair)
   if as == :pair
      return dehomogenize(LP.polymake_lp.LINEAR_OBJECTIVE),LP.polymake_lp.LINEAR_OBJECTIVE[1]
   elseif as == :function
      (c,k) = objective_function(LP, as = :pair)
      return x -> sum(x.*c)+k
   else
       throw(ArgumentError("Unsupported `as` argument :" * string(as)))
   end
end

feasible_region(lp::LinearProgram) = lp.feasible_region


###############################################################################
###############################################################################
### Solving of linear programs
###############################################################################
###############################################################################


function minimal_vertex(lp::LinearProgram)
   mv = lp.polymake_lp.MINIMAL_VERTEX
   if mv != nothing
      return dehomogenize(mv)
   else
      return nothing
   end
end

function maximal_vertex(lp::LinearProgram)
   mv = lp.polymake_lp.MAXIMAL_VERTEX
   if mv != nothing
      return dehomogenize(mv)
   else
      return nothing
   end
end

minimal_value(lp::LinearProgram) = lp.polymake_lp.MINIMAL_VALUE
maximal_value(lp::LinearProgram) = lp.polymake_lp.MAXIMAL_VALUE

"""
   solve_lp(lp::LinearProgram)

Gives a pair `(m,v)` where the optimal value `m` of the objective
 function of lp is attained at `v` (if `m` exists). If the optimum
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


###############################################################################
###############################################################################
### Construction of linear programs
###############################################################################
###############################################################################


"""
   primal_program(c, A, b)

Constructs the primal linear program max{dot(c,x) | Ax<= b}.

see Def. 4.10
"""
primal_program(c, A, b) = LinearProgram(Polyhedron(A,b), c; convention = :max)

# """
#    DualProgram(c, A, b)
#
# Constructs the dual linear program min{yb | yA=c, y>=0}.
#
# see Def. 4.10
# """
# DualProgram(c, A, b) = LinearProgram(DualPolyhedron(A,c), b)
