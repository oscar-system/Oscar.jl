function _internal_save_lp(io::IO, p::Polymake.BigObject, lp::Polymake.BigObject, max::Bool)
  # we cannot pass an IO object, so we use shell_execute to capture the stdout and write
  # the string to the IO object
  # to pass the object to the polymake shell we use a randomly generated variable name
  rstr = randstring(['A':'Z'; 'a':'z'], 12)
  Polymake.call_function(:User, :set_shell_scalar, rstr*"_poly", p)
  Polymake.call_function(:User, :set_shell_scalar, rstr*"_lp", lp)
  out, err = Polymake.shell_execute("""polytope::poly2lp(\$$(rstr)_poly, \$$(rstr)_lp, $(max ? 1 : 0));""")
  write(io, out)
  Polymake.call_function(:User, :set_shell_scalar, rstr*"_poly", 0)
  Polymake.call_function(:User, :set_shell_scalar, rstr*"_lp", 0)
  nothing
end

function _internal_save_mps(io::IO, p::Polymake.BigObject, lp::Polymake.BigObject)
  rstr = randstring(['A':'Z'; 'a':'z'], 12)
  Polymake.call_function(:User, :set_shell_scalar, rstr*"_poly", p)
  Polymake.call_function(:User, :set_shell_scalar, rstr*"_lp", lp)
  out, err = Polymake.shell_execute("""polytope::poly2mps(\$$(rstr)_poly, \$$(rstr)_lp);""")
  write(io, out)
  Polymake.call_function(:User, :set_shell_scalar, rstr*"_poly", 0)
  Polymake.call_function(:User, :set_shell_scalar, rstr*"_lp", 0)
  nothing
end

function _internal_save_lp(filename::String, p::Polymake.BigObject, lp::Polymake.BigObject, max::Bool)
  Polymake.polytope.poly2lp(p, lp, max, filename)
  nothing
end

function _internal_save_mps(filename::String, p::Polymake.BigObject, lp::Polymake.BigObject)
  Polymake.polytope.poly2mps(p, lp, Polymake.Set{Int}(), filename)
  nothing
end



@doc Markdown.doc"""
    save_lp(target::Union{String, IO}, lp::Union{MixedIntegerLinearProgram,LinearProgram})

Save a (mixed integer) linear program to an `.lp` file using the LP file format.

# Examples
Take the square $[-1/2,3/2]^2$ with objective $[1,1]$ and one integer variable.
Print the object in LP format to stdout:

```jldoctest
julia> c = cube(2, -1//2, 3//2)
A polyhedron in ambient dimension 2

julia> milp = MixedIntegerLinearProgram(c, [1,1], integer_variables=[1])
A mixed integer linear program

julia> save_lp(stdout, milp)
MAXIMIZE
  obj: +1 x1 +1 x2
Subject To
  ie0: +2 x1 >= -1
  ie1: -2 x1 >= -3
  ie2: +2 x2 >= -1
  ie3: -2 x2 >= -3
BOUNDS
  x1 free
  x2 free
GENERAL
  x1
END

"""
function save_lp(target::Union{String,IO}, lp::Union{MixedIntegerLinearProgram{fmpq},LinearProgram{fmpq}})
  _internal_save_lp(target,
                    pm_object(feasible_region(lp)),
                    pm_object(lp),
                    lp.convention == :max)
end

@doc Markdown.doc"""
    save_mps(target::String, lp::Union{MixedIntegerLinearProgram,LinearProgram})

Save a (mixed integer) linear program to an `.mps` file using the MPS file format.

# Examples
Create the square $[-1/2,3/2]^2$ with objective $[1,1]$ and one integer variable.
Print the object in MPS format to stdout:

```jldoctest
julia> c = cube(2, -1//2, 3//2)
A polyhedron in ambient dimension 2

julia> milp = MixedIntegerLinearProgram(c, [1,1], integer_variables=[1])
A mixed integer linear program

julia> save_mps(stdout, milp)
* Class:	MIP
* Rows:		5
* Columns:	2
* Format:	MPS
*
Name          unnamed#0
ROWS
 N  C0000000
 G  R0000000
 G  R0000001
 G  R0000002
 G  R0000003
COLUMNS
    M0000000  'MARKER'                 'INTORG'
    x1        C0000000  1                        R0000000  1
    x1        R0000001  -1                       
    M0000000  'MARKER'                 'INTEND'
    x2        C0000000  1                        R0000002  1
    x2        R0000003  -1                       
RHS
    B         R0000000  -0.5                     R0000001  -1.5
    B         R0000002  -0.5                     R0000003  -1.5
BOUNDS
 FR BND       x1  
 FR BND       x2  
ENDATA

"""
function save_mps(target::Union{String,IO}, lp::Union{MixedIntegerLinearProgram{fmpq},LinearProgram{fmpq}})
  _internal_save_mps(target,
                     pm_object(feasible_region(lp)),
                     pm_object(lp))
end

@doc Markdown.doc"""
    load_mps(file::String)

Load a (mixed integer) linear program from an `.mps` file using the MPS file format.
"""
function load_mps(file::String)
  poly = Polymake.polytope.mps2poly(file)::Polymake.BigObject
  if Polymake.exists(poly, "LP")
    lp = poly.LP
    Polymake.attach(lp, "convention", "max")
    return LinearProgram(Polyhedron(poly), lp, :max)
  elseif Polymake.exists(poly, "MILP")
    milp = poly.MILP
    Polymake.attach(milp, "convention", "max")
    return MixedIntegerLinearProgram(Polyhedron(poly), milp, :max)
  else
    throw(ErrorException("load_mps: cannot find LP or MILP subobject in polymake object"))
  end
end


@doc Markdown.doc"""
    load_lp(file::String)

Load a (mixed integer) linear program from an `.lp` file using the LP file format.
"""
function load_lp(file::String)
  poly = Polymake.polytope.lp2poly(file; nocheck=true)::Polymake.BigObject
  if Polymake.exists(poly, "LP")
    lp = poly.LP
    convention = occursin("MAXIMIZE", String(Polymake.getdescription(lp))) ? :max : :min
    Polymake.attach(lp, "convention", String(convention))
    return LinearProgram(Polyhedron(poly), lp, convention)
  elseif Polymake.exists(poly, "MILP")
    milp = poly.MILP
    convention = occursin("MAXIMIZE", String(Polymake.getdescription(milp))) ? :max : :min
    Polymake.attach(milp, "convention", String(convention))
    return MixedIntegerLinearProgram(Polyhedron(poly), milp, convention)
  else
    throw(ErrorException("load_lp: cannot find LP or MILP subobject in polymake object"))
  end
end

