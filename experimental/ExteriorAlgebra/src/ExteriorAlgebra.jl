export exterior_algebra

#  See also:
#   * tests in Oscar.jl/test/Experimental/ExteriorAlgebra-test.jl
#   * doc tests in Oscar.jl/docs/src/NoncommutativeAlgebra/PBWAlgebras/quotients.md

#---------------------- MAIN FUNCTION ----------------------

# Returns 2 components: ExtAlg, list of the gens/variables

# DEVELOPER DOC
#   This impl is "inefficient": it must create essentially 2 copies of
#   the exterior algebra.  One copy is so that Oscar knows the structure
#   of the ext alg (as a quotient of a PBWAlg); the other copy is a
#   Singular implementation (seen as a "black box" by Oscar) which
#   actually does the arithmetic (quickly).
#
#   To make this work I had to make changes to struct PBWAlgQuo: (see the source)
#   (*) previously a PBWAlgQuo had a single datum, namely the (2-sided) ideal
#       used to make the quotient -- the base ring could be derived from the ideal
# 
#   (*) now a PBWAlgQuo has an extra data field "sring" which refers to
#       the underlying Singular ring structure which actually performs
#       the arithmetic.  For exterior algebras "sring" refers to a
#       specific Singular ring "exteriorAlgebra"; in other cases "sring"
#       refers to the Singular(plural) ring which implements the PBW alg
#       (i.e. IGNORING the fact that the elems are in a quotient)

#   PBWAlgQuoElem did not need to change.  Its "data" field just refers
#   to a PBWAlgElem (namely some representative of the class).

@doc raw"""
    exterior_algebra(K::Ring, nvars::Int)
    exterior_algebra(K::Ring, varnames::AbstractVector{<:VarName})

Given a coefficient ring `K` and variable names, say `varnames = [:x1, :x2, ...]`,
return a tuple `E, [x1, x2, ...]` consisting of the exterior algebra `E` over
the polynomial ring `R[x1, x2, ...]` and its generators `x1, x2, ...`.

If `K` is a field, this function will use a special implementation in Singular.

!!! note
    Creating an `exterior_algebra` with many variables will create an object
    occupying a lot of memory (probably cubic in `nvars`).

# Examples
```jldoctest
julia> E, (x1,x2)  =  exterior_algebra(QQ, 2);

julia> x2*x1
-x1*x2

julia> (x1+x2)^2  # over fields, result is automatically reduced!
0

julia> E, (x,y)  =  exterior_algebra(QQ, ["x","y"]);

julia> y*x
-x*y
```
"""
function exterior_algebra(K::Ring, varnames::Vector{Symbol})
  @req !isempty(varnames) "no variables given"
  numVars = length(varnames)

  R, indets = polynomial_ring(K, varnames; cached=false)
  M = zero_matrix(R, numVars, numVars)
  for i in 1:(numVars - 1)
    for j in (i + 1):numVars
      M[i, j] = -indets[i] * indets[j]
    end
  end
  PBW, PBW_indets = pbw_algebra(R, M, degrevlex(indets); check=false) # disable check since we know it is OK!
  I = two_sided_ideal(PBW, PBW_indets .^ 2)

  if K isa Field
    # Now construct the fast exteriorAlgebra in Singular;
    # get var names from PBW in case it had "mangled" them.
    K_singular = singular_coeff_ring(coefficient_ring(R))
    R_singular, _ = Singular.polynomial_ring(K_singular, symbols(PBW))
    SINGULAR_PTR = Singular.libSingular.exteriorAlgebra(
      Singular.libSingular.rCopy(R_singular.ptr)
    )
    ExtAlg_singular = Singular.create_ring_from_singular_ring(SINGULAR_PTR)
    # Create Quotient ring with special implementation:
    ExtAlg, _ = quo(PBW, I; special_impl=ExtAlg_singular)  # 2nd result is a QuoMap, apparently not needed
    generators = gens(ExtAlg)
  else
    ExtAlg, QuoMap = quo(PBW, I)
    generators = QuoMap.(PBW_indets)
  end
  set_attribute!(ExtAlg, :show, show_exterior_algebra)
  return ExtAlg, generators
end

AbstractAlgebra.@varnames_interface exterior_algebra(K::Ring, varnames) macros = :no

function show_exterior_algebra(io::IO, E::PBWAlgQuo)
  x = symbols(E)
  io = pretty(io)
  print(io, "Exterior algebra over ", Lowercase(), coefficient_ring(E))
  print(io, " in (")
  join(io, x, ", ")
  print(io, ")")
end

# # BUGS/DEFICIENCIES (2023-02-13):
# # (1)  Computations with elements DO NOT AUTOMATICALLY REDUCE
# #      modulo the squares of the generators.
# # (2)  Do we want/need a special printing function?  (show/display)
