export exterior_algebra  # MAIN EXPORT!

# Commented out old impl:  exterior_algebra_PBWAlgQuo  (allows coeffs in a non-field)


#--------------------------------------------
# Two implementations of exterior algebras:
# (1) delegating everything to Singular  -- fast, coeff ring must be a field
# (2) as a quotient of PBW algebra       -- slower, coeff ring must be commutative

# ADVICE: avoid impl (2) which deliberately has an awkward name.

#  See also:
#   * tests in Oscar.jl/test/Experimental/ExteriorAlgebra-test.jl
#   * doc tests in Oscar.jl/docs/src/NoncommutativeAlgebra/PBWAlgebras/quotients.md

# -------------------------------------------------------
# Exterior algebra: delegating everything to Singular.


#---------------------- MAIN FUNCTION ----------------------

# exterior_algebra constructor:  args are
#  - underlying coeff FIELD and
#  - number of indets (or list of indet names)

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


# Attach docstring to "abstract" function exterior_algebra, so that
# it is automatically "inherited" by the methods.

@doc Markdown.doc"""
    exterior_algebra(K::Field, numVars::Int)
    exterior_algebra(K::Field, listOfVarNames::Union{AbstractVector{<:AbstractString},
                                                     AbstractVector{Symbol},
                                                     AbstractVector{Char}})

The *first form* returns an exterior algebra with coefficient field `K` and
`numVars` variables: `numVars` must be positive, and the variables are
 called `e1, e2, ...`.

The *second form* returns an exterior algebra with coefficient field `K`, and
variables named as specified in `listOfVarNames` (which must be non-empty).

NOTE: Creating an `exterior_algebra` with many variables will create an object
occupying a lot of memory (probably cubic in `numVars`).


# Examples
```jldoctest
julia> ExtAlg, (e1,e2)  =  exterior_algebra(QQ, 2);

julia> e2*e1
-e1*e2

julia> (e1+e2)^2  # result is automatically reduced!
0

julia> ExtAlg, (x,y)  =  exterior_algebra(QQ, ["x","y"]);

julia> y*x
-x*y
```
"""
function exterior_algebra end


# ---------------------------------
# -- Method where caller specifies just number of variables

function exterior_algebra(K::Field, numVars::Int)
    if numVars < 1
        throw(ArgumentError("numVars must be strictly positive, but numVars=$numVars"))
    end
    return exterior_algebra(K,  (1:numVars) .|> (k -> "e$k"))
end

#---------------------------------
# Method where caller specifies name of variables.

function exterior_algebra(K::Field, listOfVarNames::Union{AbstractVector{<:AbstractString},
                                                          AbstractVector{Symbol},
                                                          AbstractVector{Char}})
    numVars = length(listOfVarNames)
    if numVars == 0
        throw(ArgumentError("no variables/indeterminates given"))
    end
#    if (!allunique(VarNames))
#        throw(ArgumentError("variable names must be distinct"))
#    end

    R, indets = polynomial_ring(K, listOfVarNames)
    SameCoeffRing = singular_coeff_ring(coefficient_ring(R))
    M = zero_matrix(R, numVars, numVars)
    for i in 1:numVars-1
        for j in i+1:numVars
            M[i,j] = -indets[i]*indets[j]
        end
    end
    PBW, PBW_indets = pbw_algebra(R, M, degrevlex(indets); check=false) # disable check since we know it is OK!
    I = two_sided_ideal(PBW, PBW_indets.^2)
    # Now construct the fast exteriorAlgebra in Singular; get var names from
    # PBW in case it had "mangled" them.
    P, _ = Singular.polynomial_ring(SameCoeffRing, string.(symbols(PBW)))
    SINGULAR_PTR = Singular.libSingular.exteriorAlgebra(Singular.libSingular.rCopy(P.ptr))
    ExtAlg_singular = Singular.create_ring_from_singular_ring(SINGULAR_PTR)
    # ***WORKAROUND***
    # Singular is "too smart" when there is just 1 indet, so we error out
    # When Singular is fixed, a @test_broken in the test suite will fail!!
    # When that happens remove this comment and the if stmt below (& fix the test)
    if supertype(typeof(ExtAlg_singular)) != AbstractAlgebra.NCRing
        throw(NotImplementedError(:exterior_algebra, "1 variable not yet supported  (requires Singular update)"))
    end
    # ***END OF WORKAROUND***
    # Create Quotient ring with special implementation:
    ExtAlg,_ = quo(PBW, I;  SpecialImpl = ExtAlg_singular)  # 2nd result is a QuoMap, apparently not needed
    return ExtAlg, gens(ExtAlg)
end



# COMMENTED OUT "OLD IMPLEMENTATION" (so as not to lose the code)

# #--------------------------------------------
# # Exterior algebra implementation as a quotient of a PBW algebra;
# # **PREFER** exterior_algebra over this SLOW implementation!

# # Returns 2 components: ExtAlg, list of the gens/variables in order (e1,..,en)


# @doc Markdown.doc"""
#     exterior_algebra_PBWAlgQuo(coeffRing::Ring, numVars::Int)
#     exterior_algebra_PBWAlgQuo(coeffRing::Ring, listOfVarNames::Vector{String})

# Use `exterior_algebra` in preference to this function when `coeffRing` is a field.

# The first form returns an exterior algebra with given `coeffRing` and `numVars` variables;
# the variables are called `e1, e2, ...`.  The value `numVars` must be positive; be aware that
# large values will create an object occupying a lot of memory (probably cubic in `numVars`).

# The second form returns an exterior algebra with given `coeffRing`, and variables named
# as specified in `listOfVarNames` (which must be non-empty).


# # Examples
# ```jldoctest
# julia> ExtAlg, (e1,e2)  =  exterior_algebra_PBWAlgQuo(QQ, 2);

# julia> e2*e1
# -e1*e2

# julia> is_zero((e1+e2)^2)
# true

# julia> ExtAlg, (x,y)  =  exterior_algebra_PBWAlgQuo(QQ, ["x","y"]);

# julia> y*x
# -x*y
# ```
# """
# function exterior_algebra_PBWAlgQuo(K::Ring, numVars::Int)
#     if (numVars < 1)
#         throw(ArgumentError("numVars must be strictly positive: numVars=$numVars"))
#     end
#     return exterior_algebra_PBWAlgQuo(K,  (1:numVars) .|> (k -> "e$k"))
# end

# function exterior_algebra_PBWAlgQuo(K::Ring, listOfVarNames::Union{AbstractVector{<:AbstractString}, AbstractVector{Symbol}, AbstractVector{Char}})
#     numVars = length(listOfVarNames)
#     if (numVars == 0)
#         throw(ArgumentError("no variables/indeterminates given"))
#     end
#     # if (!allunique(listOfVarNames))
#     #     throw(ArgumentError("variable names must be distinct"))
#     # end
#     R, indets = polynomial_ring(K, listOfVarNames)
#     M = zero_matrix(R, numVars, numVars)
#     for i in 1:numVars-1
#         for j in i+1:numVars
#             M[i,j] = -indets[i]*indets[j]
#         end
#     end
#     PBW, PBW_indets = pbw_algebra(R, M, degrevlex(indets);  check = false) # disable check since we know it is OK!
#     I = two_sided_ideal(PBW, PBW_indets.^2)
#     ExtAlg,QuoMap = quo(PBW, I)
#     return ExtAlg, QuoMap.(PBW_indets)
# end



# # BUGS/DEFICIENCIES (2023-02-13):
# # (1)  Computations with elements DO NOT AUTOMATICALLY REDUCE
# #      modulo the squares of the generators.
# # (2)  Do we want/need a special printing function?  (show/display)
