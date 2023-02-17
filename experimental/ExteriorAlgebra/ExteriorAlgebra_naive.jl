export exterior_algebra_naive;  # constructor
export exterior_algebra_singular;

#--------------------------------------------


# Returns 2 components: ExtAlg, list of the gens/variables in order (e1,..,en)

Markdown.@doc doc"""
    exterior_algebra_naive(CoeffRing::Ring, NumVars::Int)
    exterior_algebra_naive(CoeffRing::Ring, ListOfVarNames::Vector{String})

The first form returns an exterior algebra with given `CoeffRing` and `NumVars` variables;
the variables are called `e1, e2, ...`.  The value `NumVars` must be positive; be aware that
large values will create an object occupying a lot of memory (probably cubic in `NumVars`).

The second form returns an exterior algebra with given `CoeffRing`, and variables named
as specified in `ListOfVarNames` (which must be non-empty).


# Examples
```jldoctest
julia> ExtAlg, (e1,e2)  =  exterior_algebra_naive(QQ, 2);

julia> e2*e1
-e1*e2

julia> is_zero((e1+e2)^2)
true

julia> ExtAlg, (x,y)  =  exterior_algebra_naive(QQ, ["x","y"]);

julia> y*x
-x*y
```
"""
function exterior_algebra_naive(K::Ring, NumVars::Int)
  fn_name = "exterior_algebra_naive";
  if (NumVars < 1)  throw(ArgumentError("$fn_name ctor: (arg2) NumVars must be strictly positive"));  end;
  return exterior_algebra_naive(K, (1:NumVars) .|> (k -> "e$k"));
end

function exterior_algebra_naive(K::Ring, xs::Union{AbstractVector{<:AbstractString}, AbstractVector{Symbol}, AbstractVector{Char}})
    fn_name = "exterior_algebra_naive";
    NumVars = length(xs)
    if (NumVars == 0)  throw(ArgumentError("$fn_name: (arg2) no variables/indeterminates given"));  end;
    if (!allunique(xs))  throw(ArgumentError("$fn_name: (arg2) variable names must be distinct"));  end;
    R, indets = PolynomialRing(K, xs);
    M = zero_matrix(R, NumVars, NumVars);
    for i in 1:NumVars-1  for j in i+1:NumVars  M[i,j] = -indets[i]*indets[j];  end;  end;
    PBW, PBW_indets = pbw_algebra(R, M, degrevlex(indets); check=false); # disable check since we know it is OK!
    I = two_sided_ideal(PBW, PBW_indets.^2);
    ExtAlg,QuoMap = quo(PBW, I);  ## does it make sense to cache here???
    return ExtAlg, QuoMap.(PBW_indets);
end;



# BUGS/DEFICIENCIES (2023-02-13):
# (1)  Computations with elements DO NOT AUTOMATICALLY REDUCE
#      modulo the squares of the generators.
# (2)  Do we want/need a special printing function?  (show/display)



### EXTERNAL FILES:
#   * tests in Oscar.jl/test/Experimental/ExteriorAlgebra-test.jl


# -------------------------------------------------------
# Exterior algebra: delegating everything to Singular.


# This fn NOT EXPORTED:  auxiliary function since Singular.PolynomialRing wants a list of strings

"""
    var_names_for_singular(V)  where V is Vector of Strings, Symbols or Chars
    
Return a `Vector{String}` by converting the elements of `V` successively into strings.

NEEDED ONLY to prepare input for `Singular.PolynomialRing`;
maybe later `Singular.PolynomialRing` will accept more general vectors
for its second argument, and this function becomes redundant.
"""
function var_names_for_singular(xs::Union{AbstractVector{<:AbstractString},
                                          AbstractVector{Symbol},
                                          AbstractVector{Char}})
   NumVars = length(xs);
# This commented out code appears to be a bit slower.
##   L = [""  for _ in 1:NumVars];
##   for i = 1:NumVars   L[i] = string(xs[i]);  end;
   L = String[];
   for var in xs   push!(L, string(var));  end;
   return L;
end;


#----------------------

# exterior_algebra constructor:
#  args should be underlying CoeffRing and number of indets (or list of indet names)
#  optional keyword arg: indet_prefix (must be roman-alphabetic"


# Returns 2 components: ExtAlg, list of the gens/variables in order (e1,..,en)
Markdown.@doc doc"""
    exterior_algebra_singular(CoeffRing::Field, NumVars::Int)
    exterior_algebra_singular(CoeffRing::Field, ListOfVarNames::Vector{String})

The first form returns an exterior algebra with given `CoeffRing` and `NumVars` variables;
the variables are called `e1, e2, ...`.  The value of `NumVars` must be positive; be aware that
large values will create an object occupying a lot of memory (probably cubic in `NumVars`).

The second form returns an exterior algebra with given `CoeffRing`, and variables named
as specified in `ListOfVarNames` (which must be non-empty).


# Examples
```jldoctest
julia> ExtAlg, (e1,e2)  =  exterior_algebra_singular(QQ, 2);

julia> e2*e1
-e1*e2

julia> (e1+e2)^2  # result is automatically reduced!
0

julia> is_zero((e1+e2)^2)
true

julia> ExtAlg, (x,y)  =  exterior_algebra_singular(QQ, ["x","y"]);

julia> y*x
-x*y
```
"""
function exterior_algebra_singular(K::Field, NumVars::Int)
  fn_name = "exterior_algebra_singular";
  if (NumVars < 1)  throw(ArgumentError("$fn_name ctor: (arg2) NumVars must be strictly positive"));  end;
  return exterior_algebra_singular(K, (1:NumVars) .|> (k -> "e$k"));
end

function exterior_algebra_singular(K::Field, xs::Union{AbstractVector{<:AbstractString}, AbstractVector{Symbol}, AbstractVector{Char}})
    fn_name = "exterior_algebra_singular";
    NumVars = length(xs)
    if (NumVars == 0)
        throw(ArgumentError("$fn_name: (arg2) no variables/indeterminates given"));
    end;
    VarNames = var_names_for_singular(xs);
    if (!allunique(VarNames))  throw(ArgumentError("$fn_name: (arg2) variable names must be distinct"));  end;
    P, _ = Singular.PolynomialRing(K, VarNames)
    SINGULAR_PTR = Singular.libSingular.exteriorAlgebra(Singular.libSingular.rCopy(P.ptr));
    ExtAlg = Singular.create_ring_from_singular_ring(SINGULAR_PTR);
    return ExtAlg, gens(ExtAlg);
end;
