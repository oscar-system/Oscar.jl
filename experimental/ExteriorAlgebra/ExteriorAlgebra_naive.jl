export exterior_algebra_naive;  # constructor
export exterior_algebra_singular;

#export ContainsDuplicates;   # export just for developing/testing !!!!

#--------------------------------------------

# ContainsDuplicates is auxiliary function -- worth keeping?
# The application here is to check whether the supplied variable names are distinct.

"""
    ContainsDuplicates(L)

Assumes `L` is indexable from 1 to `length(L)`.
Assumes the elements of `L` are ordered -- internally uses `sort`.
Returns `true` if `L` contains 2 elements which test as equal; otherwise `false`.
"""
function ContainsDuplicates(L)  # L must be indexable (1:n), and its elems must be sortable
  n = length(L);
  if (n < 2)  return false;  end;
  if (n == 2)  return (L[1] == L[2]);  end;
  sorted = sort(L);
  for i=1:n-1  if (sorted[i] == sorted[i+1])  return true;  end;  end;
  return false;
end;

### BUGS and DEFICIENCIES
# (1)  Requires copying and sorting: wastes space, values must be sortable!
# (2)  Arg type is currently unspecified -- should it be Vector{Any} ???


#--------------------------------------------

# exterior_algebra constructor:
#  args should be underlying CoeffRing and number of indets (or list of indet names)
#  optional keyword arg: indet_prefix (must be roman-alphabetic"


# Returns 2 components: ExtAlg, list of the gens/variables in order (e1,..,en)
# [[depends on what pbw_algebra returns!]]
@doc Markdown.doc"""
    exterior_algebra_naive(CoeffRing, NumVars)
    exterior_algebra_naive(CoeffRing, NumVars; var_prefix="x")
    exterior_algebra_naive(CoeffRing, ListOfVarNames)

First two forms return an exterior algebra with given `CoeffRing` and number of variables.
By default the variables are called `e1, e2, ...`; the optional keyword parameter can
specify an alphabetic prefix to use instead of `e`.

The third form returns an exterior algebra with given `CoeffRing`, and variables named
as specified in `ListOfVarNames` (which must be non-empty).

**Note:** in the first two forms the parameter `NumVars` must be positive; be aware that
large values will create an object occupying a lot of memory (probably cubic in `NumVars`).

# Examples
```jldoctest
julia> ExtAlg,(e1,e2) = exterior_algebra_naive(QQ, 2);

julia> e2*e1
-e1*e2

julia> is_zero((e1+e2)^2)
true

julia> ExtAlg,(x1,x2) = exterior_algebra_naive(QQ, 2; var_prefix="x");

julia> x2*x1
-x1*x2

julia> ExtAlg,(x,y) = exterior_algebra_naive(QQ, ["x","y"]);

julia> y*x
-x*y
```
"""
function exterior_algebra_naive(K::Ring, NumVars::Int64;  var_prefix="e")
  fn_name = "exterior_algebra_naive";
  if (NumVars < 1)  throw(ArgumentError("$fn_name ctor: (arg2) NumVars must be strictly positive"));  end;
  if (var_prefix == "" || !isascii(var_prefix) || findfirst(!isletter,var_prefix) != nothing)
    throw(ArgumentError("$fn_name ctor: var_prefix must be (roman)alphabetic"));
  end;
  return exterior_algebra_naive(K, (1:NumVars) .|> (k -> "$var_prefix$k"));
end

function exterior_algebra_naive(K::Ring, xs::Union{AbstractVector{<:AbstractString}, AbstractVector{Symbol}, AbstractVector{Char}})
    fn_name = "exterior_algebra_naive";
    NumVars = length(xs)
    if (NumVars == 0)  throw(ArgumentError("$fn_name: (arg2) no variables/indeterminates given"));  end;
    if (ContainsDuplicates(xs))  throw(ArgumentError("$fn_name: (arg2) variable names must be distinct"));  end;
    R, indets = PolynomialRing(K, xs);  # ordering??
    M = zero_matrix(R,NumVars,NumVars);
    for i=1:NumVars-1  for j=i+1:NumVars  M[i,j] = -indets[i]*indets[j];  end;  end;
    PBW, PBW_indets = pbw_algebra(R, M, degrevlex(indets); check=false); # disable check since we know it is OK!
    I = two_sided_ideal(PBW, PBW_indets.^2);
    ExtAlg,QuoMap = quo(PBW, I);  ## does it make sense to cache here???
    return ExtAlg, QuoMap.(PBW_indets);
##???  return ExtAlg, QuoMap.(PBW_indets)::Vector{PBWAlgQuoElem{T1, T2}} where { T1 < :Any, T2 < :Any};
end;



# BUGS/DEFICIENCIES (2023-02-13):
# (1)  By default chooses names "e1", "e2", ... up to "e$n";
#      with 3rd arg the prefix "e" is replaced by the specified var_prefix.
# (2)  Computations with elements do not automatically
#      reduce modulo the squares of the generators.
# (3)  Do we want a special printing function?  (show/display)
# (4)  Should the restriction to roman-alphabetic prefixes be slackened?


### EXTERNAL FILES:
#   * tests in Oscar.jl/test/Experimental/ExteriorAlgebra-test.jl


# -------------------------------------------------------
# Exterior algebra: delegating everything to Singular.

# Auxiliary function: since Singular.PolynomialRing wants a list of strings

"""
    var_names_for_singular(V)
    
Return a `Vector{String}` by conversion from an `AbstractVector`  of
`AbstractString`  or  `Symbol` or  `Char`.

NEEDED ONLY to prepare input for `Singular.PolynomialRing`;
maybe later `Singular.PolynomialRing` will accept more general vectors
for its second argument, and this function becomes redundant.
"""
function var_names_for_singular(xs::Union{AbstractVector{<:AbstractString},
                                          AbstractVector{Symbol},
                                          AbstractVector{Char}})
   NumVars = length(xs);
##   L = [""  for _ in 1:NumVars];
##   for i = 1:NumVars   L[i] = string(xs[i]);  end;
   L = String[];
   for var = xs   push!(L, string(var));  end;
   return L;
end;


#----------------------

# exterior_algebra constructor:
#  args should be underlying CoeffRing and number of indets (or list of indet names)
#  optional keyword arg: indet_prefix (must be roman-alphabetic"


# Returns 2 components: ExtAlg, list of the gens/variables in order (e1,..,en)
@doc Markdown.doc"""
    exterior_algebra_singular(CoeffRing, NumVars)
    exterior_algebra_singular(CoeffRing, NumVars; var_prefix="x")
    exterior_algebra_singular(CoeffRing, ListOfVarNames)

First two forms return an exterior algebra with given `CoeffRing` and number of variables.
By default the variables are called `e1, e2, ...`; the optional keyword parameter can
specify an alphabetic prefix to use instead of `e`.

The third form returns an exterior algebra with given `CoeffRing`, and variables named
as specified in `ListOfVarNames` (which must be non-empty).

**Note:** in the first two forms the parameter `NumVars` must be positive; be aware that
large values will create an object occupying a lot of memory (probably cubic in `NumVars`).

# Examples
```jldoctest
julia> ExtAlg,(e1,e2) = exterior_algebra_singular(QQ, 2);

julia> e2*e1
-e1*e2

julia> (e1+e2)^2  # result is automatically reduced!
0

julia> is_zero((e1+e2)^2)
true

julia> ExtAlg,(x1,x2) = exterior_algebra_singular(QQ, 2; var_prefix="x");

julia> x2*x1
-x1*x2

julia> ExtAlg,(x,y) = exterior_algebra_singular(QQ, ["x","y"]);

julia> y*x
-x*y
```
"""
function exterior_algebra_singular(K::Ring, NumVars::Int64;  var_prefix="e")
  fn_name = "exterior_algebra_singular";
  if (NumVars < 1)  throw(ArgumentError("$fn_name ctor: (arg2) NumVars must be strictly positive"));  end;
  if (var_prefix == "" || !isascii(var_prefix) || findfirst(!isletter,var_prefix) != nothing)
    throw(ArgumentError("$fn_name ctor: var_prefix must be (roman)alphabetic"));
  end;
  return exterior_algebra_singular(K, (1:NumVars) .|> (k -> "$var_prefix$k"));
end

function exterior_algebra_singular(K::AbstractAlgebra.Ring, xs::Union{AbstractVector{<:AbstractString}, AbstractVector{Symbol}, AbstractVector{Char}})
    fn_name = "exterior_algebra_singular";
    NumVars = length(xs)
    if (NumVars == 0)
        throw(ArgumentError("$fn_name: (arg2) no variables/indeterminates given"));
    end;
    VarNames = var_names_for_singular(xs);
    if (ContainsDuplicates(VarNames))  throw(ArgumentError("$fn_name: (arg2) variable names must be distinct"));  end;
    P, _ = Singular.PolynomialRing(K, VarNames)
    SINGULAR_PTR = Singular.libSingular.exteriorAlgebra(Singular.libSingular.rCopy(P.ptr));
    # what should I do with the pointer??
    ExtAlg = Singular.create_ring_from_singular_ring(SINGULAR_PTR);
    return ExtAlg, gens(ExtAlg);
end;
