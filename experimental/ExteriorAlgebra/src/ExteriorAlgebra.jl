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

@doc raw"""
    exterior_algebra(K::Field, numVars::Int)
    exterior_algebra(K::Field, listOfVarNames::AbstractVector{<:VarName})

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
  return exterior_algebra(K, (k -> "e$k").((1:numVars)))
end

#---------------------------------
# Method where caller specifies name of variables.

function exterior_algebra(K::Field, listOfVarNames::AbstractVector{<:VarName})
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
  for i in 1:(numVars - 1)
    for j in (i + 1):numVars
      M[i, j] = -indets[i] * indets[j]
    end
  end
  PBW, PBW_indets = pbw_algebra(R, M, degrevlex(indets); check=false) # disable check since we know it is OK!
  I = two_sided_ideal(PBW, PBW_indets .^ 2)
  # Now construct the fast exteriorAlgebra in Singular;
  # get var names from PBW in case it had "mangled" them.
  P, _ = Singular.polynomial_ring(SameCoeffRing, string.(symbols(PBW)))
  SINGULAR_PTR = Singular.libSingular.exteriorAlgebra(Singular.libSingular.rCopy(P.ptr))
  ExtAlg_singular = Singular.create_ring_from_singular_ring(SINGULAR_PTR)
  # Create Quotient ring with special implementation:
  ExtAlg, _ = quo(PBW, I; SpecialImpl=ExtAlg_singular)  # 2nd result is a QuoMap, apparently not needed
  ######    set_attribute!(ExtAlg, :is_exterior_algebra, :true)  ### DID NOT WORK (see PBWAlgebraQuo.jl)  Anyway, the have_special_impl function suffices.
  return ExtAlg, gens(ExtAlg)
end

# COMMENTED OUT "OLD IMPLEMENTATION" (so as not to lose the code)

# #--------------------------------------------
# # Exterior algebra implementation as a quotient of a PBW algebra;
# # **PREFER** exterior_algebra over this SLOW implementation!

# # Returns 2 components: ExtAlg, list of the gens/variables in order (e1,..,en)

# @doc raw"""
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

# function exterior_algebra_PBWAlgQuo(K::Ring, listOfVarNames::AbstractVector{<:VarName})
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


### Manual realization of exterior algebras

include("Types.jl")
function ExtAlgElem(E::ExteriorAlgebra{T}, a::Vector{Tuple{Int, SRow{T}}}; check::Bool=true) where {T}
  return ExtAlgElem(E, Dict{Int, SRow{T}}(i => v for (i, v) in a if !iszero(v)))
end

function ExtAlgElem(E::ExteriorAlgebra{T}, p::Int, v::SRow{T}; check::Bool=true) where {T}
  return ExtAlgElem(E, Dict{Int, SRow{T}}(i => v for (i, v) in a if !iszero(v)))
end

function mul!(a::T, w::ExtAlgElem{T}) where {T}
  @assert parent(a) === base_ring(parent(w))
  for i in keys(components(w))
    components(w)[i]*=a
    iszero(components(w)[i]) && delete!(components(w), i)
  end
  return w
end

function *(a::T, w::ExtAlgElem{T}) where {T <:NCRingElem}
  @assert parent(a) === base_ring(parent(w))
  nc = [(i, a*v) for (i, v) in components(w)]
  return ExtAlgElem(parent(w),
                    Dict{Int, SRow{T}}(i => v for (i, v) in nc if !iszero(v));
                    check=false
                   )
end

function *(a, w::ExtAlgElem{T}) where {T<:NCRingElem}
  return base_ring(w)(a)*w
end

function *(a::NCRingElement, w::ExtAlgElem{T}) where {T<:NCRingElem}
  return base_ring(w)(a)*w
end

function add!(v::ExtAlgElem{T}, w::ExtAlgElem{T}) where {T <: NCRingElem}
  for (i, a) in components(w)
    if i in keys(components(v))
      components(v)[i] += components(w)[i]
      iszero(components(v)[i]) && delete!(components(v), i)
    else
      components(v)[i] = components(w)[i]
    end
  end
  return v
end

function Base.deepcopy_internal(w::ExtAlgElem, id::IdDict)
  return ExtAlgElem(parent(w), deepcopy_internal(components(w), id); check=false)
end

function -(v::ExtAlgElem{T}, w::ExtAlgElem{T}) where {T}
  return v + (-w)
end

function +(v::ExtAlgElem{T}, w::ExtAlgElem{T}) where {T}
  return add!(copy(v), w)
end

function -(w::ExtAlgElem)
  return -one(base_ring(w))*w
end

base_ring(w::ExtAlgElem) = base_ring(parent(w))

parent(w::ExtAlgElem) = w.parent

rank(E::ExteriorAlgebra) = E.rank
base_ring(E::ExteriorAlgebra) = E.base_ring

zero(E::ExteriorAlgebra{T}) where {T} = ExtAlgElem(E)
zero(w::ExtAlgElem) = zero(parent(w))

function components(w::ExtAlgElem{T}) where {T}
  if !isdefined(w, :components)
    w.components = Dict{Int, SRow{T}}()
  end
  return w.components
end

function *(v::ExtAlgElem{T}, w::ExtAlgElem{T}) where {T<:NCRingElem}
  E = parent(v)
  R = base_ring(E)
  n = rank(E)
  result = zero(E)
  new_elem = Dict{Int, T}()
  for (q, b) in components(w)
    for (p, a) in components(v)
      r = p + q
      r > n && continue
      #new_elem = sizehint!(Tuple{Int, T}[], length(a)*length(b))
      empty!(new_elem)
      ht = multiplication_hash_table(E, p, q)
      for (i, c) in a
        for (j, d) in b
          s, k = ht[i, j]
          is_zero(s) && continue
          res = s * c * d
          if !iszero(res)
            if haskey(new_elem, k)
              new_elem[k] += res
            else
              new_elem[k] = res
            end
          end
          #!iszero(res) && push!(new_elem, (k, res))
        end
      end
      for (i, v) in new_elem
        iszero(v) && delete!(new_elem, i)
      end
      if  !isempty(new_elem)
        if haskey(components(result), r)
          components(result)[r] += sparse_row(R, [(i, v) for (i, v) in new_elem])
        else
          components(result)[r] = sparse_row(R, [(i, v) for (i, v) in new_elem])
        end
      end
    end
  end
  for (i, v) in components(result)
    iszero(v) && delete!(components(result), i)
  end
  return result
end
          
function print_symbols(E::ExteriorAlgebra, p::Int)
  if !isdefined(E, :print_symbols)
    E.print_symbols = Dict{Int, Vector{Symbol}}()
  end
  roof = is_unicode_allowed() ? "∧" : "^"
  if !haskey(E.print_symbols, p)
    E.print_symbols[p] = [Symbol(join([string(E.symbols[i]) for i in indices(I)], roof)) for I in OrderedMultiIndexSet(p, rank(E))]
  end
  return E.print_symbols[p]
end

function gens(E::ExteriorAlgebra{T}) where {T}
  return [ExtAlgElem(E, Dict{Int, SRow{T}}(1 => sparse_row(base_ring(E), 
                                                           [(i, one(base_ring(E)))]))) 
          for i in 1:rank(E)]
end

one(E::ExteriorAlgebra{T}) where {T} = ExtAlgElem(E, Dict{Int, SRow{T}}(0 => sparse_row(base_ring(E), [(1, one(base_ring(E)))])); check=false)

function Base.show(io::IO, w::ExtAlgElem)
  E = parent(w)
  if is_empty(components(w))
    println(io, "0")
    #println(io, "0 : 0")
    return
  end

  parts = String[]
  for p in 0:rank(E)
    !haskey(components(w), p) && continue
    v = components(w)[p]
    if iszero(p)
      push!(parts, "$(first(v.values))")
      #print(io, "0 : $(first(v.values))")
      continue
    end
    #print(io, "$p : ")
    p_symb = print_symbols(E, p)
    push!(parts, join(["$(c)*"*string(p_symb[i]) for (i, c) in v], " + "))
    #println(io, join(["$(c)*"*string(p_symb[i]) for (i, c) in v], " + "))
  end
  print(io, join(parts, " + "))
end

function ==(v::ExtAlgElem, w::ExtAlgElem)
  E = parent(v)
  @assert E === parent(w)
  return components(v) == components(w)
end

function multiplication_hash_table(E::ExteriorAlgebra, p::Int, q::Int)
  if !isdefined(E, :multiplication_hash_tables)
    E.multiplication_hash_tables = Dict{Tuple{Int, Int}, Matrix{Tuple{Int, Int}}}()
  end
  n = rank(E)
  if !haskey(E.multiplication_hash_tables, (p, q))
    A = [_wedge(ordered_multi_index(i, p, n), ordered_multi_index(j, q, n)) for i in 1:binomial(n, p), j in 1:binomial(n, q)]
    B = [(s, linear_index(k)) for (s, k) in A]
    E.multiplication_hash_tables[(p, q)] = B
  end
  return E.multiplication_hash_tables[(p, q)]
end

ExteriorAlgebra(R::NCRing, n::Int) = ExteriorAlgebra(R, [Symbol("e$i") for i in 1:n])

is_graded(E::ExteriorAlgebra) = true
function grading_group(E::ExteriorAlgebra)
  if !isdefined(E, :grading_group)
    E.grading_group = free_abelian_group(1)
  end
  return E.grading_group
end

elem_type(::Type{ExteriorAlgebra{T}}) where {T} = ExtAlgElem{T}
elem_type(E::ExteriorAlgebra) = elem_type(typeof(E))
parent_type(::Type{ExtAlgElem{T}}) where {T} = ExteriorAlgebra{T}
parent_type(w::ExtAlgElem) = parent_type(typeof(w))

function Base.show(io::IO, E::ExteriorAlgebra)
  print(io, "exterior algebra over $(base_ring(E)) in "*join(string.(E.symbols), ", "))
end

function degree(w::ExtAlgElem)
  isempty(components(w)) && return nothing
  return maximum(collect(keys(components(w))))*grading_group(parent(w))[1]
end

function is_homogeneous(w::ExtAlgElem)
  return length(components(w)) < 2
end

function (E::ExteriorAlgebra)(a::NCRingElem)
  return E(base_ring(E)(a))
end

function (E::ExteriorAlgebra{T})(a::ExtAlgElem{T}) where {T}
  @assert parent(a) === E
  return a
end

function (E::ExteriorAlgebra{T})(a::T) where {T}
  return a*one(E)
end

function (E::ExteriorAlgebra)(a::Int)
  return a*one(E)
end

function (E::ExteriorAlgebra)()
  return zero(E)
end

function _degree_fast(w::ExtAlgElem)
  return degree(w)
end

function getindex(E::ExteriorAlgebra, i::Int)
  if !isdefined(E, :graded_parts)
    E.graded_parts = Dict{Int, FreeMod}()
  end
  parts = E.graded_parts
  R = base_ring(E)
  n = rank(E)
  if !haskey(parts, i)
    if i == 1
      F = FreeMod(R, E.symbols)
      F.S = E.symbols
      parts[i] = F
    elseif (i >= 0 && i <= n)
      parts[i] = exterior_power(E[1], i)[1]
    else
      parts[i] = FreeMod(R, 0)
    end
  end
  return parts[i]::FreeMod{elem_type(R)}
end

@attr Dict{Int, FreeMod{T}} function _graded_parts(F::FreeMod{ExtAlgElem{T}}) where T
  return Dict{Int, FreeMod{T}}()
end

function _graded_part(F::FreeMod{ExtAlgElem{T}}, p::Int) where T
  graded_parts = _graded_parts(F)
  if !haskey(graded_parts, p)
    E = base_ring(F)
    R = base_ring(E)
    parts = [E[p - Int(degree(F[k])[1])] for k in 1:ngens(F)]
    result = direct_sum(parts...)[1]
    graded_parts[p] = result
  end
  return graded_parts[p]
end


@attr Dict{Int, FreeModuleHom} function _graded_parts(phi::FreeModuleHom{ModuleType, ModuleType, Nothing}) where {ModuleType <: FreeMod{<:ExtAlgElem}}
  return Dict{Int, FreeModuleHom}()
end

function getindex(phi::FreeModuleHom{ModuleType, ModuleType, Nothing}, p::Int) where {ModuleType <: FreeMod{<:ExtAlgElem}}
  @assert is_homogeneous(phi) && iszero(degree(phi)) "morphisms must be homogeneous of degree zero"
  parts = _graded_parts(phi)
  if !haskey(parts, p)
    dom = domain(phi)
    cod = codomain(phi)

    dom_p = _graded_part(dom, p)
    cod_p = _graded_part(cod, p)
    #=
    img_gens = gens(dom_p)
    img_gens = [dom(g, p) for g in img_gens]
    img_gens = phi.(img_gens)
    =#
    img_gens = elem_type(cod_p)[cod_p(phi(dom(g, p))) for g in gens(dom_p)]
    phi_p = hom(dom_p, cod_p, img_gens)
    parts[p] = phi_p
  end
  return parts[p]::FreeModuleHom
end

function (F::FreeMod{ExtAlgElem{T}})(v::FreeModElem{T}, p::Int) where {T}
  E = base_ring(F)
  F_p = parent(v)
  R = base_ring(F_p)
  @assert F_p === _graded_part(F, p)
  g_E = gens(F)
  d = degree.(gens(F))
  prs = canonical_projections(F_p)
  v_comps = [p(v) for p in prs]
  d = [Int(a[1]) for a in degree.(gens(F))]
  v_comp_coords = [iszero(c.coords) ? E() : ExtAlgElem(E, Dict{Int, SRow{T}}([p-d[i]=>c.coords])) for (i, c) in enumerate(v_comps)]
  return sum(a*g for (a, g) in zip(v_comp_coords, g_E); init=zero(F))
end

function (F_p::FreeMod{T})(v::FreeModElem{ExtAlgElem{T}}; check::Bool=true) where {T}
  @check is_homogeneous(v) "elements must be homogeneous for conversion"
  iszero(v) && return zero(F_p)
  p = Int(degree(v)[1])
  F = parent(v)
  @assert F_p === _graded_part(F, p)
  c_v = coordinates(v)::SRow
  inj = canonical_injections(F_p)
  result = zero(F_p)
  d = [Int(a[1]) for a in degree.(gens(F))]
  G = grading_group(F)
  for (i, c) in c_v
    inc = inj[i]
    v_i = domain(inc)(c.components[p - d[i]])
    result += inc(v_i)
  end
  return result
end

