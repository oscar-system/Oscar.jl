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
  return add!(deepcopy(v), w)
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

function getindex(w::ExtAlgElem, p::Int)
  !haskey(w.components, p) && return sparse_row(base_ring(parent(w)))
  return w.components[p]
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
    print(io, "0")
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
  (!isdefined(v, :components) && !isdefined(w, :components)) && (return !isdefined(v, :components) && !isdefined(w, :components))
  #return components(v) == components(w)
  for (p, c) in components(v)
    p in keys(components(w)) || return false
    c == components(w)[p] || return false
  end
  return all(q in keys(components(v)) for q in keys(components(w)))
  for (q, cc) in components(w)
    q in keys(components(v)) || return false
  end
  return true
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

function degree(w::ExtAlgElem; check::Bool=true)
  (!isdefined(w, :components) || isempty(components(w))) && return zero(grading_group(parent(w)))
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

function (E::ExteriorAlgebra{T})(a::T) where {T <: NCRingElem}
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
    parts = [E[p - Int(degree(g)[1])] for g in gens(F)]
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
    #img_gens = elem_type(cod_p)[cod_p(phi(dom(g, p)); check=false) for g in gens(dom_p)]
    img_gens = cod_p(phi.(dom(gens(dom_p), p)); check=false)
    phi_p = hom(dom_p, cod_p, img_gens)
    parts[p] = phi_p
  end
  return parts[p]::FreeModuleHom
end

function (F::FreeMod{ExtAlgElem{T}})(v::Vector{FreeModElem{T}}, p::Int) where {T}
  E = base_ring(F)
  isempty(v) && return elem_type(F)[]
  F_p = parent(first(v))
  R = base_ring(F_p)
  @assert F_p === _graded_part(F, p)
  g_E = gens(F)
  #=
  # old code which is correct, but slow
  prs = canonical_projections(F_p)
  v_comp_list = [[p(w) for p in prs] for w in v]
  =#
  ranges = get_attribute(F_p, :ranges)::Vector{UnitRange{Int}}
  F_p_parts = get_attribute(F_p, :direct_product)::NTuple
  v_comp_list = [[G(coordinates(w)[r]) for (G, r) in zip(F_p_parts, ranges)] for w in v]
  #d = [Int(degree(g; check=false)[1]) for g in gens(F)]
  d = [Int(d[1]) for d in degrees_of_generators(F; check=false)]
  v_comp_coords_list = [[iszero(c.coords) ? E() : ExtAlgElem(E, Dict{Int, SRow{T}}([p-d[i]=>c.coords])) for (i, c) in enumerate(v_comps)] for v_comps in v_comp_list]
  return elem_type(F)[sum(a*g for (a, g) in zip(v_comp_coords, g_E); init=zero(F)) for v_comp_coords in v_comp_coords_list]
end

function (F::FreeMod{ExtAlgElem{T}})(v::FreeModElem{T}, p::Int) where {T}
  E = base_ring(F)
  F_p = parent(v)
  R = base_ring(F_p)
  @assert F_p === _graded_part(F, p)
  g_E = gens(F)
  #=
  prs = canonical_projections(F_p)
  v_comps = [p(v) for p in prs]
  =#
  ranges = get_attribute(F_p, :ranges)::Vector{UnitRange{Int}}
  F_p_parts = get_attribute(F_p, :direct_product)::NTuple
  v_comps = [G(coordinates(v)[r]) for (G, r) in zip(F_p_parts, ranges)]
  #d = [Int(degree(g; check=false)[1]) for g in gens(F)]
  d = [Int(d[1]) for d in degrees_of_generators(F; check=false)]
  v_comp_coords = [iszero(c.coords) ? E() : ExtAlgElem(E, Dict{Int, SRow{T}}([p-d[i]=>c.coords])) for (i, c) in enumerate(v_comps)]
  return sum(a*g for (a, g) in zip(v_comp_coords, g_E); init=zero(F))
end

function (F_p::FreeMod{T})(v::FreeModElem{ExtAlgElem{T}}; check::Bool=true) where {T}
  @check is_homogeneous(v) "elements must be homogeneous for conversion"
  iszero(v) && return zero(F_p)
  p = Int(degree(v; check=false)[1])
  F = parent(v)
  @assert F_p === _graded_part(F, p)
  c_v = coordinates(v)::SRow
  # Slow code:
  result = zero(F_p)
  #d = [Int(degree(g; check=false)[1]) for g in gens(F)]
  d = [Int(d[1]) for d in degrees_of_generators(F; check=false)]
  G = grading_group(F)
  ranges = get_attribute(F_p, :ranges)::Vector{UnitRange{Int}}
  for (i, c) in c_v
    w = deepcopy(c.components[p-d[i]])
    offset = first(ranges[i])-1
    w.pos.+=offset
    result += F_p(w)
  end
  return result
  # old code left here to show what the intention was.
  # the above is a hack to speed stuff up, but it messes with unstable internals
  inj = canonical_injections(F_p)
  for (i, c) in c_v
    inc = inj[i]
    v_i = domain(inc)(c.components[p - d[i]])
    result += inc(v_i)
  end
  return result
end

function (F_p::FreeMod{T})(a::Vector{FreeModElem{ExtAlgElem{T}}}; check::Bool=true) where {T}
  @check is_homogeneous(v) "elements must be homogeneous for conversion"
  isempty(a) && return elem_type(F_p)[]
  j = findfirst(!iszero(v) for v in a)
  j === nothing && return [zero(F_p) for i in 1:length(a)]
  p = Int(degree(a[j]; check=false)[1])
  F = parent(first(a))
  @assert F_p === _graded_part(F, p)
  c_v_list = [coordinates(v)::SRow for v in a]
  inj = canonical_injections(F_p)
  result = elem_type(F_p)[]
  #d = [Int(degree(g; check=false)[1]) for g in gens(F)]
  d = [Int(d[1]) for d in degrees_of_generators(F; check=false)]
  G = grading_group(F)
  for c_v in c_v_list
    res = zero(F_p)
    for (i, c) in c_v
      inc = inj[i]
      v_i = domain(inc)(c.components[p - d[i]])
      res += inc(v_i)
    end
    push!(result, res)
  end
  return result
end

@attr Dict{Int, SubModuleOfFreeModule{T}} function _graded_parts(I::SubModuleOfFreeModule{ExtAlgElem{T}}) where {T}
  return Dict{Int, SubModuleOfFreeModule{T}}()
end

function _graded_part(I::SubModuleOfFreeModule{ExtAlgElem{T}}, p::Int) where {T}
  graded_parts = _graded_parts(I)
  if !haskey(graded_parts, p)
    F = ambient_free_module(I)
    pp = grading_group(F)([p])
    E = base_ring(F)
    R = base_ring(E)
    iszero(F) && return SubModuleOfFreeModule(F, elem_type(F)[])
    #d = [Int(degree(g; check=false)[1]) for g in gens(F)]
    d = [Int(d[1]) for d in degrees_of_generators(F; check=false)]
    d0 = minimum(d)
    d_max = maximum(d)
    F_p = _graded_part(F, p)
    if p < d0 || p > d_max + rank(E)
      result = SubModuleOfFreeModule(F_p, elem_type(F_p)[])
      graded_parts[p] = result
      return result
    elseif p == d0
      g_0 = elem_type(F)[g for g in gens(I) if degree(g; check=false) == pp]
      g_0_p = elem_type(F_p)[F_p(g) for g in g_0]
      result = SubModuleOfFreeModule(F_p, g_0_p)
      set_attribute!(result, :_mapping_dict=>IdDict{elem_type(F_p), Tuple{elem_type(F), elem_type(E)}}(g_0_p[i] => (g_0[i], one(E)) for i in 1:length(g_0)))
      graded_parts[p] = result
      return result
    end

    # we can assume that d0 < p <= d_max + rank(E) and by induction 
    # that I_{p-1} has already been computed
    F_p = _graded_part(F, p)
    gens_p = elem_type(F_p)[]
    I_q = _graded_part(I, p-1)
    f = gens(I_q)
    # mapping the generators of F_p to ω ⋅ g for generators g of I and monomials ω ∈ E
    map_dict_p = IdDict{elem_type(F_p), Tuple{elem_type(F), elem_type(E)}}()
    map_dict_q = _mapping_dict(I_q)
    mult_gens = elem_type(F)[]
    for f_q in gens(I_q)
      f, w = map_dict_q[f_q]
      wf = w*f
      for (i, a) in enumerate(gens(E))
        g = a*wf
        (iszero(g) || g in mult_gens || -g in mult_gens) && continue
        push!(mult_gens, g)
        aw = a*w
        g_p = F_p(g)
        push!(gens_p, g_p)
        map_dict_p[g_p] = (f, aw)
      end
    end
    gens_ext = [g for g in gens(I) if degree(g; check=false) == pp]
    gens_p_ext = elem_type(F_p)[F_p(g) for g in gens_ext]
    result = SubModuleOfFreeModule(F_p, vcat(gens_p, gens_p_ext))
    for i in 1:length(gens_ext)
      map_dict_p[gens_p_ext[i]] = (gens_ext[i], one(E))
    end
    set_attribute!(result, :_mapping_dict=>map_dict_p)
    graded_parts[p] = result
  end
  return graded_parts[p]
end

@attr IdDict function _mapping_dict(I::SubModuleOfFreeModule)
  error("this attribute needs to be set manually")
end

function in(v::FreeModElem{ExtAlgElem{T}}, I::SubModuleOfFreeModule{ExtAlgElem{T}}) where {T}
  @assert is_homogeneous(v)
  iszero(v) && return true
  p = Int(degree(v; check=false)[1])
  I_p = _graded_part(I, p)
  F = ambient_free_module(I)
  F_p = _graded_part(F, p)
  return F_p(v) in I_p
end

function coordinates(v::FreeModElem{ExtAlgElem{T}}, I::SubModuleOfFreeModule{ExtAlgElem{T}}) where {T}
  @assert is_homogeneous(v)
  F = ambient_free_module(I)
  @assert parent(v) === F
  E = base_ring(F)
  R = base_ring(E)
  iszero(v) && return sparse_row(E)
  p = Int(degree(v; check=false)[1])
  I_p = _graded_part(I, p)
  F_p = _graded_part(F, p)
  c = coordinates(F_p(v), I_p)
  map_dict = _mapping_dict(I_p)
  result_list = Tuple{Int, elem_type(E)}[]
  for (i, a) in c
    g_p = I_p[i]
    g, w = map_dict[g_p]
    j = findfirst(j->g===I[j], 1:ngens(I))
    j === nothing && error("generator not found")
    push!(result_list, (j, a*w))
  end
  return sparse_row(E, result_list)
end

@attr Dict{Int, SubModuleOfFreeModule{ExtAlgElem{T}}} function _kernel_parts(
    phi::FreeModuleHom{ModuleType, ModuleType, Nothing}
  ) where {T, ModuleType <: FreeMod{ExtAlgElem{T}}}
  return Dict{Int, SubModuleOfFreeModule{ExtAlgElem{T}}}()
end

function _kernel_part(
    phi::FreeModuleHom{ModuleType, ModuleType, Nothing}, p::Int
  ) where {T, ModuleType <: FreeMod{ExtAlgElem{T}}}
  kernel_parts = _kernel_parts(phi)
  @vprint :ExteriorAlgebras 2 "computing kernel up to degree $p\n"
  if !haskey(kernel_parts, p)
    F = domain(phi)
    F_p = _graded_part(F, p)
    iszero(F) && return SubModuleOfFreeModule(F, elem_type(F)[])
    gens_deg = [Int(degree(g; check=false)[1]) for g in gens(F)]
    d0 = minimum(gens_deg)
    E = base_ring(F)
    n = rank(E)
    dmax = maximum(gens_deg) + n

    if p < d0 || p > dmax
      result = SubModuleOfFreeModule(F, elem_type(F)[])
      kernel_parts[d0-1] = result
      @vprint :ExteriorAlgebras 2 "done with kernel computation up to degree $p\n"
      return result
    elseif p == d0
      @vprint :ExteriorAlgebras 3 "extracting the graded part of degree $p\n"
      phi_p = phi[p]
      @vprint :ExteriorAlgebras 3 "restricted map has $(ngens(domain(phi_p))) generators in the domain and $(ngens(codomain(phi_p))) in the codomain\n"
      K_p, _ = kernel(phi_p)
      gens_p = filter!(!iszero, ambient_representatives_generators(K_p))
      g = elem_type(F)[F(v, p) for v in gens_p]
      Z = SubModuleOfFreeModule(F, g)
      kernel_parts[p] = Z
      @vprint :ExteriorAlgebras 2 "done with kernel computation up to degree $p\n"
      return Z
    end

    @vprint :ExteriorAlgebras 3 "extracting the graded part of degree $p\n"
    phi_p = phi[p]
    K_p, _ = kernel(phi_p)
    @vprint :ExteriorAlgebras 3 "restricted map has $(ngens(domain(phi_p))) generators in the domain and $(ngens(codomain(phi_p))) in the codomain\n"
    B = _kernel_part(phi, p-1)
    B_p = _graded_part(B, p)
    if all(g in B_p for g in ambient_representatives_generators(K_p))
      kernel_parts[p] = B
      @vprint :ExteriorAlgebras 3 "no new generators\n"
      @vprint :ExteriorAlgebras 2 "done with kernel computation up to degree $p\n"
      return B
    end
    M = SubquoModule(K_p.sub, B_p)
    #M2, to_M = simplify_light(M)
    M2, to_M = prune_with_map(M)
    new_gens_p = ambient_representative.(to_M.(gens(M2)))
    @vprint :ExteriorAlgebras 3 "$(length(new_gens_p)) new generators\n"
    new_gens = elem_type(F)[F(v, p) for v in new_gens_p]
    Z = B + SubModuleOfFreeModule(F, new_gens)
    kernel_parts[p] = Z
  end
  @vprint :ExteriorAlgebras 2 "done with kernel computation up to degree $p\n"
  return kernel_parts[p]
end




function kernel(
    phi::FreeModuleHom{ModuleType, ModuleType, Nothing}; 
    degree_bound::Union{Int, Nothing}=0
  ) where {ModuleType <: FreeMod{<:ExtAlgElem}}
  @vprint :ExteriorAlgebras 1 "computing kernel of map with degrees\n"
  @vprint :ExteriorAlgebras 1 "$([degree(g; check=false) for g in gens(domain(phi))])\n"
  @vprint :ExteriorAlgebras 1 "in the domain and degrees\n"
  @vprint :ExteriorAlgebras 1 "$([degree(g; check=false) for g in gens(codomain(phi))])\n"
  @vprint :ExteriorAlgebras 1 "in the codomain up to degree $(degree_bound)\n"
  F = domain(phi)
  iszero(F) && return sub(F, elem_type(F)[])
  d0 = minimum([Int(degree(g; check=false)[1]) for g in gens(F)])
  E = base_ring(F)
  n = rank(E)
  dmax = maximum([Int(degree(g; check=false)[1]) for g in gens(F)]) + n
  if degree_bound !== nothing
    dmax > degree_bound && (dmax = degree_bound)
  end
  Z = _kernel_part(phi, dmax)
  return sub(F, gens(Z))
end

function change_base_ring(R::NCRing, E::ExteriorAlgebra{T}) where {T<:NCRingElem}
  @assert R === base_ring(E) "the ring needs to be the base ring of the exterior algebra"
  return R, x->x[0][1]
end

function change_base_ring(R::Ring, F::FreeMod{ExtAlgElem{T}}) where {T<:RingElem}
  E = base_ring(F)
  @assert R === base_ring(E) "the ring needs to be the base ring of the exterior algebra"

  indices = [i for i in 1:ngens(F) if iszero(degree(F[i]; check=false))]
  FoR = FreeMod(R, length(indices))
  #set_attribute!(FoR, :degrees_of_generators=>degrees_of_generators(F))
  _, red0 = change_base_ring(R, E)
  img_gens = elem_type(FoR)[]
  j = 1
  for (i, g) in enumerate(gens(F))
    if !iszero(degree(g; check=false))
      push!(img_gens, zero(FoR))
    else
      push!(img_gens, gen(FoR, j))
      j += 1
    end
  end
  #img_gens = [i in indices ? FoR[indices[i]] : zero(FoR) for i in 1:ngens(F)]
  red = hom(F, FoR, img_gens, red0)
  return FoR, red
end

function change_base_ring(R::Ring, f::ModuleFPHom{MT, MT, Nothing};
    domain_base_change::Map = change_base_ring(R, domain(f))[2],
    codomain_base_change::Map = change_base_ring(R, codomain(f))[2]
  ) where {MT <: ModuleFP{<:ExtAlgElem}}
  bc_dom = codomain(domain_base_change)
  bc_cod = codomain(codomain_base_change)
  @assert domain(f) === domain(domain_base_change)
  @assert codomain(f) === domain(codomain_base_change)
  indices = [i for i in 1:ngens(domain(f)) if iszero(degree(domain(f)[i]; check=false))]
  img_gens = gens(domain(f))[indices]
  img_gens = f.(img_gens)
  img_gens = codomain_base_change.(img_gens)
  img_gens = [codomain_base_change(f(domain(f)[i])) for i in indices]
  return hom(bc_dom, bc_cod, img_gens)
end

function _strand(F::FreeMod{T}, d) where {T<:MPolyDecRingElem}
  @assert is_z_graded(F)
  C = ZeroDimensionalComplex(F)
  C0, map_to_orig = strand(C, d)
  return C0[()], map_to_orig[()]
end

### Partially copied from the StrandComplexes. TODO: Clean up to avoid duplication!
function strand(F::FreeMod{<:MPolyDecRingElem}, d::Int)
  @assert is_z_graded(F)
  S = base_ring(F)
  R = base_ring(S)
  Fd = FreeMod(R, length(all_exponents(F, d)))
  
  # Use a dictionary for fast mapping of the monomials to the 
  # generators of `Fd`.
  inv_exp_dict = Dict{Tuple{Vector{Int}, Int}, elem_type(Fd)}(m=>Fd[k] for (k, m) in enumerate(all_exponents(F, d)))
  # Hashing of FreeModElem's can not be assumed to be non-trivial. Hence we use the exponents directly.
  img_gens = elem_type(F)[]
  x = gens(S)
  for (e, i) in all_exponents(F, d) # iterate through the generators of `F`
    m = prod(x^k for (x, k) in zip(x, e); init=one(S))*F[i]
    push!(img_gens, m)
  end
  Fd_to_F = hom(Fd, F, img_gens, S)
  function my_map(v::FreeModElem)
    iszero(v) && return zero(Fd)
    @assert Int(degree(v; check=false)[1]) == d "input must be homogeneous of degree $d"
    w = zero(Fd)
    for (i, b) in coordinates(v)
      w += sum(c*inv_exp_dict[(n, i)] for (c, n) in zip(coefficients(b), exponents(b)); init=zero(Fd))
    end
    return w
  end
  F_to_Fd = MapFromFunc(F, Fd, my_map)
  return Fd, Fd_to_F, F_to_Fd
end

function _ext_module_map(
    F::FreeMod{T}, d::Int;
    domain_strand::Tuple{<:FreeMod, <:Map, <:Map}=strand(F, d),
    codomain_strand::Tuple{<:FreeMod, <:Map, <:Map}=strand(F, d+1),
  ) where {T<:MPolyDecRingElem}
  @assert is_z_graded(F)
  Fd, Fd_to_F, F_to_Fd = domain_strand
  e = d+1
  Fe, Fe_to_F, F_to_Fe = codomain_strand
  S = base_ring(F)
  R = base_ring(S)
  E = _exterior_algebra(S)
  G = grading_group(E)
  Fd_E = graded_free_module(E, [-d*G[1] for i in 1:ngens(Fd)])
  Fe_E = graded_free_module(E, [-e*G[1] for i in 1:ngens(Fe)])
  Fe_map = hom(Fe, Fe_E, gens(Fe_E), E)
  Fd_map = hom(Fd, Fd_E, gens(Fd_E), E)
  img_gens = elem_type(Fe_E)[]
  for (i, g) in enumerate(gens(Fd))
    img = zero(Fe_E)
    for j in 1:nvars(S)
      h = g
      h = Fd_to_F(h)
      h = gen(S, j)*h
      h = F_to_Fe(h)
      h = Fe_map(h)
      h = gen(E, j)*h
      img += h
    end
    push!(img_gens, img)
  end
  return hom(Fd_E, Fe_E, img_gens), Fd_map, Fe_map
end

function _ext_module_map(M::SubquoModule{T}, d::Int) where {T<:MPolyDecRingElem}
  pres = presentation(M)
  F0 = pres[0]
  F1 = pres[1]
  b = map(pres, 1) # F1 -> F0
  F0_d, to_F0, to_F0_d = strand(F0, d)
  F1_d, to_F1, to_F1_d = strand(F1, d)
  e = d + 1
  F0_e, to_F0_2, to_F0_e = strand(F0, e)
  F1_e, to_F1_2, to_F1_e = strand(F1, e)

  dom_rels, inc_dom = sub(F0_d, to_F0_d.(b.(to_F1.(gens(F1_d)))))
  cod_rels, inc_cod = sub(F0_e, to_F0_e.(b.(to_F1_2.(gens(F1_e)))))
  phi, F0_d_map, F0_e_map = _ext_module_map(F0, d; domain_strand = (F0_d, to_F0, to_F0_d),
                                            codomain_strand = (F0_e, to_F0_2, to_F0_e))
  @assert domain(F0_d_map) === F0_d
  @assert domain(F0_e_map) === F0_e
  @assert domain(phi) === codomain(F0_d_map)
  @assert codomain(phi) === codomain(F0_e_map)
  
  M_dom, pr_dom = quo(domain(phi), F0_d_map.(ambient_representatives_generators(dom_rels)))
  M_cod, pr_cod = quo(codomain(phi), F0_e_map.(ambient_representatives_generators(cod_rels)))
  @assert domain(pr_cod) === codomain(phi)
  
  psi = hom(M_dom, M_cod, pr_cod.(phi.(F0_d_map.(gens(F0_d)))))

  return psi, MapFromFunc(M_dom, M, v->map(pres, 0)(to_F0(repres(v)))), MapFromFunc(M_cod, M, v->map(pres, 0)(to_F0_2(repres(v))))
end

function gen(E::ExteriorAlgebra{T}, i::Int) where {T}
  return ExtAlgElem(E, Dict{Int, SRow{T}}(1 => sparse_row(base_ring(E), [(i, one(base_ring(E)))]))) 
end
  


@attr ExteriorAlgebra{T} function _exterior_algebra(S::MPolyDecRing{T}) where {T}
  @assert is_z_graded(S)
  R = base_ring(S)
  new_symb = Symbol.([string(s)*" ̌" for s in symbols(S)])
  return ExteriorAlgebra(R, new_symb)
end

function _derived_pushforward_BGG(M::ModuleFP{T}; regularity::Int=_regularity_bound(M)) where {T <: MPolyDecRingElem{<:RingElem}}
  phi, _, _ = _ext_module_map(M, regularity+1)
  S = base_ring(M)
  R = base_ring(S)
  I, inc = kernel(phi)
  res, aug = free_resolution(SimpleFreeResolution, I)
  res0, red = change_base_ring(R, res)
  return res0
end


