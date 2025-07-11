###############################################################################
# SubquoModuleElem constructors
###############################################################################

@doc raw"""
    SubquoModuleElem(v::SRow{R}, SQ::SubquoModule) where {R}

Return the element $\sum_i v[i] \cdot SQ[i]$.
"""
SubquoModuleElem(v::SRow{R}, SQ::SubquoModule) where {R} = SubquoModuleElem{R}(v, SQ)

@doc raw"""
    SubquoModuleElem(a::FreeModElem{R}, SQ::SubquoModule) where {R}

Construct an element $v \in SQ$ that is represented by $a$.
"""
SubquoModuleElem(a::FreeModElem{R}, SQ::SubquoModule; is_reduced::Bool=false) where {R} = SubquoModuleElem{R}(a, SQ; is_reduced)

elem_type(::Type{SubquoModule{T}}) where {T} = SubquoModuleElem{T}
parent_type(::Type{SubquoModuleElem{T}}) where {T} = SubquoModule{T}

function in(v::SubquoModuleElem, M::SubquoModule)
  ambient_free_module(parent(v)) === ambient_free_module(M) || return false
  return represents_element(repres(v), M)
end

@doc raw"""
    getindex(v::SubquoModuleElem, i::Int)

Let $v \in M$ with $v = \sum_i a[i] \cdot M[i]$. Return $a[i]$
"""
function getindex(v::SubquoModuleElem, i::Int)
  if isempty(coordinates(v))
    return zero(base_ring(v.parent))
  end
  return coordinates(v)[i]
end

#######################################################
@doc raw"""
    coordinates(m::SubquoModuleElem)

Given an element `m` of a subquotient $M$ over a ring $R$, say,
return the coefficients of an $R$-linear combination of the generators of $M$
which gives $m$.

Return the coefficients of `m` with respect to the basis of standard unit vectors.

The result is returned as a sparse row.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; x*y; y^2; z^4]
[x^2]
[x*y]
[y^2]
[z^4]

julia> M = subquotient(A, B);

julia> m = z*M[1] + M[2]
(x*z + y)*e[1]

julia> coordinates(m)
Sparse row with positions [1, 2] and values QQMPolyRingElem[z, 1]
```
"""
function coordinates(m::SubquoModuleElem)
  if !isdefined(m, :coeffs)
    @assert isdefined(m, :repres) "neither coeffs nor repres is defined on a SubquoModuleElem"
    m.coeffs = coordinates(repres(m), parent(m))
  end
  return m.coeffs
end

#########################################################

@doc raw"""
    repres(v::SubquoModuleElem)

Return a free module element that is a representative of `v`.
"""
function repres(v::SubquoModuleElem)
  if !isdefined(v, :repres)
    @assert isdefined(v, :coeffs) "neither coeffs nor repres is defined on a SubquoModuleElem"
    M = parent(v)
    v.repres = sum(a*M.sub[i] for (i, a) in coordinates(v); init=zero(M.sub))
  end
  return v.repres
end

#######################################################

# simplify modifies the representative v of el as follows:
#
#  - if el is zero, v is zero
#  - if el is homogeneous, but the current representative is not
#    then a homogeneous representative is returned.
#  - it sets the field is_reduced to true.
function simplify(el::SubquoModuleElem{<:MPolyRingElem{T}}) where {T<:Union{<:FieldElem, <:ZZRingElem}}
  el.is_reduced && return el
  if is_zero(el) # We have to do this check because otherwise the coordinates of the representative are not reset.
    result = zero(parent(el))
    result.is_reduced = true # Todo: Should be done in zero(...)
    return result
  end
  if !isdefined(parent(el), :quo) || is_zero(parent(el).quo)
    el.is_reduced = true
    return el
  end
  !isdefined(el, :repres) && repres(el) # Make sure the field is filled
  reduced = reduce(el.repres, parent(el).quo)
  result = SubquoModuleElem(reduced, parent(el), is_reduced=true)
  return result
end

function simplify!(el::SubquoModuleElem{<:MPolyRingElem{T}}) where {T<:Union{<:FieldElem, <:ZZRingElem}}
  el.is_reduced && return el
  if is_zero(el) # We have to do this check because otherwise the coordinates of the representative are not reset.
    result = zero(parent(el))
    result.is_reduced = true # Todo: Should be done in zero(...)
    return result
  end
  if !isdefined(parent(el), :quo) || is_zero(parent(el).quo)
    el.is_reduced = true
    return el
  end
  !isdefined(el, :repres) && repres(el) # Make sure the field is filled
  el.repres = reduce(el.repres, parent(el).quo)
  el.is_reduced = true
  return el
end

# The default only checks whether an element is zero.
function simplify(el::SubquoModuleElem)
  el.is_reduced && return el
  if is_zero(el)
    result = zero(parent(el))
    result.is_reduced = true # Todo: Should be done in zero(...)
    return result
  end
  el.is_reduced = true
  return el
end

function simplify!(el::SubquoModuleElem)
  el.is_reduced && return el
  if is_zero(el)
    el.coeffs = sparse_row(base_ring(parent(el)))
    el.repres = zero(ambient_free_module(parent(el)))
    el.is_reduced = true
    return el
  end
  el.is_reduced = true
  return el
end


#######################################################
@doc raw"""
    ambient_representative(m::SubquoModuleElem)

Given an element `m` of a subquotient $M$, say, return

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; x*y; y^2; z^4]
[x^2]
[x*y]
[y^2]
[z^4]

julia> M = subquotient(A, B);

julia> m = z*M[1] + M[2]
(x*z + y)*e[1]

julia> typeof(m)
SubquoModuleElem{QQMPolyRingElem}

julia> fm = ambient_representative(m)
(x*z + y)*e[1]

julia> typeof(fm)
FreeModElem{QQMPolyRingElem}

julia> parent(fm) == ambient_free_module(M)
true
```
"""
ambient_representative(m::SubquoModuleElem) = repres(m)

# another method for compatibility in generic code
ambient_representative(a::FreeModElem) = a

#######################################################

@doc raw"""
    Vector(v::SubquoModuleElem)

Return the coefficients of a representative of `v` as a Vector.
"""
function Vector(v::SubquoModuleElem)
  return Vector(repres(v))
end

@doc raw"""
    standard_basis(F::ModuleGens{T}, reduced::Bool=false) where {T <: MPolyRingElem}

Return a standard basis of `F` as an object of type `ModuleGens`.
If `reduced` is set to `true` and the ordering of the underlying ring is global,
a reduced Gröbner basis is computed.
"""
function standard_basis(F::ModuleGens{T}, reduced::Bool=false) where {T <: MPolyRingElem}
  @req is_exact_type(elem_type(base_ring(F))) "This functionality is only supported over exact fields."
  if reduced
    @assert has_global_singular_ordering(F)
  end
  if singular_generators(F).isGB && !reduced
    return F
  end
  return ModuleGens(F.F, Singular.std(singular_generators(F), complete_reduction=reduced))
end

@doc raw"""
    lift_std(M::ModuleGens{T}) where {T <: MPolyRingElem}

Return a standard basis `G` of `F` as an object of type `ModuleGens` along with
a transformation matrix `T` such that `T*matrix(M) = matrix(G)`.
"""
function lift_std(M::ModuleGens{T}) where {T <: MPolyRingElem}
  R = base_ring(M)
  G,Trans_mat = Singular.lift_std(singular_generators(M)) # When Singular supports reduction add it also here
  mg = ModuleGens(M.F, G)
  mg.isGB = true
  mg.S.isGB = true
  mg.ordering = default_ordering(M.F)
  mat = map_entries(R, transpose(Trans_mat))
  set_attribute!(mg, :transformation_matrix => mat)
  return mg, mat
end

@doc raw"""
    lift_std(M::ModuleGens{T}, ordering::ModuleOrdering) where {T <: MPolyRingElem}

Return a standard basis `G` of `F` with respect to the given `ordering`
as an object of type `ModuleGens` along with a transformation
matrix `T` such that `T*matrix(M) = matrix(G)`.
"""
function lift_std(M::ModuleGens{T}, ordering::ModuleOrdering) where {T <: MPolyRingElem}
  M = ModuleGens(oscar_generators(M), M.F, ordering)
  mg, mat = lift_std(M)
  mg.ordering = ordering
  return mg, mat
end

@doc raw"""
    leading_monomials(F::ModuleGens)

Return the leading monomials of `F` as an object of type `ModuleGens`.
The leading module is with respect to the ordering defined on the Singular side.
"""
function leading_monomials(F::ModuleGens)
  # TODO
  # The following doesn't work yet. When comparison / lead for module elements
  # is implemented this should be uncommented.
  #if !isdefined(F, :S)
    #return ModuleGens(F.F, [lead(g) for g in oscar_generators(F)])
  #end
  singular_gens = singular_generators(F)
  return ModuleGens(oscar_free_module(F), Singular.lead(singular_gens))
end

function show(io::IO, b::SubquoModuleElem)
  print(io, repres(b))
end

@doc raw"""
    parent(b::SubquoModuleElem)

Let $b \in M$. Return $M$.
"""
parent(b::SubquoModuleElem) = b.parent

@doc raw"""
    (M::SubquoModule{T})(f::FreeModElem{T}) where T

Given an element `f` of the ambient free module of `M` which represents an element of `M`,
return the represented element.
"""
function (M::SubquoModule{T})(f::FreeModElem{T}) where T
  coords = coordinates(f, M)
  if coords === nothing
    error("not in the module")
  end
  return SubquoModuleElem(coords, M)
end

@doc raw"""
    (M::SubquoModule{T})(c::SRow{T}) where T

Return the subquotient element $\sum_i a[i] \cdot M[i]\in M$.
"""
function (R::SubquoModule)(a::SRow)
  return SubquoModuleElem(a, R)
end

@doc raw"""
    SubquoModuleElem(c::Vector{T}, parent::SubquoModule{T}) where T

Return the element of  `parent`  defined as a linear combination
of the generators of $parent$ with coefficients given by the entries of `c`.
"""
function SubquoModuleElem(c::Vector{T}, parent::SubquoModule{T}) where T
  @assert length(c) == ngens(parent)
  sparse_coords = sparse_row(base_ring(parent), collect(1:ngens(parent)), c)
  return SubquoModuleElem{T}(sparse_coords,parent)
end


@doc raw"""
    (M::SubquoModule{T})(c::Vector{T}) where T

Return the element of `M` defined as a linear combination
of the generators of $M$ with coefficients given by the entries of `c`.
"""
function (M::SubquoModule{T})(c::Vector{T}) where T
 return SubquoModuleElem(c, M)
end

@doc raw"""
    (R::SubquoModule)(a::SubquoModuleElem)

Return `a` if it lives in `R`.
"""
function (R::SubquoModule)(a::SubquoModuleElem)
  if parent(a) == R
    return a
  end
  error("illegal coercion")
end

@doc raw"""
    index_of_gen(v::SubquoModuleElem)

Let $v \in G$ with $v$ the `i`th generator of $G$. Return `i`.
"""
function index_of_gen(v::SubquoModuleElem)
  @assert length(coordinates(v).pos) == 1
  @assert isone(coordinates(v).values[1])
  return coordinates(v).pos[1]
end

# function to check whether two module elements are in the same module
function check_parent(a::Union{AbstractFreeModElem,SubquoModuleElem}, b::Union{AbstractFreeModElem,SubquoModuleElem})
  if parent(a) !== parent(b)
    error("elements not compatible")
  end
end

function +(a::SubquoModuleElem, b::SubquoModuleElem)
  check_parent(a,b)
  if isdefined(a, :coeffs) && isdefined(b, :coeffs)
    return SubquoModuleElem(coordinates(a)+coordinates(b), a.parent)
  else
    return SubquoModuleElem(repres(a) + repres(b), parent(a))
  end
end

function -(a::SubquoModuleElem, b::SubquoModuleElem)
  check_parent(a,b)
  if isdefined(a, :coeffs) && isdefined(b, :coeffs)
    return SubquoModuleElem(coordinates(a)-coordinates(b), a.parent)
  else
    return SubquoModuleElem(repres(a) - repres(b), parent(a))
  end
end

-(a::SubquoModuleElem) = SubquoModuleElem(-coordinates(a), a.parent)

function *(a::MPolyDecRingElem, b::SubquoModuleElem)
  if parent(a) !== base_ring(parent(b))
    return base_ring(parent(b))(a)*b # this will throw if conversion is not possible
  end
  isdefined(b, :coeffs) && return SubquoModuleElem(a*coordinates(b), b.parent)
  return SubquoModuleElem(a*repres(b), b.parent)
end

function *(a::MPolyRingElem, b::SubquoModuleElem)
  if parent(a) !== base_ring(parent(b))
    return base_ring(parent(b))(a)*b # this will throw if conversion is not possible
  end
  isdefined(b, :coeffs) && return SubquoModuleElem(a*coordinates(b), b.parent)
  return SubquoModuleElem(a*repres(b), b.parent)
end

function *(a::NCRingElem, b::SubquoModuleElem)
  if parent(a) !== base_ring(parent(b))
    return base_ring(parent(b))(a)*b # this will throw if conversion is not possible
  end
  isdefined(b, :coeffs) && return SubquoModuleElem(a*coordinates(b), b.parent)
  return SubquoModuleElem(a*repres(b), b.parent)
end

*(a::Int, b::SubquoModuleElem) = SubquoModuleElem(a*coordinates(b), b.parent)
*(a::Integer, b::SubquoModuleElem) = SubquoModuleElem(a*coordinates(b), b.parent)
*(a::QQFieldElem, b::SubquoModuleElem) = SubquoModuleElem(a*coordinates(b), b.parent)

function (==)(a::SubquoModuleElem, b::SubquoModuleElem)
  if parent(a) !== parent(b)
    return false
  end
  return iszero(a-b)
end

function Base.hash(a::SubquoModuleElem)
  b = 0xaa2ba4a32dd0b431 % UInt
  h = hash(typeof(a), h)
  h = hash(parent(a), h)
  return xor(h, b)
end

function Base.hash(a::SubquoModuleElem{<:MPolyRingElem{<:FieldElem}}, h::UInt)
  b = 0xaa2ba4a32dd0b431 % UInt
  h = hash(typeof(a), h)
  h = hash(parent(a), h)
  simplify!(a)
  return hash(a.repres, h)
end

function Base.deepcopy_internal(a::SubquoModuleElem, dict::IdDict)
  return SubquoModuleElem(deepcopy_internal(coordinates(a), dict), a.parent)
end

@doc raw"""
    sub(F::FreeMod{T}, V::Vector{<:FreeModElem{T}}; cache_morphism::Bool=false) where {T}

Given a vector `V` of (homogeneous) elements of `F`, return a pair `(I, inc)`
consisting of the (graded) submodule `I` of `F` generated by these elements
and its inclusion map `inc : I ↪ F`.

When `cache_morphism` is set to true, then `inc` will be cached and available
for `transport` and friends.

If only the submodule itself is desired, use `sub_object` instead.
"""
function sub(F::FreeMod{T}, V::Vector{<:FreeModElem{T}}; cache_morphism::Bool=false) where {T}
  s = SubquoModule(F, V)
  emb = hom(s, F, V; check=false)
  set_attribute!(s, :canonical_inclusion => emb)
  cache_morphism && register_morphism!(emb)
  return s, emb
end

function sub_object(F::FreeMod{T}, V::Vector{<:FreeModElem{T}}) where {T}
  return SubquoModule(F, V)
end

@doc raw"""
    sub(F::FreeMod{T}, A::MatElem{T}; cache_morphism::Bool=false) where {T}

Given a (homogeneous) matrix `A` interpret the rows of `A` as elements
of the free module `F` and return a pair `(I, inc)`
consisting of the (graded) submodule `I` of `F` generated by these row vectors,
together with its inclusion map `inc : I ↪ F`.

When `cache_morphism` is set to true, then `inc` will be cached and available
for `transport` and friends.

If only the submodule itself is desired, use `sub_object` instead.
"""
function sub(F::FreeMod{T}, A::MatElem{T}; cache_morphism::Bool=false) where {T}
  M = SubquoModule(SubModuleOfFreeModule(F, A))
  #M = SubquoModule(F, A, zero_matrix(base_ring(F), 1, rank(F)))
  emb = hom(M, F, ambient_representatives_generators(M); check=false)
  emb.matrix = A
  set_attribute!(M, :canonical_inclusion => emb) # TODO: Can this be removed?
  cache_morphism && register_morphism!(emb)
  return M, emb
end

function sub_object(F::FreeMod{T}, A::MatElem{T}) where {T}
  return SubquoModule(SubModuleOfFreeModule(F, A))
end

@doc raw"""
    sub(F::FreeMod{T}, O::Vector{<:SubquoModuleElem{T}}; cache_morphism::Bool=false) where T

Suppose the `ambient_free_module` of the `parent` `M` of the elements `v_i`
in `O` is `F` and `M` is a submodule (i.e. no relations are present).
Then this returns a pair `(I, inc)` consisting of the submodule `I`
generated by the elements in `O` in `F`, together with its inclusion
morphism `inc : I ↪ F`.

When `cache_morphism` is set to true, then `inc` will be cached and available
for `transport` and friends.

If only the submodule itself is desired, use `sub_object` instead.
"""
function sub(F::FreeMod{T}, O::Vector{<:SubquoModuleElem{T}}; cache_morphism::Bool=false) where T
  s = SubquoModule(F, [repres(x) for x = O])
  return sub(F, s; cache_morphism)
end

function sub_object(F::FreeMod{T}, O::Vector{<:SubquoModuleElem{T}}) where T
  return SubquoModule(F, [repres(x) for x = O])
end

@doc raw"""
    sub(F::FreeMod{T}, M::SubquoModule{T}; cache_morphism::Bool=false) where T

Return `M` as a submodule of `F`, together with its inclusion morphism
`inc : M ↪ F`.

When `cache_morphism` is set to true, then `inc` will be cached and available
for `transport` and friends.

The `ambient_free_module` of `M` needs to be `F` and `M` has to have no
relations.

If only the submodule itself is desired, use `sub_object` instead.
"""
function sub(F::FreeMod{T}, s::SubquoModule{T}; cache_morphism::Bool=false) where T
  @assert !isdefined(s, :quo)
  @assert s.F === F
  emb = hom(s, F, elem_type(F)[repres(x) for x in gens(s)]; check=false)
  set_attribute!(s, :canonical_inclusion => emb)
  cache_morphism && register_morphism!(emb)
  return s, emb
end

function sub_object(F::FreeMod{T}, s::SubquoModule{T}) where T
  @assert !isdefined(s, :quo)
  @assert s.F === F
  return s
end

@doc raw"""
    sub(M::SubquoModule{T}, V::Vector{<:SubquoModuleElem{T}}; cache_morphism::Bool=false) where T

Given a vector `V` of (homogeneous) elements of `M`, return the (graded) submodule `I` of `M` generated by these elements
together with its inclusion map `inc : I ↪ M.

When `cache_morphism` is set to true, then `inc` will be cached and available
for `transport` and friends.

If only the submodule itself is desired, use `sub_object` instead.
"""
function sub(M::SubquoModule{T}, V::Vector{<:SubquoModuleElem{T}}; cache_morphism::Bool=false) where T
  @assert all(x -> x.parent === M, V)
  t = SubquoModule(M.F, FreeModElem[repres(x) for x in V])
  if isdefined(M, :quo)
    t.quo = M.quo
    t.sum = sum(t.sub, t.quo)
  end
  emb = hom(t, M, V; check=false)
  set_attribute!(t, :canonical_inclusion => emb) # TODO: Can this be removed?
  cache_morphism && register_morphism!(emb)
  return t, emb
end

function sub_object(M::SubquoModule{T}, V::Vector{<:SubquoModuleElem{T}}) where T
  @assert all(x -> x.parent === M, V)
  t = SubquoModule(M.F, FreeModElem[repres(x) for x in V])
  if isdefined(M, :quo)
    t.quo = M.quo
    t.sum = sum(t.sub, t.quo)
  end
  return t
end

@doc raw"""
    sub(M::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}; cache_morphism::Bool=false) where T

Given a vector `V` of (homogeneous) elements of `M`, return the (graded) submodule `I` of `M` generated by these elements
together with its inclusion map `inc : I ↪ M.

When `cache_morphism` is set to true, then `inc` will be cached and available
for `transport` and friends.

If only the submodule itself is desired, use `sub_object` instead.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> F = free_module(R, 1);

julia> V = [x^2*F[1]; y^3*F[1]; z^4*F[1]];

julia> N, incl = sub(F, V);

julia> N
Submodule with 3 generators
  1: x^2*e[1]
  2: y^3*e[1]
  3: z^4*e[1]
represented as subquotient with no relations

julia> incl
Module homomorphism
  from N
  to F
```
"""
function sub(M::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}; cache_morphism::Bool=false) where T
 error("sub is not implemented for the given types.")
end

#=
@doc raw"""
    return_sub_wrt_task(M::SubquoModule, emb::SubQuoHom, task::Symbol)

This helper function returns `M`, `emb` or both according
to `task`.
"""
function return_sub_wrt_task(M::SubquoModule, emb::SubQuoHom, task::Symbol)
  (task == :none || task == :module) && return M
  task == :cache_morphism && register_morphism!(emb)
  task == :only_morphism && return emb
  (task == :cache_morphism || task == :both || task == :with_morphism) && return M, emb
  error("No valid option for task.")
end
=#

@doc raw"""
    quo(F::FreeMod{T}, V::Vector{<:FreeModElem{T}}; cache_morphism::Bool=false) where T

Given a vector `V` of (homogeneous) elements of `F`, return a pair `(M, pr)` consisting
of the quotient `M = F/⟨V⟩` and the projection map `pr : F → M`.

If one is only interested in the actual object `M`, but not the map, use `quo_object` instead.

If `cache_morphism` is set to `true`, the projection is cached and available to `transport` and friends.
"""
function quo(F::FreeMod{T}, V::Vector{<:FreeModElem{T}}; cache_morphism::Bool=false) where T
  S = SubquoModule(F, basis(F))
  Q = SubquoModule(S, V)
  phi = hom(F, Q, gens(Q), check=false)
  cache_morphism && register_morphism!(phi)
  return Q, phi
end

function quo_object(F::FreeMod{T}, V::Vector{<:FreeModElem{T}}) where T
  S = SubquoModule(F, basis(F))
  Q = SubquoModule(S, V)
  return Q
end

@doc raw"""
    quo(F::FreeMod{T}, A::MatElem{T}; cache_morphism::Bool=false) where {T}

Given a matrix `A`, interpret the row vectors `v_i` of `A` as elements of `F`
and return a pair `(M, pr)` consisting of the quotient `M = F/I` of `F` by the
submodule `I ⊂ F` generated by the rows of `A`, together with the projection
map `pr : F → M`.

If one is only interested in the actual object `M`, but not the map, use `quo_object` instead.

If `cache_morphism` is set to `true`, the projection is cached and available to `transport` and friends.
"""
function quo(F::FreeMod{T}, A::MatElem{T}; cache_morphism::Bool=false) where {T}
  E = identity_matrix(base_ring(F), rank(F))
  Q = SubquoModule(F, E, A)

  phi = hom(F, Q, gens(Q), check=false)
  cache_morphism && register_morphism!(phi)
  return Q, phi
end

function quo_object(F::FreeMod{T}, A::MatElem{T}) where {T}
  E = identity_matrix(base_ring(F), rank(F))
  Q = SubquoModule(F, E, A)
  return Q
end

@doc raw"""
    quo(F::FreeMod{T}, O::Vector{<:SubquoModuleElem{T}}; cache_morphism::Bool=false) where T

Given a vector `O` of (homogeneous) elements of some submodule `I` of `F`,
return a pair `(M, pr)` consisting of the quotient `M = F/⟨V⟩` and
the projection map `pr : F → M`.

If one is only interested in the actual object `M`, but not the map, use `quo_object` instead.
Compute $F / T$, where $T$ is generated by $O$.

Note that the submodule `I` must have `F` as its `ambient_free_module`.

If `cache_morphism` is set to `true`, the projection is cached and available to `transport` and friends.
"""
function quo(F::FreeMod{T}, O::Vector{<:SubquoModuleElem{T}}; cache_morphism::Bool=false) where T
  S = SubquoModule(F, basis(F))
  Q = SubquoModule(S, [repres(x) for x = O])

  phi = hom(F, Q, gens(Q), check=false)
  cache_morphism && register_morphism!(phi)
  return Q, phi
end

function quo_object(F::FreeMod{T}, O::Vector{<:SubquoModuleElem{T}}) where T
  S = SubquoModule(F, basis(F))
  Q = SubquoModule(S, [repres(x) for x = O])
  return Q
end

@doc raw"""
    quo(S::SubquoModule{T}, O::Vector{<:FreeModElem{T}}; cache_morphism::Bool=false) where T

Given a vector `V` of (homogeneous) elements of `S`, return a pair `(M, pr)` consisting
of the quotient `M = S/⟨V⟩` and the projection map `pr : S → M`.

If one is only interested in the actual object `M`, but not the map, use `quo_object` instead.

Note that the elements of `O` must belong to the `ambient_free_module` of `S` and represent
elements of `S`.

If `cache_morphism` is set to `true`, the projection is cached and available to `transport` and friends.
"""
function quo(F::SubquoModule{T}, O::Vector{<:FreeModElem{T}}; cache_morphism::Bool=false) where T
  if length(O) > 0
    @assert parent(O[1]) === F.F
  end
  if isdefined(F, :quo)
    SF = singular_freemodule(F.quo.gens)
    s = Singular.Module(base_ring(SF), [SF(x) for x in [O; oscar_generators(F.quo.gens)]]...)
    Q = SubquoModule(F.F, singular_generators(F.sub.gens), s)
    phi = hom(F, Q, gens(Q), check=false)
    cache_morphism && register_morphism!(phi)
    return Q, phi
  end
  Q = SubquoModule(F, O)
  phi = hom(F, Q, gens(Q), check=false)
  cache_morphism && register_morphism!(phi)
  return Q, phi
end

function quo_object(F::SubquoModule{T}, O::Vector{<:FreeModElem{T}}) where T
  if length(O) > 0
    @assert parent(O[1]) === F.F
  end
  if isdefined(F, :quo)
    SF = singular_freemodule(F.quo.gens)
    s = Singular.Module(base_ring(SF), [SF(x) for x in [O; oscar_generators(F.quo.gens)]]...)
    return SubquoModule(F.F, singular_generators(F.sub.gens), s)
  end
  return SubquoModule(F, O)
end

@doc raw"""
    quo(S::SubquoModule{T}, V::Vector{<:SubquoModuleElem{T}}; cache_morphism::Bool=false) where T

Given a vector `V` of (homogeneous) elements of `S`, return a pair `(M, pr)` consisting
of the quotient `M = S/⟨V⟩` and the projection map `pr : S → M`.

If one is only interested in the actual object `M`, but not the map, use `quo_object` instead.

If `cache_morphism` is set to `true`, the projection is cached and available to `transport` and friends.
"""
function quo(M::SubquoModule{T}, V::Vector{<:SubquoModuleElem{T}}; cache_morphism::Bool=false) where T
  return quo(M, [repres(x) for x = V]; cache_morphism)
end

function quo_object(M::SubquoModule{T}, V::Vector{<:SubquoModuleElem{T}}) where T
  return quo_object(M, [repres(x) for x = V])
end

@doc raw"""
    quo(M::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}; cache_morphism::Bool=false) where T

Given a vector `V` of (homogeneous) elements of `M`, return a pair `(N, pr)` consisting
of the quotient `N = M/⟨V⟩` and the projection map `pr : M → N`.

If one is only interested in the actual object `M`, but not the map, use `quo_object` instead.

If `cache_morphism` is set to `true`, the projection is cached and available to `transport` and friends.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> F = free_module(R, 1);

julia> V = [x^2*F[1]; y^3*F[1]; z^4*F[1]];

julia> N, proj = quo(F, V);

julia> N
Subquotient of submodule with 1 generator
  1: e[1]
by submodule with 3 generators
  1: x^2*e[1]
  2: y^3*e[1]
  3: z^4*e[1]

julia> proj
Module homomorphism
  from F
  to N
```
"""
function quo(M::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}; cache_morphism::Bool=false) where T
 error("quo is not implemented for the given types.")
end


@doc raw"""
    quo(M::SubquoModule{T}, U::SubquoModule{T}) where T

Return a pair `(N, pr)` consisting of the quotient $N = M / U$ together with the projection
map `pr : M → N`.

If one is only interested in the actual object `N`, but not the map, use `quo_object` instead.

If `cache_morphism` is set to `true`, the projection is cached and available to `transport` and friends.
"""
function quo(M::SubquoModule{T}, U::SubquoModule{T}; cache_morphism::Bool=false) where T
  if isdefined(M, :quo) && isdefined(U, :quo)
    F = ambient_free_module(M)
    @assert F === ambient_free_module(U)
    # We can not assume that the SubModuleOfFreeModule layer is implemented in general,
    # so we deflect to the Subquo-layer instead.
    @assert SubquoModule(F, relations(M)) == SubquoModule(F, relations(U))
  else
    @assert !isdefined(M, :quo) && !isdefined(U, :quo)
  end
  Q = SubquoModule(M, gens(U.sub))
  pr = hom(M, Q, gens(Q), check=false)
  cache_morphism && register_morphism!(pr)
  return Q, pr
end

function quo_object(M::SubquoModule{T}, U::SubquoModule{T}) where T
  if isdefined(M, :quo) && isdefined(U, :quo)
    F = ambient_free_module(M)
    @assert F === ambient_free_module(U)
    # We can not assume that the SubModuleOfFreeModule layer is implemented in general,
    # so we deflect to the Subquo-layer instead.
    @assert SubquoModule(F, relations(M)) == SubquoModule(F, relations(U))
  else
    @assert !isdefined(M, :quo) && !isdefined(U, :quo)
  end
  return SubquoModule(M, gens(U.sub))
end

function quo(M::SubquoModule{T}, U::SubquoModule{T}; cache_morphism::Bool=false) where {T<:MPolyRingElem}
  if isdefined(M, :quo) && isdefined(U, :quo)
    @assert M.quo == U.quo
  else
    @assert !isdefined(M, :quo) && !isdefined(U, :quo)
  end
  Q = SubquoModule(M, oscar_generators(U.sub.gens))
  pr = hom(M, Q, gens(Q), check=false)
  cache_morphism && register_morphism!(pr)
  return Q, pr
end

function quo_object(M::SubquoModule{T}, U::SubquoModule{T}) where {T<:MPolyRingElem}
  if isdefined(M, :quo) && isdefined(U, :quo)
    @assert M.quo == U.quo
  else
    @assert !isdefined(M, :quo) && !isdefined(U, :quo)
  end
  return SubquoModule(M, oscar_generators(U.sub.gens))
end

@doc raw"""
    quo(F::FreeMod{R}, T::SubquoModule{R}; cache_morphism::Bool=false) where R

Return a pair `(N, pr)` consisting of the quotient $N = F / T$ together with the projection
map `pr : F → N`.

If one is only interested in the actual object `N`, but not the map, use `quo_object` instead.

If `cache_morphism` is set to `true`, the projection is cached and available to `transport` and friends.
"""
function quo(F::FreeMod{R}, T::SubquoModule{R}; cache_morphism::Bool=false) where R
  @assert !isdefined(T, :quo)
  return quo(F, gens(T); cache_morphism)
end

function quo_object(F::FreeMod{R}, T::SubquoModule{R}) where R
  @assert !isdefined(T, :quo)
  return quo_object(F, gens(T))
end

@doc raw"""
    syzygy_module(F::ModuleGens; sub = FreeMod(base_ring(F.F), length(oscar_generators(F))))
"""
function syzygy_module(F::ModuleGens{T}; sub = FreeMod(base_ring(F.F), length(oscar_generators(F)))) where {T <: MPolyRingElem}
  # TODO Obtain the Gröbner basis and cache it
  s = Singular.syz(singular_generators(F))
  return SubquoModule(sub, s)
end

@doc raw"""
    gens(M::SubquoModule{T}) where T

Return the generators of `M`.
"""
function gens(M::SubquoModule{T}) where T
  R = base_ring(M)
  e = R(1)
  return [SubquoModuleElem{T}(sparse_row(R, [i], [e]), M) for i in 1:ngens(M)]
end

@doc raw"""
    gen(M::SubquoModule{T}, i::Int) where T

Return the `i`th generator of `M`.
"""
function gen(M::SubquoModule{T}, i::Int) where T
  R = base_ring(M)
  v = sparse_row(R, [i], [R(1)])
  return SubquoModuleElem{T}(v, M)
end

@doc raw"""
    number_of_generators(M::SubquoModule)

Return the number of generators of `M`.
"""
number_of_generators(M::SubquoModule) = number_of_generators(M.sub)

@doc raw"""
    base_ring(M::SubquoModule)

Given an `R`-module `M`, return `R`.
"""
base_ring(M::SubquoModule) = base_ring(M.F)::base_ring_type(M)

base_ring_type(::Type{SubquoModule{T}}) where {T} = base_ring_type(FreeMod{T})

@doc raw"""
    zero(M::SubquoModule)

Return the zero element of `M`.
"""
zero(M::SubquoModule) = SubquoModuleElem(sparse_row(base_ring(M)), M)

@doc raw"""
    is_zero(M::SubquoModule)

Return `true` if `M` is the zero module, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> F = free_module(R, 1)
Free module of rank 1 over R

julia> A = R[x^2+y^2;]
[x^2 + y^2]

julia> B = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = subquotient(F, A, B)
Subquotient of submodule with 1 generator
  1: (x^2 + y^2)*e[1]
by submodule with 3 generators
  1: x^2*e[1]
  2: y^3*e[1]
  3: z^4*e[1]

julia> is_zero(M)
false
```
"""
@attr Bool function is_zero(M::SubquoModule)
  return all(iszero, gens(M))
end

function iterate(F::ModuleGens, i::Int = 1)
  if i>length(F)
    return nothing
  else
    return F[i], i+1
  end
end
Base.eltype(::Type{ModuleGens{T}}) where {T} = FreeModElem{T}

#??? A scalar product....
function *(a::FreeModElem, b::Vector{FreeModElem})
  @assert rank(parent(a)) == length(b)
  s = zero(parent(a))
  for (p,v) in coordinates(a)
    s += v*b[p]
  end
  return s
end

@doc raw"""
    is_zero(m::SubquoModuleElem)

Return `true` if `m` is zero, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> F = free_module(R, 1)
Free module of rank 1 over R

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = subquotient(F, A, B)
Subquotient of submodule with 2 generators
  1: x*e[1]
  2: y*e[1]
by submodule with 3 generators
  1: x^2*e[1]
  2: y^3*e[1]
  3: z^4*e[1]

julia> is_zero(M[1])
false

julia> is_zero(x*M[1])
true
```

```jldoctest
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> F = graded_free_module(Rg, 1)
Graded free module Rg^1([0]) of rank 1 over Rg

julia> A = Rg[x; y]
[x]
[y]

julia> B = Rg[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = subquotient(F, A, B)
Graded subquotient of graded submodule of F with 2 generators
  1: x*e[1]
  2: y*e[1]
by graded submodule of F with 3 generators
  1: x^2*e[1]
  2: y^3*e[1]
  3: z^4*e[1]

julia> is_zero(M[1])
false

julia> is_zero(x*M[1])
true

```
"""
function is_zero(m::SubquoModuleElem)
  is_zero(ambient_representative(m)) && return true
  m.is_reduced && return false
  isdefined(parent(m), :quo) || return false
  return (ambient_representative(m) in parent(m).quo)
end

function is_zero(m::SubquoModuleElem{<:MPolyRingElem{T}}) where {T<:Union{ZZRingElem, <:FieldElem}}
  C = parent(m)
  isdefined(m, :coeffs) && is_zero(m.coeffs) && return true
  is_zero(repres(m)) && return true
  if !isdefined(C, :quo)
    return false
  end
  x = reduce(repres(m), C.quo)
  return iszero(x)
end


