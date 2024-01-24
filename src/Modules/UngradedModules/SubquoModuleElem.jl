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
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; x*y; y^2; z^4]
[x^2]
[x*y]
[y^2]
[z^4]

julia> M = SubquoModule(A, B);

julia> m = z*M[1] + M[2]
(x*z + y)*e[1]

julia> coordinates(m)
Sparse row with positions [1, 2] and values QQMPolyRingElem[z, 1]
```
"""
function coordinates(m::SubquoModuleElem)
  if !isdefined(m, :coeffs)
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
    M = parent(v)
    v.repres = sum(a*M.sub[i] for (i, a) in coordinates(v); init=zero(M.sub))
  end
  return v.repres
end

#######################################################

function simplify(el::SubquoModuleElem{<:MPolyRingElem{<:FieldElem}})
  el.is_reduced && return el
  !isdefined(parent(el), :quo) && return el
  iszero(parent(el).quo) && return el
  !isdefined(el, :repres) && repres(el) # Make sure the field is filled
  reduced = reduce(el.repres, parent(el).quo)
  result = SubquoModuleElem(reduced, parent(el), is_reduced=true)
  return result
end

function simplify!(el::SubquoModuleElem{<:MPolyRingElem{<:FieldElem}})
  el.is_reduced && return el
  !isdefined(parent(el), :quo) && return el
  iszero(parent(el).quo) && return el
  !isdefined(el, :repres) && repres(el) # Make sure the field is filled
  el.repres = reduce(el.repres, parent(el).quo)
  el.is_reduced = true
  return el
end

#######################################################
@doc raw"""
    ambient_representative(m::SubquoModuleElem)

Given an element `m` of a subquotient $M$, say, return 

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; x*y; y^2; z^4]
[x^2]
[x*y]
[y^2]
[z^4]

julia> M = SubquoModule(A, B);

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
  singular_assure(F)
  if reduced
    @assert Singular.has_global_ordering(base_ring(F.SF))
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
  singular_assure(M)
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
  M = ModuleGens(M.O, M.F, ordering)
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
    #return ModuleGens(F.F, [lead(g) for g in F.O])
  #end
  singular_assure(F)
  singular_gens = singular_generators(F)
  return ModuleGens(F.F, Singular.lead(singular_gens))
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

function *(a::RingElem, b::SubquoModuleElem) 
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
    sub(F::FreeMod{T}, V::Vector{<:FreeModElem{T}}, task::Symbol = :with_morphism) where T

Given a vector `V` of (homogeneous) elements of `F`, return the (graded) submodule of `F` generated by these elements.

Put more precisely, return the submodule as an object of type `SubquoModule`. 

Additionally, if `N` denotes this object, 

- return the inclusion map `N` $\to$ `F` if `task = :with_morphism` (default),
- return and cache the inclusion map `N` $\to$ `F` if `task = :cache_morphism`,
- do none of the above if `task = :none`.

If `task = :only_morphism`, return only the inclusion map.
"""
function sub(F::FreeMod{T}, V::Vector{<:FreeModElem{T}}, task::Symbol = :with_morphism) where T
  s = SubquoModule(F, V)
  emb = hom(s, F, V; check=false)
  set_attribute!(s, :canonical_inclusion => emb)
  return return_sub_wrt_task(s, emb, task)
end

@doc raw"""
    sub(F::FreeMod{T}, A::MatElem{T}, task::Symbol = :with_morphism) where {T} 

Given a (homogeneous) matrix `A`, return the (graded) submodule of `F` generated by the rows of `A`.

Put more precisely, return this submodule as an object of type `SubquoModule`. 

Additionally, if `N` denotes this submodule, 

- return the inclusion map `N` $\to$ `F` if `task = :with_morphism` (default),
- return and cache the inclusion map `N` $\to$ `F` if `task = :cache_morphism`,
- do none of the above if `task = :none`.

If `task = :only_morphism`, return only the inclusion map.
"""
function sub(F::FreeMod{T}, A::MatElem{T}, task::Symbol = :with_morphism) where {T}
  M = SubquoModule(SubModuleOfFreeModule(F, A)) 
  #M = SubquoModule(F, A, zero_matrix(base_ring(F), 1, rank(F)))
  emb = hom(M, F, ambient_representatives_generators(M); check=false)
  emb.matrix = A
  set_attribute!(M, :canonical_inclusion => emb)
  return return_sub_wrt_task(M, emb, task)
end

@doc raw"""
    sub(F::FreeMod{T}, O::Vector{<:SubquoModuleElem{T}}, task::Symbol = :with_morphism) where T

Return `S` as a submodule of `F`, where `S` is generated by `O`.
The embedding module of the parent of the elements of `O` must be `F`.
If `task` is set to `:none` or to `:module` return only `S`.
If `task` is set to `:with_morphism` (default option) or to `:both` return also the canonical injection morphism
$S \to F$.
If `task` is set to `:cache_morphism` the morphism is also cached.
If `task` is set to `:only_morphism` return only the morphism.
"""
function sub(F::FreeMod{T}, O::Vector{<:SubquoModuleElem{T}}, task::Symbol = :with_morphism) where T
  s = SubquoModule(F, [repres(x) for x = O])
  return sub(F, s, task)
end

@doc raw"""
    sub(F::FreeMod{T}, s::SubquoModule{T}, task::Symbol = :with_morphism) where T

Return `s` as a submodule of `F`, that is the embedding free module of `s` must 
be `F` and `s` has no relations.
If `task` is set to `:none` or to `:module` return only `s`.
If `task` is set to `:with_morphism` (default option) or to `:both` return also the canonical injection morphism
$s \to F$.
If `task` is set to `:cache_morphism` the morphism is also cached.
If `task` is set to `:only_morphism` return only the morphism.
"""
function sub(F::FreeMod{T}, s::SubquoModule{T}, task::Symbol = :with_morphism) where T
  @assert !isdefined(s, :quo)
  @assert s.F === F
  emb = hom(s, F, elem_type(F)[repres(x) for x in gens(s)]; check=false)
  #emb = hom(s, F, [FreeModElem(x.repres.coords, F) for x in gens(s)])
  set_attribute!(s, :canonical_inclusion => emb)
  return return_sub_wrt_task(s, emb, task)
end

@doc raw"""
    sub(M::SubquoModule{T}, V::Vector{<:SubquoModuleElem{T}}, task::Symbol = :with_morphism) where T

Given a vector `V` of (homogeneous) elements of `M`, return the (graded) submodule of `M` generated by these elements.

Put more precisely, return this submodule as an object of type `SubquoModule`. 

Additionally, if `N` denotes this object,

- return the inclusion map `N` $\to$ `M` if `task = :with_morphism` (default),
- return and cache the inclusion map `N` $\to$ `M` if `task = :cache_morphism`,
- do none of the above if `task = :none`.

If `task = :only_morphism`, return only the inclusion map.
"""
function sub(M::SubquoModule{T}, V::Vector{<:SubquoModuleElem{T}}, task::Symbol = :with_morphism) where T
  @assert all(x -> x.parent === M, V)
  t = SubquoModule(M.F, FreeModElem[repres(x) for x in V])
  if isdefined(M, :quo)
    t.quo = M.quo
    t.sum = sum(t.sub, t.quo)
  end
  emb = hom(t, M, V; check=false)
  set_attribute!(t, :canonical_inclusion => emb)
  return return_sub_wrt_task(t, emb, task)
end

@doc raw"""
    sub(M::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}, task::Symbol = :with_morphism) where T

Given a vector `V` of (homogeneous) elements of `M`, return the (graded) submodule of `M` generated by these elements.

Put more precisely, return this submodule as an object of type `SubquoModule`. 

Additionally, if `N` denotes this object,

- return the inclusion map `N` $\to$ `M` if `task = :with_morphism` (default),
- return and cache the inclusion map `N` $\to$ `M` if `task = :cache_morphism`,
- do none of the above if `task = :none`.

If `task = :only_morphism`, return only the inclusion map.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

julia> F = free_module(R, 1);

julia> V = [x^2*F[1]; y^3*F[1]; z^4*F[1]];

julia> N, incl = sub(F, V);

julia> N
Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]
represented as subquotient with no relations.

julia> incl
Map with following data
Domain:
=======
Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]
represented as subquotient with no relations.
Codomain:
=========
Free module of rank 1 over Multivariate polynomial ring in 3 variables over QQ
```
"""
function sub(M::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}, task::Symbol = :with_morphism) where T
 error("sub is not implemented for the given types.")
end

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

@doc raw"""
    quo(F::FreeMod{T}, V::Vector{<:FreeModElem{T}}, task::Symbol = :with_morphism) where T

Given a vector `V` of (homogeneous) elements of `F`, return the quotient of `F` by the (graded) submodule of `F` which is generated by these elements.

Put more precisely, return this quotient as an object of type `SubquoModule`. 

Additionally, if `N` denotes this object,

- return the projection map `F` $\to$ `N` if `task = :with_morphism` (default),
- return and cache the projection map `F` $\to$ `N` if `task = :cache_morphism`,
- do none of the above if `task = :none` or `task = :module`.

If `task = :only_morphism`, return only the projection map.
"""
function quo(F::FreeMod{T}, V::Vector{<:FreeModElem{T}}, task::Symbol = :with_morphism) where T
  S = SubquoModule(F, basis(F))
  Q = SubquoModule(S, V)

  return return_quo_wrt_task(F, Q, task)
end

@doc raw"""
    quo(F::FreeMod{T}, A::MatElem{T}, task::Symbol = :with_morphism) where {T}

Given a (homogeneous) matrix `A`, return the quotient of `F` by the graded submodule of `F` which is generated by 
the rows of `A`.

Put more precisely, return this quotient as an object of type `SubquoModule`. 

Additionally, if `N` denotes this object,

- return the projection map `F` $\to$ `N` if `task = :with_morphism` (default),
- return and cache the projection map `F` $\to$ `N` if `task = :cache_morphism`,
- do none of the above if `task = :none` or `task = :module`.

If `task = :only_morphism`, return only the projection map.
"""
function quo(F::FreeMod{T}, A::MatElem{T}, task::Symbol = :with_morphism) where {T}
  E = identity_matrix(base_ring(F), rank(F))
  Q = SubquoModule(F, E, A)

  return return_quo_wrt_task(F, Q, task)
end

@doc raw"""
    quo(F::FreeMod{T}, O::Vector{<:SubquoModuleElem{T}}, task::Symbol = :with_morphism) where T

Compute $F / T$, where $T$ is generated by $O$.
The embedding free module of the parent of the elements of `O` must be `F`.
If `task` is set to `:with_morphism` (default option) or to `:both` return also the 
canonical projection morphism $F \to F/T$.
If `task` is set to `:cache_morphism` the morphism is also cached.
If `task` is set to `:only_morphism` return only the morphism.
If `task` is set to `:none` or `:module` return only the module.
"""
function quo(F::FreeMod{T}, O::Vector{<:SubquoModuleElem{T}}, task::Symbol = :with_morphism) where T
  S = SubquoModule(F, basis(F))
  Q = SubquoModule(S, [repres(x) for x = O])
  
  return return_quo_wrt_task(F, Q, task)
end

@doc raw"""
    quo(F::SubquoModule{T}, O::Vector{<:FreeModElem{T}}, task::Symbol = :with_morphism) where T

Compute $F / T$, where $T$ is generated by $O$.
The elements of `O` must be elements of the embedding free module of `S`.
If `task` is set to `:with_morphism` (default) or to `:both` return also the 
canonical projection morphism $F \to F/T$.
If `task` is set to `:cache_morphism` the morphism is also cached.
If `task` is set to `:only_morphism` return only the morphism.
If `task` is set to `:none` or `:module` return only the module.
"""
function quo(F::SubquoModule{T}, O::Vector{<:FreeModElem{T}}, task::Symbol = :with_morphism) where T
  if length(O) > 0
    @assert parent(O[1]) === F.F
  end
  if isdefined(F, :quo)
    oscar_assure(F.quo.gens)
    singular_assure(F.quo.gens)
    s = Singular.Module(base_ring(F.quo.gens.SF), [F.quo.gens.SF(x) for x = [O; oscar_generators(F.quo.gens)]]...)
    Q = SubquoModule(F.F, singular_generators(F.sub.gens), s)
    return return_quo_wrt_task(F, Q, task)
  end
  Q = SubquoModule(F, O)
  return return_quo_wrt_task(F, Q, task)
end

@doc raw"""
    quo(M::SubquoModule{T}, V::Vector{<:SubquoModuleElem{T}}, task::Symbol = :with_morphism) where T

Given a vector `V` of (homogeneous) elements of `M`, return the quotient of `M` by the (graded) submodule of `M` which is generated by these elements.

Put more precisely, return the quotient as an object of type `SubquoModule`. 

Additionally, if `N` denotes this object,

- return the projection map `M` $\to$ `N` if `task = :with_morphism` (default),
- return and cache the projection map `M` $\to$ `N` if `task = :cache_morphism`,
- do none of the above if `task = :none` or `task = :module`.

If `task = :only_morphism`, return only the projection map.
"""
function quo(M::SubquoModule{T}, V::Vector{<:SubquoModuleElem{T}}, task::Symbol = :with_morphism) where T
  return quo(M, [repres(x) for x = V], task)
end

@doc raw"""
    quo(M::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}, task::Symbol = :with_morphism) where T

Given a vector `V` of (homogeneous) elements of `M`, return the quotient of `M` by the (graded) submodule of `M` which is generated by these elements.

Put more precisely, return the quotient as an object of type `SubquoModule`. 

Additionally, if `N` denotes this object,

- return the projection map `M` $\to$ `N` if `task = :with_morphism` (default),
- return and cache the projection map `M` $\to$ `N` if `task = :cache_morphism`,
- do none of the above if `task = :none` or `task = :module`.

If `task = :only_morphism`, return only the projection map.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

julia> F = free_module(R, 1);

julia> V = [x^2*F[1]; y^3*F[1]; z^4*F[1]];

julia> N, proj = quo(F, V);

julia> N
Subquotient of Submodule with 1 generator
1 -> e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]

julia> proj
Map with following data
Domain:
=======
Free module of rank 1 over Multivariate polynomial ring in 3 variables over QQ
Codomain:
=========
Subquotient of Submodule with 1 generator
1 -> e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]
```
"""
function quo(M::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}, task::Symbol = :with_morphism) where T
 error("quo is not implemented for the given types.")
end


@doc raw"""
    quo(M::SubquoModule{T}, U::SubquoModule{T}, task::Symbol = :with_morphism) where T

Return the quotient $M / U$.

Put more precisely, if `N` denotes this quotient, return `N` as an object of type `SubquoModule`. Additionally,

- return the projection map `M` $\to$ `N` if `task = :with_morphism` (default),
- return and cache the projection map `M` $\to$ `N` if `task = :cache_morphism`,
- do none of the above if `task = :none` or `task = :module`.

If `task = :only_morphism`, return only the projection map.
"""
function quo(M::SubquoModule{T}, U::SubquoModule{T}, task::Symbol = :with_morphism) where T
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
  return return_quo_wrt_task(M, Q, task)
end

function quo(M::SubquoModule{T}, U::SubquoModule{T}, task::Symbol = :with_morphism) where {T<:MPolyRingElem}
  if isdefined(M, :quo) && isdefined(U, :quo)
    @assert M.quo == U.quo
  else
    @assert !isdefined(M, :quo) && !isdefined(U, :quo)
  end
  Q = SubquoModule(M, oscar_generators(U.sub.gens))
  return return_quo_wrt_task(M, Q, task)
end

@doc raw"""
    quo(F::FreeMod{R}, T::SubquoModule{R}, task::Symbol = :with_morphism) where R

Compute $F / T$.
If `task` is set to `:with_morphism` (default option) or to `:both` return also the 
canonical projection morphism $F \to F/T$.
If `task` is set to `:cache_morphism` the morphism is also cached.
If `task` is set to `:only_morphism` return only the morphism.
If `task` is set to `:none` or `:module` return only the module.
"""
function quo(F::FreeMod{R}, T::SubquoModule{R}, task::Symbol = :with_morphism) where R
  @assert !isdefined(T, :quo)
  return quo(F, gens(T), task)
end

@doc raw"""
    return_quo_wrt_task(M::ModuleFP, Q::ModuleFP, task)

This helper function returns the module `Q = M / N` for some `N` 
along with the canonical projection morphism $M \to Q$ according to the 
given `task`.
"""
function return_quo_wrt_task(M::ModuleFP, Q::ModuleFP, task)
  if task == :none || task == :module
    return Q
  else
    pro = hom(M, Q, gens(Q); check=false)
    pro.generators_map_to_generators = true # Makes evaluation of the inclusion easier
    task == :cache_morphism && register_morphism!(pro)
    task == :only_morphism && return pro
    return Q, pro
  end
end

@doc raw"""
    syzygy_module(F::ModuleGens; sub = FreeMod(base_ring(F.F), length(oscar_generators(F))))
"""
function syzygy_module(F::ModuleGens{T}; sub = FreeMod(base_ring(F.F), length(oscar_generators(F)))) where {T <: MPolyRingElem}
  singular_assure(F)
  # TODO Obtain the Gröbner basis and cache it
  s = Singular.syz(singular_generators(F))
  return SubquoModule(sub, s)
end

@doc raw"""
    gens(M::SubquoModule{T}) where T

Return the generators of `M`.
"""
function gens(M::SubquoModule{T}) where T
  return SubquoModuleElem{T}[gen(M,i) for i=1:ngens(M)]
end

@doc raw"""
    gen(M::SubquoModule{T}, i::Int) where T

Return the `i`th generator of `M`.
"""
function gen(M::SubquoModule{T}, i::Int) where T
  R = base_ring(M)
  v::SRow{T} = sparse_row(R)
  v.pos = [i]
  v.values = [R(1)]
  return SubquoModuleElem{T}(v, M)
end

@doc raw"""
    ngens(M::SubquoModule)

Return the number of generators of `M`.
"""
ngens(M::SubquoModule) = ngens(M.sub)

@doc raw"""
    base_ring(M::SubquoModule)

Given an `R`-module `M`, return `R`.
"""
base_ring(M::SubquoModule) = base_ring(M.F)::base_ring_type(M.F)

@doc raw"""
    zero(M::SubquoModule)

Return the zero element of `M`.
"""
zero(M::SubquoModule) = SubquoModuleElem(SRow(base_ring(M)), M)

@doc raw"""
    is_zero(M::SubquoModule)

Return `true` if `M` is the zero module, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> F = free_module(R, 1)
Free module of rank 1 over Multivariate polynomial ring in 3 variables over QQ

julia> A = R[x^2+y^2;]
[x^2 + y^2]

julia> B = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = SubquoModule(F, A, B)
Subquotient of Submodule with 1 generator
1 -> (x^2 + y^2)*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]

julia> is_zero(M)
false
```
"""
function is_zero(M::SubquoModule)
  return all(iszero, gens(M))
end

@doc raw"""
    getindex(F::SubquoModule, i::Int)

Return the `i`th generator of `F`.
"""
function getindex(F::SubquoModule, i::Int)
  i == 0 && return zero(F)
  return gen(F, i)
end

function iterate(F::ModuleGens, i::Int = 1)
  if i>length(F)
    return nothing
  else
    return F[i], i+1
  end
end
eltype(::ModuleGens{T}) where {T} = FreeModElem{T} 

#??? A scalar product....
function *(a::FreeModElem, b::Vector{FreeModElem})
  @assert dim(parent(a)) == length(b)
  s = zero(parent(a))
  for (p,v) = a.coords
    s += v*b[p]
  end
  return s
end

@doc raw"""
    is_zero(m::SubquoModuleElem)

Return `true` if `m` is zero, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> F = free_module(R, 1)
Free module of rank 1 over Multivariate polynomial ring in 3 variables over QQ

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = SubquoModule(F, A, B)
Subquotient of Submodule with 2 generators
1 -> x*e[1]
2 -> y*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]

julia> is_zero(M[1])
false

julia> is_zero(x*M[1])
true
```

```jldoctest
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> F = graded_free_module(Rg, 1)
Graded free module Rg^1([0]) of rank 1 over Rg

julia> A = Rg[x; y]
[x]
[y]

julia> B = Rg[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = SubquoModule(F, A, B)
Graded subquotient of submodule of F generated by
1 -> x*e[1]
2 -> y*e[1]
by submodule of F generated by
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]

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

function iszero(m::SubquoModuleElem{<:MPolyRingElem})
  C = parent(m)
  if !isdefined(C, :quo)
    return iszero(repres(m))
  end
  x = reduce(repres(m), C.quo)
  return iszero(x)
end


