@doc raw"""
    is_equal_with_morphism(M::SubquoModule{T}, N::SubquoModule{T}, task::Symbol = :none) where {T}

If $M = N$ (mathematically, but with (possibly) different generating systems), return $\phi : M \to N$ 
which is mathematically the identity. 
If `task == :inverse` also the inverse map is computed and cached (in the morphism).
If `task == :cache_morphism` the inverse map is also cached in `M` and `N`.

# Examples
```jldoctest
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> F = graded_free_module(Rg, 1);

julia> A = Rg[x; y];

julia> B = Rg[x^2; y^3; z^4];

julia> M = SubquoModule(F, A, B);

julia> is_equal_with_morphism(M, M)
Homogeneous module homomorphism
  from M
  to M
defined by
  x*e[1] -> x*e[1]
  y*e[1] -> y*e[1]
```
"""
function is_equal_with_morphism(M::SubquoModule{T}, N::SubquoModule{T}, task::Symbol = :none) where {T}
  @assert M == N

  M_to_N = hom(M, N, Vector{elem_type(N)}([SubquoModuleElem(coordinates(repres(m), N), N) for m in gens(M)]))

  if task == :cache_morphism || task == :inverse
    N_to_M = hom(N, M, Vector{elem_type(M)}([SubquoModuleElem(coordinates(repres(n), M), M) for n in gens(N)]))
    M_to_N.inverse_isomorphism = N_to_M
    N_to_M.inverse_isomorphism = M_to_N

    if task == :cache_morphism
      register_morphism!(M_to_N) 
      register_morphism!(N_to_M)
    end
  end
  
  return M_to_N
end


function Hecke.ring(I::MPolyIdeal)
  return base_ring(I)
end

# We can not use the signature with T because the MPolyQuoIdeals are 
# not parametrized by the element type of their ring.
#function *(I::Ideal{T}, M::ModuleFP{T}) where {T<:RingElem}
function *(I::Ideal, M::ModuleFP)
  base_ring(I) === base_ring(M) || error("ideal and module are not defined over the same ring")
  return sub(M, elem_type(M)[g*e for g in gens(I) for e in gens(M)])
end

@doc raw"""
    ideal_to_module(I::MPolyIdeal{T}, F::FreeMod{T})

Convert the ideal `I` into a submodule of `F`. If the rank of `F`
is not equal to 1 an error is thrown.
"""
function ideal_to_module(I::MPolyIdeal{T}, F::FreeMod{T}) where T
  @assert rank(F) == 1
  return sub_object(F, [x*gen(F,1) for x in gens(I)])
end

@doc raw"""
    ideal_to_module(I::MPolyIdeal)

Convert the ideal `I` into a submodule.
"""
function ideal_to_module(I::MPolyIdeal)
  F = FreeMod(base_ring(I), 1)
  return ideal_to_module(I, F)
end

#############################
# TODO ?
#############################
@doc raw"""
    find_sequence_of_morphisms(N::SubquoModule, M::SubquoModule)

Compute a path from `N` to `M` in the graph of cached (canonical) morphisms.
Return it as Vector of maps where the first element has domain `N` and the last 
has codomain `M`.
If there exists no path an error is thrown.
"""
function find_sequence_of_morphisms(N::SubquoModule, M::SubquoModule)
  if M===N
    return [id_hom(M)]
  end
  parent_hom = IdDict{SubquoModule, ModuleFPHom}()
  modules = [M]
  found_N = false
  for A in modules
    if N in keys(A.incoming)
      parent_hom[N] = _recreate_morphism(N, A, A.incoming[N])
      found_N = true
      break
    end

    for (I, f) in A.incoming
      I === A && continue
      parent_hom[I] = _recreate_morphism(I, A, f)
      push!(modules, I)
    end
  end

  if !found_N
    throw(DomainError("There is no path of canonical homomorphisms between the modules!"))
  end
  morphisms = Vector{ModuleFPHom}()
  A = N
  while A !== M
    f = parent_hom[A]
    push!(morphisms, f)
    A = codomain(f)
  end
  return morphisms
end

function _recreate_morphism(dom::ModuleFP, cod::ModuleFP, t::Tuple{<:SMat, <:Any})
  A, bc = t
  if bc === nothing
    return hom(dom, cod, [sum(a*cod[i] for (i, a) in v; init=zero(cod)) for v in A], check=false)
  else
    return hom(dom, cod, [sum(a*cod[i] for (i, a) in v; init=zero(cod)) for v in A], bc, check=false)
  end
end

@doc raw"""
    transport(M::SubquoModule, v::SubquoModuleElem)

Map the element `v` to an element of the module `M` using cached 
canonical homomorphisms between the parent Module of `v` and `M`.
If this is not possible an error is thrown.
"""
function transport(M::SubquoModule, v::SubquoModuleElem)
  N = parent(v)
  morphisms = find_sequence_of_morphisms(N, M)

  return foldl(|>, morphisms; init=v)
end

@doc raw"""
    find_morphism(M::SubquoModule, N::SubquoModule)

Find a morphism from `M` to `N` in the graph of cached (canonical) morphisms
and return it.
If there exists no such morphism an error is thrown.
"""
function find_morphism(M::SubquoModule, N::SubquoModule)
  morphisms = find_sequence_of_morphisms(M, N)
  return reduce(*, morphisms)
end

@doc raw"""
    find_morphisms(N::SubquoModule, M::SubquoModule)

Traverse the graph of all cached morphisms originating in `N` and 
return the compositions of all loop-free paths from `N` to `M`.
"""
function find_morphisms(N::SubquoModule, M::SubquoModule)
  # from N to M

  all_paths = []

  function helper_dfs!(U::SubquoModule, D::SubquoModule, visited::Vector{<:ModuleFP}, path::Vector)
    if U === D
      push!(all_paths, path)
      return
    end
    for (cod, neighbor_morphism) in U.outgoing
      any(x->x===cod, visited) && continue
      helper_dfs!(cod, D, push!(visited, cod), union(path, [_recreate_morphism(U, cod, neighbor_morphism)]))
    end
  end

  helper_dfs!(N, M, Vector{ModuleFP}(), [])

  morphisms = Vector{ModuleFPHom}()
  for path in all_paths
    phi = id_hom(N)
    for h in path
      phi = phi*h
    end
    push!(morphisms, phi)
  end

  return morphisms
end

#############################
# Useful functions
#############################

@doc raw"""
    register_morphism!(f::ModuleFPHom)

Cache the morphism `f` in the corresponding caches of the domain and codomain of `f`.
"""
function register_morphism!(f::ModuleFPHom)
  dom = domain(f)
  cod = codomain(f)
  dom.outgoing[cod] = sparse_matrix(f), ring_map(f)
  cod.incoming[dom] = sparse_matrix(f), ring_map(f)
  f
end

# Some missing methods for the above to work
function sparse_matrix(f::ModuleFPHom)
  dom = domain(f)
  R = base_ring(codomain(f))
  result = sparse_matrix(R, 0, ngens(codomain(f)))
  for v in gens(dom)
    push!(result, coordinates(f(v)))
  end
  return result
end

ring_map(f::FreeModuleHom{<:AbstractFreeMod, <:ModuleFP, Nothing}) = nothing
ring_map(f::FreeModuleHom) = f.ring_map

ring_map(f::SubQuoHom{<:AbstractFreeMod, <:ModuleFP, Nothing}) = nothing
ring_map(f::SubQuoHom) = f.ring_map

function default_ordering(F::FreeMod)
  if !isdefined(F, :default_ordering)
    if iszero(F)
      F.default_ordering = default_ordering(base_ring(F))*ModuleOrdering(F, Orderings.ModOrdering(Vector{Int}(), :lex))
    else
      F.default_ordering = default_ordering(base_ring(F))*lex(gens(F))
    end
  end
  return F.default_ordering::ModuleOrdering{typeof(F)}
end

function default_ordering(F::FreeMod{T}) where {T<:Union{ZZRingElem, FieldElem}}
    if !isdefined(F, :default_ordering)
        if iszero(F)
            F.default_ordering = ModuleOrdering(F, Orderings.ModOrdering(Int[], :lex))
        else
            F.default_ordering = lex(gens(F))
        end
    end
    return F.default_ordering::ModuleOrdering{typeof(F)}
end


##############################
#TODO: move to Singular.jl ?

function _reduce(a::Singular.smodule, b::Singular.smodule)
  @assert b.isGB
  p = Singular.libSingular.p_Reduce(a.ptr, b.ptr, base_ring(b).ptr)
  return Singular.Module(base_ring(b), p)
end

function _reduce(a::Singular.svector, b::Singular.smodule)
  @assert b.isGB
  p = _reduce(Singular.Module(base_ring(b), a), b)[1] # TODO When available in Singular.jl use reduce(::svector, ::smodule)
  return Singular.Module(base_ring(b), p)[1]
end

#TODO: tensor_product from Raul's H is broken

#############################################
# Test hom
#############################################
function to_julia_matrix(A::Union{MatElem})
  return eltype(A)[A[i, j] for i = 1:nrows(A), j = 1:ncols(A)]
end

function copy_and_reshape(M::MatElem, n, m)
  julia_matrix = to_julia_matrix(M)
  julia_matrix = reshape(julia_matrix, n, m)
  R = base_ring(M)
  return matrix(R, n, m, julia_matrix)
end

function hom_matrices_helper(f1::MatElem{T}, g1::MatElem{T}) where T
  R = base_ring(f1)
  s1, s0 = size(f1)
  t1, t0 = size(g1)

  g2 = matrix_kernel(g1)
  n = s1*t0
  m = s0*t0 + s1*t1
  delta::MatrixElem{T} = zero_matrix(R, m,n)
  for j=1:s0*t0
    b_vector::MatrixElem{T} = zero_matrix(R, 1,s0*t0)
    b_vector[1,j] = R(1)
    A = copy_and_reshape(b_vector, s0, t0)
    res = f1*A
    delta[j,:] = copy_and_reshape(res, 1, n)
  end
  for j=s0*t0+1:m
    b_vector::MatrixElem{T} = zero_matrix(R, 1,m-s0*t0)
    b_vector[1,j-s0*t0] = R(1)
    B = copy_and_reshape(b_vector, s1, t1)
    res = -B*g1
    delta[j,:] = copy_and_reshape(res, 1,n)
  end

  gamma = matrix_kernel(delta)
  gamma = gamma[:,1:s0*t0]

  rho::MatrixElem{T} = zero_matrix(R, s0*t1, s0*t0)

  for j=1:s0*t1
    b_vector = zero_matrix(R, 1,s0*t1)
    b_vector[1,j] = R(1)
    C = copy_and_reshape(b_vector, s0, t1)
    res1 = C*g1
    rho[j,1:length(res1)] = copy_and_reshape(res1, 1,length(res1))
  end

  M = SubquoModule(gamma, rho)

  function convert_to_matrix(v::SubquoModuleElem{T})
    if parent(v) !== M
      throw(DomainError("v does not represent a homomorphism"))
    end
    R = base_ring(M)
    c = coordinates(repres(v))
    A = copy_and_reshape(dense_row(c[1:s0*t0], s0*t0), s0, t0)
    return A
  end

  return M, convert_to_matrix
end


# Hom and hom_matrices implement the same mathematical algorithm, but the implementations 
# differ a lot e.g. in the data structures. As a result performance differs depending 
# on the example favoring the one or the other. So it makes sense to offer both. 
# With option :matrices in hom() hom_matrices is used.
@doc raw"""
    hom_matrices(M::SubquoModule{T},N::SubquoModule{T}) where T

Return a subquotient $S$ such that $\text{Hom}(M,N) \cong S$
"""
function hom_matrices(M::SubquoModule{T},N::SubquoModule{T},simplify_task=true) where T
  f1 = generator_matrix(present_as_cokernel(M).quo)
  g1 = generator_matrix(present_as_cokernel(N).quo)
  R = base_ring(M)
  SQ, convert_to_matrix = hom_matrices_helper(f1,g1)
  if simplify_task
    SQ2, i, p = simplify(SQ)
    to_homomorphism = function(elem::SubquoModuleElem{T})
      elem2 = i(elem)
      A = convert_to_matrix(elem2)
      return SubQuoHom(M,N,A; check=false)
    end
    to_subquotient_elem = function(H::ModuleFPHom)
      m = length(matrix(H))
      v = copy_and_reshape(matrix(H),1,m)
      v = FreeModElem(sparse_row(v), FreeMod(R, length(v)))
      return p(SQ(v))
    end

    to_hom_map = MapFromFunc(SQ2, Hecke.MapParent(M, N, "homomorphisms"), to_homomorphism, to_subquotient_elem)
    set_attribute!(SQ2, :hom => (M, N), :module_to_hom_map => to_hom_map)

    return SQ2, to_hom_map
  else
    to_subquotient_elem = function(H::ModuleFPHom)
      m = length(matrix(H))
      v = copy_and_reshape(matrix(H),1,m)
      v = FreeModElem(sparse_row(v), ambient_free_module(SQ))
      return SQ(v)
    end
    to_homomorphism = function(elem::SubquoModuleElem{T})
      A = convert_to_matrix(elem)
      return SubQuoHom(M,N,A; check=false)
    end

    to_hom_map = MapFromFunc(SQ, Hecke.MapParent(M, N, "homomorphisms"), to_homomorphism, to_subquotient_elem)
    set_attribute!(SQ, :hom => (M, N), :module_to_hom_map => to_hom_map)

    return SQ, to_hom_map
  end
end

function change_base_ring(S::Ring, F::FreeMod)
  R = base_ring(F)
  r = ngens(F)
  FS = is_graded(F) ? graded_free_module(S, degrees_of_generators(F)) : FreeMod{elem_type(S)}(ngens(F), S, F.S) # the symbols of F
  map = hom(F, FS, gens(FS), MapFromFunc(R, S, S))
  return FS, map
end

function change_base_ring(f::Map{DomType, CodType}, F::FreeMod) where {DomType<:Ring, CodType<:Ring}
  domain(f) == base_ring(F) || error("ring map not compatible with the module")
  S = codomain(f)
  r = ngens(F)
  FS = is_graded(F) ? graded_free_module(S, degrees_of_generators(F)) : FreeMod{elem_type(S)}(ngens(F), S, F.S)
  map = hom(F, FS, gens(FS), f)
  return FS, map
end

function change_base_ring(S::Ring, M::SubquoModule)
  F = ambient_free_module(M)
  R = base_ring(M)
  FS, mapF = change_base_ring(S, F)
  g = ambient_representatives_generators(M)
  rels = relations(M)
  MS = SubquoModule(FS, mapF.(g), mapF.(rels))
  map = SubQuoHom(M, MS, gens(MS), MapFromFunc(R, S, S); check=false)
  return MS, map
end

function change_base_ring(f::Map{DomType, CodType}, M::SubquoModule) where {DomType<:Ring, CodType<:Ring}
  domain(f) === base_ring(M) || error("ring map not compatible with the module")
  S = codomain(f)
  F = ambient_free_module(M)
  R = base_ring(M)
  FS, mapF = change_base_ring(f, F)
  g = ambient_representatives_generators(M)
  rels = relations(M)
  MS = SubquoModule(FS, mapF.(g), mapF.(rels))
  map = SubQuoHom(M, MS, gens(MS), f; check=false)
  return MS, map
end

function change_base_ring(phi::Any, f::ModuleFPHom; 
    domain_base_change=change_base_ring(phi, domain(f))[2], 
    codomain_base_change=change_base_ring(phi, codomain(f))[2]
  )
  new_dom = codomain(domain_base_change)
  new_cod = codomain(codomain_base_change)
  return hom(new_dom, new_cod, codomain_base_change.(f.(gens(domain(f)))))
end


### Duals of modules
@doc raw"""
    dual(M::ModuleFP; codomain::Union{FreeMod, Nothing}=nothing)

Return a pair ``(M*, i)`` consisting of the dual of ``M`` and its 
interpretation map ``i``, turning an element ``φ`` of ``M*`` into 
a homomorphism ``M → R``. 

The optional argument allows to specify a free module of rank ``1`` 
for the codomain of the dualizing functor.
"""
function dual(M::ModuleFP; codomain::Union{FreeMod, Nothing}=nothing)
  R = base_ring(M)
  codomain = codomain === nothing ? (is_graded(M) ? graded_free_module(R, 1) : FreeMod(R, 1)) : codomain
  base_ring(codomain) === R && rank(codomain) == 1 || error("codomain must be free of rank one over the base ring of the first argument")
  return hom(M, codomain)
end

@doc raw"""
    double_dual(M::ModuleFP)

For a finite ``R``-module ``M`` return a pair ``(M**, ϕ)`` consisting of 
its double dual ``M** = Hom(Hom(M, R), R)`` together with the canonical 
map ``ϕ : M → M**, v ↦ (φ ↦ φ(v)) ∈ Hom(M*, R)``.
"""
function double_dual(M::FreeMod{T}; codomain::Union{FreeMod, Nothing}=nothing, check::Bool=true) where T
  R = base_ring(M)
  codomain = codomain === nothing ? (is_graded(M) ? graded_free_module(R, 1) : FreeMod(R, 1)) : codomain
  M_dual, _ = dual(M, codomain=codomain)
  M_double_dual, _ = dual(M_dual, codomain=codomain)
  if length(gens(M_dual)) == 0
    psi_gens = [zero(M_double_dual) for _ in gens(M)]
  else
    psi_gens = [
      homomorphism_to_element(
        M_double_dual,
        FreeModuleHom(M_dual, codomain, [element_to_homomorphism(phi)(x) for phi in gens(M_dual)]; check)
      )
      for x in gens(M)
    ]
  end
  psi = FreeModuleHom(M, M_double_dual, psi_gens; check)
  return M_double_dual, psi
end

function double_dual(M::SubquoModule{T}; codomain::Union{FreeMod, Nothing}=nothing, check::Bool=true) where T
  R = base_ring(M)
  codomain = codomain === nothing ? (is_graded(M) ? graded_free_module(R, 1) : FreeMod(R, 1)) : codomain
  M_dual, _ = dual(M, codomain=codomain)
  M_double_dual, _ = dual(M_dual, codomain=codomain)
  if length(gens(M_dual)) == 0
    psi_gens = [zero(M_double_dual) for _ in gens(M)]
  else
    psi_gens = [
      homomorphism_to_element(
        M_double_dual,
        SubQuoHom(M_dual, codomain, [element_to_homomorphism(phi)(x) for phi in gens(M_dual)]; check)
      )
      for x in gens(M)
    ]
  end
  psi = SubQuoHom(M, M_double_dual, psi_gens; check)
  return M_double_dual, psi
end


@doc raw"""
    dual(f::ModuleFPHom; codomain::FreeMod)

Given a morphism of modules ``f : M → N``, return the morphism
``fᵀ : N* → M*, φ ↦ (v ↦ φ(f(v)))`` induced on the duals.

The optional argument allows to specify a free module of rank one over the 
base ring of ``f`` for building the duals of ``M`` and ``N``.
"""
function dual(f::ModuleFPHom{<:ModuleFP, <:ModuleFP, Nothing}; # Third parameter assures same base ring
    codomain::FreeMod=FreeMod(base_ring(domain(f)), 1), 
    domain_dual::ModuleFP=dual(Oscar.domain(f), codomain=codomain)[1],
    codomain_dual::ModuleFP=dual(Oscar.codomain(f), codomain=codomain)[1]
  )
  M = Oscar.domain(f)
  N = Oscar.codomain(f)
  R = base_ring(Oscar.domain(f))
  R === base_ring(N) || error("modules must be defined over the same rings")

  M_dual = domain_dual
  N_dual = codomain_dual

  return hom(N_dual, M_dual, 
             [homomorphism_to_element(M_dual, 
                                      hom(M, codomain, 
                                          [element_to_homomorphism(phi)(f(v)) for v in gens(M)]
                                         )
                                     )
              for phi in gens(N_dual)])
end

##########################################################################
## Functionality for modules happening to be finite dimensional vector
## spaces
#
# The general user facing signature is 
#
#   vector_space[_dimension/_basis](kk::Field, M::ModuleFP; check::Bool)
#
# where we assume that either 1) the `base_ring` of `M` is already the 
# field `kk`, or 2) the `base_ring` `R` of `M` is an algebra over the 
# field `kk`. Other cases are possible, but must be caught and delegated 
# to a custom implementation. 
#
# For convenience we also allow 
#
#   vector_space[_dimension/_basis](M::ModuleFP; check::Bool, cached::Bool=true)
#
# which picks the field `kk` according to the above assumptions. Note 
# that this has the default option to cache the result (as it does not depend 
# on a second argument `kk` anymore).
#
# Internally we need to first make sure that the module is finitely presented, 
# i.e. that its `.sub` part is actually the whole `ambient_free_module`. Thus 
# the generic code passes to a presentation and converts the result back if 
# necessary. 
#
# Once we are sure that the module is presented, we delegate to respective 
# internal methods 
#
#   _vector_space[_dimension/_basis](kk::Field, M::ModuleFP; check::Bool)
#
# These might still do internal checks, like e.g. finiteness over `kk`.
#
# For the non-graded polynomial case there are methods to filter out a 
# graded part w.r.t the total degree via 
#
#   vector_space[_dimension/_basis](kk::Field, M::ModuleFP, d::Int; check::Bool)
#
# These will eventually be deprecated to internal methods, but they are part of 
# the system for now. In the graded case we aim to have 
#   
#   vector_space[_dimension/_basis](kk::Field, M::ModuleFP, d::FinGenAbGroupElem; check::Bool)
#
# for which there are stubs at the moment, but only partial implementations.
##########################################################################
@doc raw"""
    vector_space_dim(M::SubquoModule)

Let ``R`` be a `Ring` over a `coefficient_ring` ``k`` which is a field and let 
``M`` be a subquotient module over ``R``. Return the dimension of `M`, seen 
as a the ``k``-vectorspace. 

If `R` itself is a `Field`, this simply returns the dimension of `M`.

# Examples:
```jldoctest
julia> R,(x,y,z,w) = QQ[:x, :y, :z, :w];

julia> F = free_module(R,2);

julia> M,_ = quo(F,[1*gen(F,1),x^2*gen(F,2),y^3*gen(F,2),z*gen(F,2),w*gen(F,2)]);

julia> vector_space_dim(M,1)
2

julia> vector_space_dim(M,2)
2

julia> vector_space_dim(M,3)
1

julia> vector_space_dim(M)
6

```
"""
function vector_space_dim(M::SubquoModule; check::Bool=true, cached::Bool=true)
  # Per default we assume that the user means the vector space dimension 
  # over the `base_ring` of the ring on which `M` is defined.
  R = base_ring(M)
  @assert coefficient_ring(R) isa Field "`coefficient_ring` of the ring over which the module is defined is not a field"
  cached && has_attribute(M, :vector_space_dimension) && return get_attribute(M, :vector_space_dimension)::Union{Int, PosInf}
  res = vector_space_dim(coefficient_ring(R), M; check)
  cached && set_attribute!(M, :vector_space_dimension=>res)
  return res
end

# If the module's ring is already a field, we compute the dimension over that.
function vector_space_dim(M::SubquoModule{T}; check::Bool=true, cached::Bool=true) where {T <: FieldElem}
  return vector_space_dim(base_ring(M), M; check)
end

# Syntax to be coherent with other methods for `vector_space_dim`.
function vector_space_dim(kk::Field, M::SubquoModule; check::Bool=true)
  R = base_ring(M)
  kk === R || kk === coefficient_ring(R) || error("not implemented over fields different from the ground ring or the `coefficient_ring` thereof")

  S = base_ring(M)
  F = ambient_free_module(M)
  # We need `M` to be presented  
  if !((ngens(M) == ngens(F)) && all(repres(v) == e for (v, e) in zip(gens(M), gens(F))))
    # pres = presentation(M)
    # MM = cokernel(map(pres, 1))    # cokernel sometimes returns also the canonical projection
    MM = present_as_cokernel(M)
    return _vector_space_dim(kk, MM; check)
  end
  # At this point we may assume `M` to be presented
  return _vector_space_dim(kk, M; check)
end

# Assumes `M` to be presented
function _vector_space_dim(kk::Field, M::SubquoModule; check::Bool=true)
  _is_finite(kk, M) || return inf

  # The generic implementation just takes the length of a basis. 
  # This might not be efficient, so consider overwriting it in specific cases.
  return length(vector_space_basis(kk, M; check))
end

@doc raw"""
    _is_finite(kk::Field, M::SubquoModule)

Return whether `M` is a finite vector space over `kk`.
Note that this assumes the `ambient_representatives_generators` of `M` to 
coincide with the generators of the `ambient_free_module`. This internal 
function must only be called in such appropriate settings. 
"""
function _is_finite(kk::Field, M::SubquoModule)
  error("finiteness check for $M over $kk not implemented")
end

function _is_finite(kk::Field, M::SubquoModule{T}) where {T<:MPolyRingElem{<:FieldElem}}
  @assert kk === coefficient_ring(base_ring(M)) "not implemented for fields other than the `coefficient_ring` of the `base_ring` of the module"
  is_zero(M) && return true
  !isdefined(M, :quo) && return false
  return has_monomials_on_all_axes(leading_module(M.quo, default_ordering(M)))
end

function _is_finite(kk::Field, M::SubquoModule{<:FieldElem})
  kk === base_ring(M) && return true
  error("finiteness check for $M over $kk not implemented")
end

@doc raw""" 
    vector_space_basis(M::SubquoModule; cached::Bool=true)

Assume `M` to be defined over a `kk`-algebra `R` which is not a field. 
Return a list of elements in `M` which form a basis for `M` when considered 
as a vector space over the `coefficient_ring` `kk` of its `base_ring` `R`.
"""
function vector_space_basis(M::SubquoModule; cached::Bool=true, check::Bool=true)
  cached && has_attribute(M, :vector_space_basis) && return get_attribute(M, :vector_space_basis)::Vector{elem_type(M)}
  R = base_ring(M)
  kk = coefficient_ring(R)
  result = vector_space_basis(kk, M; check)
  if cached
    set_attribute!(M, :vector_space_basis=>result)
  end
  return result
end

function vector_space_basis(M::SubquoModule{<:FieldElem}; cached::Bool=true, check::Bool=true)
  cached && has_attribute(M, :vector_space_basis) && return get_attribute(M, :vector_space_basis)::Vector{elem_type(M)}
  kk = base_ring(M)
  result = vector_space_basis(kk, M; check)
  if cached
    set_attribute!(M, :vector_space_basis=>result)
  end
  return result
end

function vector_space_basis(kk::Field, M::SubquoModule{<:FieldElem}; check::Bool=true)
  kk === base_ring(M) || error("not implemented for other fields than the `base_ring` of the module")
  # TODO: look up the existing implementations of rank and put the relevant things here.
  error("not implemented")
end

@doc raw"""
    vector_space_basis(kk::Field, M::SubquoModule)

Return a list of elements in `M` which form a basis of `M` as a vector space over `kk`.
"""
function vector_space_basis(kk::Field, M::SubquoModule; check::Bool=true)
  R = base_ring(M)
  # At the moment we do not support vector space dimensions for modules 
  # in cases different from ``kk`` being the `coefficient_ring` of the `base_ring`
  # For widening the functionality, catch your case here and delegate 
  # to a respective internal method. 
  kk === coefficient_ring(R) || error("`vector_space_basis` not implemented over fields other than the `coefficient_ring` of the ring over which the module is defined")

  S = base_ring(M)
  F = ambient_free_module(M)
  # We need `M` to be presented  
  if !((ngens(M) == ngens(F)) && all(repres(v) == e for (v, e) in zip(gens(M), gens(F))))
    # pres = presentation(M)
    # MM = cokernel(map(pres, 1))    # cokernel sometimes returns also the canonical projection
    MM = present_as_cokernel(M)
    B = _vector_space_basis(kk, MM; check)
    # aug = map(pres, 0)
    pres0 = ambient_free_module(MM)
    aug = hom(pres0, M, gens(M))
    # return elem_type(M)[aug(pres[0](coordinates(v))) for v in B]
    return elem_type(M)[aug(pres0(coordinates(v))) for v in B]
  end
  
  # If execution gets here, `M` is presented.
  result = is_graded(M) ? _vector_space_basis_graded(kk, M; check) : _vector_space_basis(kk, M; check)
  return result
end

# compute a vector space basis over the `coefficient_ring` of the `base_ring`
function _vector_space_basis(kk::Field, M::SubquoModule; check::Bool=true)
  error("not implemented for modules of type $(typeof(M))")
end

function _vector_space_basis_graded(kk::Field, M::SubquoModule; check::Bool=true)
  error("not implemented for modules of type $(typeof(M))")
end

### the ungraded case of polynomial rings over fields
# This is an internal method which assumes `M` to be presented.
function _vector_space_basis(kk::Field, M::SubquoModule{T}; check::Bool=true) where {T <: MPolyRingElem{<:FieldElem}}
  R = base_ring(M)
  @assert kk === coefficient_ring(R) "not implemented for other fields than the coefficients of the polynomial ring"
  @check _is_finite(kk, M) "module is not finite over the given field"
  F = ambient_free_module(M)
  if !isdefined(M, :quo)
    is_zero(M) || error("vector space basis of an infinite dimensional module can not be computed")
    return elem_type(M)[]
  end

  I = M.quo

  o = default_ordering(M)
  LM = leading_module(I, o)
  
  d = 0
  inc = _vector_space_basis(kk, M, 0; check)
  result = elem_type(M)[]

  while !isempty(inc)
    result = vcat(result, inc)
    d = d + 1
    inc = _vector_space_basis(kk, M, d; check)
  end
  return result
end

function vector_space_dim(M::SubquoModule, d::Union{FinGenAbGroupElem, Int64})
  R = base_ring(M)
  kk = coefficient_ring(R)
  return vector_space_dim(kk, M, d)
end

function vector_space_dim(M::SubquoModule{<:FieldElem}, d::Union{FinGenAbGroupElem, Int64})
  kk = base_ring(M)
  return vector_space_dim(kk, M, d)
end

function vector_space_dim(kk::Field, M::SubquoModule, d::Union{FinGenAbGroupElem, Int64})
  return length(vector_space_basis(kk, M, d))
end

@doc raw"""
    vector_space_basis(M::SubquoModule, d::Union{FinGenAbGroupElem, Int64}; check::Bool=true)

Let ``R`` be the ring over which ``M`` is defined. We assume ``R`` to be 
a ``kk``-algebra over some field ``kk``. 

If ``M`` is graded, return a list of elements which form a basis over ``kk``
for the graded piece of ``M`` of degree `d`.

If ``M`` is not graded, we restrict to the case of where ``R`` is a polynomial 
ring over ``kk`` and consider the natural grading by total degree on both ``R`` 
and the `ambient_free_module` of `M`. 

# Examples:
```jldoctest
julia> R,(x,y,z,w) = QQ[:x, :y, :z, :w];

julia> F = free_module(R,2);

julia> M,_ = quo(F,[1*gen(F,1),x^2*gen(F,2),y^3*gen(F,2),z*gen(F,2),w*gen(F,2)]);

julia> vector_space_basis(M,2)
2-element Vector{SubquoModuleElem{QQMPolyRingElem}}:
 x*y*e[2]
 y^2*e[2]

julia> vector_space_basis(M,0)
1-element Vector{SubquoModuleElem{QQMPolyRingElem}}:
 e[2]

julia> vector_space_basis(M)
6-element Vector{SubquoModuleElem{QQMPolyRingElem}}:
 e[2]
 x*e[2]
 y*e[2]
 x*y*e[2]
 y^2*e[2]
 x*y^2*e[2]

```
"""
function vector_space_basis(M::SubquoModule, d::Union{FinGenAbGroupElem, Int64}; check::Bool=true)
  kk = coefficient_ring(base_ring(M))
  return vector_space_basis(kk, M, d; check)
end

function vector_space_basis(kk::Field, M::SubquoModule, d::Union{FinGenAbGroupElem, Int64}; check::Bool=true)
  @assert kk === coefficient_ring(base_ring(M)) "not implemented for fields other than the `coefficient_ring` of the `base_ring` of the module"
  return is_graded(M) ? _vector_space_basis_graded(kk, M, d; check) : _vector_space_basis(kk, M, d; check)
end

function vector_space_basis(M::SubquoModule{T}, d::Union{FinGenAbGroupElem, Int64}; check::Bool=true) where {T<:FieldElem}
  kk = base_ring(M)
  return vector_space_basis(kk, M, d; check)
end

function vector_space_basis(kk::Field, M::SubquoModule{T}, d::Union{FinGenAbGroupElem, Int64}; check::Bool=true) where {T<:FieldElem}
  @assert kk === base_ring(M) "not implemented for fields other than the `baser_ring` of the module" 
  return is_graded(M) ? _vector_space_basis_graded(kk, M, d; check) : _vector_space_basis(kk, M, d; check)
end

function _vector_space_basis(kk::Field, M::SubquoModule, d::FinGenAbGroupElem; check::Bool=true)
  error("module needs to be graded")
end

function _vector_space_basis(kk::Field, M::SubquoModule{T}, d::Int64; check::Bool=true) where {T <: MPolyRingElem{<:FieldElem}}
  R = base_ring(M)
  F = ambient_free_module(M)
  Mq,_ = sub(F,rels(M))

  o = default_ordering(M)
  LM = leading_module(Mq,o)

  mons = [a*e for (a, e) in Iterators.product(monomials_of_degree(R, d), gens(F))]
  return [M(mon) for mon in mons if !(mon in LM)]
end
  
function _vector_space_basis_graded(kk::Field, M::SubquoModule, d::FinGenAbGroupElem; check::Bool=true)
  error("not implemented")
end

function _vector_space_basis_graded(kk::Field, M::SubquoModule, d::Int64; check::Bool=true)
  error("not implemented")
end
  
# function vector_space_dim(M::SubquoModule{T}
#   ) where {T<:MPolyLocRingElem{<:Field, <:FieldElem, <:MPolyRing, <:MPolyRingElem, 
#                                <:MPolyComplementOfKPointIdeal}}
#   F = ambient_free_module(M)
#   Mq,_ = sub(F,rels(M))

#   ambient_representatives_generators(M) == gens(F) || error("not implemented for M/N with non-trivial M")

#   M_shift,_,_ = shifted_module(Mq)
#   o = negdegrevlex(base_ring(M_shift))*lex(ambient_free_module(M_shift))
#   LM = leading_module(M_shift,o)
#   return vector_space_dim(quo_object(ambient_free_module(LM),gens(LM)))
# end

function vector_space_dim(
    M::SubquoModule{T}, d::Int64;
    check::Bool=true
  ) where {T<:MPolyLocRingElem{<:Field, <:FieldElem, <:MPolyRing, <:MPolyRingElem, 
                               <:MPolyComplementOfKPointIdeal}}
  F = ambient_free_module(M)
  Mq,_ = sub(F,rels(M))

  ambient_representatives_generators(M) == gens(F) || error("not implemented for M/N with non-trivial M")

  M_shift,_,_ = shifted_module(Mq)
  o = negdegrevlex(base_ring(M_shift))*lex(ambient_free_module(M_shift))
  LM = leading_module(M_shift,o)
  return vector_space_dim(quo_object(ambient_free_module(LM),gens(LM)),d)
end

# function vector_space_dim(M::SubquoModule{T}; check::Bool=true
#   ) where {T<:MPolyLocRingElem}
#   error("only available in global case and for localization at a point")
# end

function vector_space_dim(M::SubquoModule{T}, d::Int64
  ) where {T<:MPolyLocRingElem}
  error("only available in global case and for localization at a point")
end

# function vector_space_basis(M::SubquoModule{T}
#   ) where {T<:MPolyLocRingElem{<:Field, <:FieldElem, <:MPolyRing, <:MPolyRingElem, 
#                                <:MPolyComplementOfKPointIdeal}}
#   F = ambient_free_module(M)
#   Mq,_ = sub(F,rels(M))

#   ambient_representatives_generators(M) == gens(F) || error("not implemented for M/N with non-trivial M")

#   M_shift,_,_ = shifted_module(Mq)
#   if isdefined(F,:ordering) && is_local(F.ordering)
#     o = F.ordering
#   else
#     o = negdegrevlex(base_ring(M_shift))*lex(ambient_free_module(M_shift))
#   end
#   LM = leading_module(M_shift,o)

#   return vector_space_basis(quo_object(ambient_free_module(LM),gens(LM)))
# end

function vector_space_basis(M::SubquoModule{T},d::Int64
  ) where {T<:MPolyLocRingElem{<:Field, <:FieldElem, <:MPolyRing, <:MPolyRingElem, 
                               <:MPolyComplementOfKPointIdeal}}
  F = ambient_free_module(M)
  Mq,_ = sub(F,rels(M))

  ambient_representatives_generators(M) == gens(F) || error("not implemented for M/N with non-trivial M")

  M_shift,_,_ = shifted_module(Mq)
  if isdefined(F,:ordering) && is_local(F.ordering)
    o = F.ordering
  else
    o = negdegrevlex(base_ring(M_shift))*lex(ambient_free_module(M_shift))
  end
  LM = leading_module(M_shift,o)

  return vector_space_basis(quo_object(ambient_free_module(LM),gens(LM)),d)
end

# function vector_space_basis(M::SubquoModule{T}
#   ) where {T<:MPolyLocRingElem}
#   error("only available in global case and for localization at a point")
# end

function vector_space_basis(M::SubquoModule{T},d::Int64
  ) where {T<:MPolyLocRingElem}
  error("only available in global case and for localization at a point")
end

@doc raw"""
    has_monomials_on_all_axes(M::SubquoModule)

Internal function to test whether M is finite-dimensional vector space. Do not use directly
"""
function has_monomials_on_all_axes(M::SubquoModule)
  length(rels(M)) == 0 || error("not implemented for quotients")
  return has_monomials_on_all_axes(M.sub)
end

function has_monomials_on_all_axes(M::SubModuleOfFreeModule)
  R = base_ring(M)
  
  ambient_rank = ngens(ambient_free_module(M))
  genlist = gens(M)
  explist = Tuple{Vector{Int64}, Int64, Int}[]
  for x in genlist
    tempexp = leading_exponent(x)
    tempdeg = sum(tempexp[1])
    push!(explist,(tempexp[1],tempexp[2],tempdeg))
  end
  for i in 1:ngens(R), j in 1:ambient_rank
    if !any(x -> (x[1][i] == x[3] && x[2]==j), explist)
      return false
    end
  end
  return true
end

### Some missing functionality
function Base.:*(k::Int, f::ModuleFPHom)
  return base_ring(codomain(f))(k)*f
end

function _is_tensor_product(M::ModuleFP)
  !has_attribute(M, :tensor_product) && return false, [M]
  return true, get_attribute(M, :tensor_product)::Tuple
end

function tensor_pure_function(M::ModuleFP)
  success, facs = _is_tensor_product(M)
  success || error("not a tensor product")
  return get_attribute(M, :tensor_pure_function)
end

function tensor_generator_decompose_function(M::ModuleFP)
  success, facs = _is_tensor_product(M)
  success || error("not a tensor product")
  return get_attribute(M, :tensor_generator_decompose_function)
end
  
function tensor_product(mod::Vector{<:ModuleFP})
  return tensor_product(mod...)
end

function tensor_product(f::ModuleFPHom...)
  return tensor_product(collect(f))
end

# We follow the convention to use the same function also for the 
# constructor of induced maps.
tensor_product(dom::ModuleFP, cod::ModuleFP, maps::Vector{<:ModuleFPHom}) = hom_tensor(dom, cod, maps)

function tensor_product(maps::Vector{<:ModuleFPHom}; 
        domain::ModuleFP = tensor_product([domain(f) for f in maps]), 
        codomain::ModuleFP = tensor_product([codomain(f) for f in maps])
    )
  return tensor_product(domain, codomain, maps)
end


# turn a module into a vector space
function vector_space(M::SubquoModule; check::Bool=false)
  R = base_ring(M)
  kk = coefficient_ring(R)
  return vector_space(ModuleFP, kk, M; check)
end

function vector_space(
    ::Type{ResType}, M::SubquoModule; 
    check::Bool=false
  ) where {ResType <: ModuleFP}
  R = base_ring(M)
  kk = coefficient_ring(R)
  return vector_space(ResType, kk, M; check)
end

function vector_space(
    ::Type{ResType}, kk::Field, M::SubquoModule; 
    check::Bool=false
  ) where {ResType <: ModuleFP}
  R = base_ring(M)
  @assert kk === coefficient_ring(R) "not implemented for fields other than the `coefficient_ring` of the `base_ring`"
  @check vector_space_dimension(kk, M) < inf "module is not finite dimensional over the coefficient field"
  B = vector_space_basis(kk, M)
  V = FreeMod(kk, length(B))
  m = MapFromFunc(V, M, x->sum(x[i]*b for (i, b) in enumerate(B); init=zero(M)))
  return V, m
end

# sparse modules over a field return themselves 
function vector_space(
    M::ModuleFP{T}; 
    check::Bool=false
  ) where {T<:FieldElem}
  kk = base_ring(M)
  return vector_space(kk, M; check)
end

function vector_space(
    kk::Field,
    M::ModuleFP{T}; 
    check::Bool=false
  ) where {T<:FieldElem}
  @assert kk === base_ring(M) "not implemented for fields other than the `base_ring`"
  return M, identity_map(M)
end

### functionality for modules over quotient rings
function _vector_space_basis(kk::Field, M::SubquoModule{T}; check::Bool=true) where {T <: MPolyQuoRingElem{<:MPolyRingElem{<:FieldElem}}}
  R = base_ring(M)
  @assert kk === coefficient_ring(R) "not implemented for other fields than the coefficients of the underlying polynomial ring"
  @check _is_finite(kk, M) "module is not finite over the given field"
  B = _vector_space_basis(kk, _as_poly_module(M), check=false)
  is_empty(B) && return elem_type(M)[]
  iota = _iso_with_poly_module(M)
  return iota.(B)
end

function _is_finite(kk::Field, M::SubquoModule{T}) where {T<:MPolyQuoRingElem{<:MPolyRingElem{<:FieldElem}}}
  @assert kk === coefficient_ring(base_ring(M)) "not implemented for fields other than the `coefficient_ring` of the `base_ring` of the module"
  return _is_finite(kk, _as_poly_module(M))
end

### functionality for modules over localized polynomial rings at a point
function _vector_space_basis(kk::Field, M::SubquoModule{T}; check::Bool=true) where {T<:MPolyLocRingElem{<:Field, <:FieldElem, <:MPolyRing, <:MPolyRingElem, <:MPolyComplementOfKPointIdeal}}
  R = base_ring(M)
  @assert kk === coefficient_ring(R) "not implemented for other fields than the coefficients of the underlying polynomial ring"
  @check _is_finite(kk, M) "module is not finite over the given field"
  F = ambient_free_module(M)
  Mq,_ = sub(F,relations(M))
  Mq_shift,_,_= Oscar.shifted_module(Mq)
  o = negdegrevlex(base_ring(Mq_shift))*lex(ambient_free_module(Mq_shift))
  LMq = leading_module(Mq_shift, o)
  B = _vector_space_basis(kk, quo_object(ambient_free_module(LMq), gens(LMq)), check=false)
  is_empty(B) && return elem_type(M)[]
  _,back_shift = Oscar.base_ring_shifts(R)
  iota = hom(parent(B[1]), M, gens(M), a -> R(back_shift(a)))
  return iota.(B)
end

function _vector_space_basis(kk::Field, M::SubquoModule{T}; check::Bool=true) where {T<:MPolyLocRingElem}
   error("only available in global case and for localization at a point")
end

function _is_finite(kk::Field, M::SubquoModule{T}) where {T<:MPolyLocRingElem{<:Field, <:FieldElem, <:MPolyRing, <:MPolyRingElem, 
                               <:MPolyComplementOfKPointIdeal}}
  @assert kk === coefficient_ring(base_ring(M)) "not implemented for fields other than the `coefficient_ring` of the `base_ring` of the module"
  F = ambient_free_module(M)
  Mq,_ = sub(F,rels(M))
  Mq_shift,_,_ = shifted_module(Mq)
  o = negdegrevlex(base_ring(Mq_shift))*lex(ambient_free_module(Mq_shift))
  LMq = leading_module(Mq_shift, o)
  return _is_finite(kk, quo_object(ambient_free_module(LMq), gens(LMq)))
end

function _is_finite(kk::Field, M::SubquoModule{T}) where {T<:MPolyLocRingElem}
  error("only available in global case and for localization at a point")
end

### functionality for modules over quotients of localized polynomial rings at a point
function _vector_space_basis(kk::Field, M::SubquoModule{T}; check::Bool=true) where {T<:MPolyQuoLocRingElem{<:Field, <:FieldElem, <:MPolyRing, <:MPolyRingElem, <:MPolyComplementOfKPointIdeal}}
 R = base_ring(M)
  @assert kk === coefficient_ring(R) "not implemented for other fields than the coefficients of the underlying polynomial ring"
  @check _is_finite(kk, M) "module is not finite over the given field"
  B = _vector_space_basis(kk, _as_polyloc_module(M), check=false)
  is_empty(B) && return elem_type(M)[]
  iota = _iso_with_polyloc_module(M)
  return iota.(B)
end

function _vector_space_basis(kk::Field, M::SubquoModule{T}; check::Bool=true) where {T<:MPolyQuoLocRingElem}
   error("only available in global case and for localization at a point")
end

function _is_finite(kk::Field, M::SubquoModule{T}) where {T<:MPolyQuoLocRingElem{<:Field, <:FieldElem, <:MPolyRing, <:MPolyRingElem, 
                               <:MPolyComplementOfKPointIdeal}}
  @assert kk === coefficient_ring(base_ring(M)) "not implemented for fields other than the `coefficient_ring` of the `base_ring` of the module"
  return _is_finite(kk, _as_polyloc_module(M))
end

function _is_finite(kk::Field, M::SubquoModule{T}) where {T<:MPolyQuoLocRingElem}
  error("only available in global case and for localization at a point")
end
