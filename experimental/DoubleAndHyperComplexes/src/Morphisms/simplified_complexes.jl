struct SimplifiedChainFactory{ChainType} <: HyperComplexChainFactory{ChainType}
  orig::AbsHyperComplex

  base_change_cache::Dict{Int, Tuple{<:SMat, <:SMat, <:SMat, <:SMat, Vector{Tuple{Int, Int}}}}
  maps_to_original::Dict{Int, <:Map}
  maps_from_original::Dict{Int, <:Map}

  function SimplifiedChainFactory(orig::AbsHyperComplex{ChainType}) where {ChainType}
    d = Dict{Int, Tuple{SMat, SMat}}()
    out_maps = Dict{Int, Map}()
    in_maps = Dict{Int, Map}()
    return new{ChainType}(orig, d, out_maps, in_maps)
  end
end

function (fac::SimplifiedChainFactory)(d::AbsHyperComplex, Ind::Tuple)
  i = first(Ind)
  c = original_complex(fac)

  #              g        f
  #      c_{i-1} ←   c_i  ←  c_{i+1}  --- original complex
  #         ↑↓   ↖g♯ ↑↓   ↙f♭  ↑↓
  #      e_{i-1} ←   e_i  ←  e_{i+1}  --- domains of intermediate `maps_to_original`
  #         ↑↓       ↑↓        ↑↓
  #      d_{i-1} ←   d_i  ←  d_{i+1}  --- simplified complex
  #              g'       f'
  # 
  # We start out by taking the matrix A for f, look for units 
  # and other entries using row- and column operations.
  # If there is already a preliminary simplification of c_{i+1}
  # stored in e_{i+1}, then we replace f by f♯.
  # This leads to a change of basis in c_i and c_{i+1} (resp. e_{i+1}).
  # We store the associated maps in `maps_to_original[i]` 
  # and in `maps_to_original[i+1]` and vice versa for the maps 
  # from the original complex, composing with and replacing the 
  # maps which have already been there, eventually. 
  # This gives presimplified versions for d_i and d_{i+1}. 
  #
  # Then we proceed with the representing matrix B for g. 
  # Since we already computed a presimplified version for c_i, 
  # we may replace g and its representing matrix by g♯ and 
  # even compose with the map to e_{i-1} if applicable. 
  # We do the same thing as for f and set the corresponding 
  # `maps_to/from_original`, eventually replacing them by their 
  # composition with the maps coming out of the operation on B.

  next = (direction(c, 1) == :chain ? i-1 : i+1)
  prev = (direction(c, 1) == :chain ? i+1 : i-1)

  ### Computations for the outgoing map
  if can_compute_index(d, next) && !has_index(d, next)
    # If the first was not true then the next map is simply not there.
    #
    # If the second was not true, then a simplification of the 
    # outgoing map has already been computed as the incoming map 
    # of the second one. We don't need to do it again, then. 
    outgoing = map(c, 1, Ind)
    if haskey(fac.maps_from_original, next)
      outgoing = compose(outgoing, fac.maps_from_original[next])
    end
    if haskey(fac.maps_to_original, i)
      outgoing = compose(fac.maps_to_original[i], outgoing)
    end

    M = domain(outgoing)
    N = codomain(outgoing)

    # Simplify for the outgoing morphism
    A = sparse_matrix(outgoing)
    S, Sinv, T, Tinv, ind = _simplify_matrix!(A)
    fac.base_change_cache[i] = S, Sinv, T, Tinv, ind

    m = nrows(A)
    n = ncols(A)
    I = [i for (i, _) in ind]
    I = [i for i in 1:m if !(i in I)]
    J = [j for (_, j) in ind]
    J = [j for j in 1:n if !(j in J)]

    # Create the maps to the old complex
    img_gens_dom = elem_type(M)[sum(c*M[j] for (j, c) in S[i]; init=zero(M)) for i in I]
    new_dom = _make_free_module(M, img_gens_dom)
    dom_map = hom(new_dom, M, img_gens_dom)

    if haskey(fac.maps_to_original, i)
      # This means that for the next map a partial or 
      # full simplification has already been computed
      fac.maps_to_original[i] = compose(dom_map, fac.maps_to_original[i])
    else
      fac.maps_to_original[i] = dom_map
    end


    img_gens_cod = elem_type(N)[sum(c*N[i] for (i, c) in T[j]; init=zero(N)) for j in J]
    new_cod = _make_free_module(N, img_gens_cod)
    cod_map = hom(new_cod, N, img_gens_cod)

    if haskey(fac.maps_to_original, next)
      fac.maps_to_original[next] = compose(cod_map, fac.maps_to_original[next])
    else
      fac.maps_to_original[next] = cod_map
    end

    # Create the maps from the old complex
    img_gens_dom = elem_type(new_dom)[]
    inv_map_dict = Dict{Int, Int}(I[j] => j for j in 1:length(I))
    for i in 1:m
      w = Sinv[i]
      v = zero(new_dom)
      for (real_j, a) in w
        !haskey(inv_map_dict, real_j) && continue
        j = inv_map_dict[real_j]
        v += a*new_dom[j]
      end

      # v = zero(new_dom)
      # for j in 1:length(I)
      #   success, a = _has_index(w, I[j])
      #   success && (v += a*new_dom[j])
      #   # a = w[I[j]]
      #   # !iszero(a) && (v += a*new_dom[j])
      # end
      push!(img_gens_dom, v)
    end
    dom_map_inv = hom(M, new_dom, img_gens_dom)

    if haskey(fac.maps_from_original, i)
      fac.maps_from_original[i] = compose(fac.maps_from_original[i], dom_map_inv)
    else
      fac.maps_from_original[i] = dom_map_inv
    end


    img_gens_cod = elem_type(new_cod)[]
    inv_map_dict = Dict{Int, Int}(J[j]=>j for j in 1:length(J))
    for i in 1:n
      w = Tinv[i]
      new_entries = Vector{Tuple{Int, elem_type(base_ring(w))}}()
      for (real_j, b) in w
        !haskey(inv_map_dict, real_j) && continue
        j = inv_map_dict[real_j]
        push!(new_entries, (j, b))
      end
      w_new = sparse_row(base_ring(w), new_entries)
      push!(img_gens_cod, FreeModElem(w_new, new_cod))
    end
    cod_map_inv = hom(N, new_cod, img_gens_cod)

    if haskey(fac.maps_from_original, next)
      fac.maps_from_original[next] = compose(fac.maps_from_original[next], cod_map_inv)
    else
      fac.maps_from_original[next] = cod_map_inv
    end
  end

  ### Computations for the incoming map
  if can_compute_index(d, prev) && !has_index(d, prev)
    # If the first was not true, then the incoming map is simply not there.
    # If the second was not true, then a simplification of the incoming 
    # map has already been computed as the outgoing map of the previous 
    # entry and we don't need to do it again.
    incoming = map(c, 1, (prev,))
    if haskey(fac.maps_from_original, i)
      incoming = compose(incoming, fac.maps_from_original[i])
    end
    if haskey(fac.maps_to_original, prev)
      incoming = compose(fac.maps_to_original[prev], incoming)
    end

    M = domain(incoming)
    N = codomain(incoming)

    # Simplify for the incoming morphism
    A = sparse_matrix(incoming)
    S, Sinv, T, Tinv, ind = _simplify_matrix!(A)

    m = nrows(A)
    n = ncols(A)
    I = [i for (i, _) in ind]
    I = [i for i in 1:m if !(i in I)]
    J = [j for (_, j) in ind]
    J = [j for j in 1:n if !(j in J)]

    # Create the maps to the old complex
    img_gens_dom = elem_type(M)[sum(c*M[j] for (j, c) in S[i]; init=zero(M)) for i in I]
    new_dom = _make_free_module(M, img_gens_dom)
    dom_map = hom(new_dom, M, img_gens_dom)

    if haskey(fac.maps_to_original, prev)
      fac.maps_to_original[prev] = compose(dom_map, fac.maps_to_original[prev])
    else
      fac.maps_to_original[prev] = dom_map
    end


    img_gens_cod = elem_type(N)[sum(c*N[i] for (i, c) in T[j]; init=zero(N)) for j in J]
    new_cod = _make_free_module(N, img_gens_cod)
    cod_map = hom(new_cod, N, img_gens_cod)

    if haskey(fac.maps_to_original, i)
      fac.maps_to_original[i] = compose(cod_map, fac.maps_to_original[i])
    else 
      fac.maps_to_original[i] = cod_map
    end

    # Create the maps from the old complex
    img_gens_dom = elem_type(new_dom)[]
    for i in 1:m
      w = Sinv[i]
      v = zero(new_dom)
      for j in 1:length(I)
        a = w[I[j]]
        !iszero(a) && (v += a*new_dom[j])
      end
      push!(img_gens_dom, v)
    end
    dom_map_inv = hom(M, new_dom, img_gens_dom)

    if haskey(fac.maps_from_original, prev)
      fac.maps_from_original[prev] = compose(fac.maps_from_original[prev], dom_map_inv)
    else
      fac.maps_from_original[prev] = dom_map_inv
    end


    img_gens_cod = elem_type(new_cod)[]
    for i in 1:n
      w = Tinv[i]
      new_entries = Vector{Tuple{Int, elem_type(base_ring(w))}}()
      for (real_j, b) in w
        j = findfirst(k->k==real_j, J)
        j === nothing && continue
        push!(new_entries, (j, b))
      end
      w_new = sparse_row(base_ring(w), new_entries)
      push!(img_gens_cod, FreeModElem(w_new, new_cod))
    end
    cod_map_inv = hom(N, new_cod, img_gens_cod)

    if haskey(fac.maps_from_original, i)
      fac.maps_from_original[i] = compose(fac.maps_from_original[i], cod_map_inv)
    else
      fac.maps_from_original[i] = cod_map_inv
    end
  end

  return domain(fac.maps_to_original[i])
end

function can_compute(fac::SimplifiedChainFactory, c::AbsHyperComplex, I::Tuple)
  i = first(I)
  can_compute_index(original_complex(fac), I) || return false
  I_p = (i+1,)
  I_m = (i-1,)
  if can_compute_index(original_complex(fac), I_p) 
    if direction(c, 1) == :chain
      can_compute_map(original_complex(fac), 1, I_p) || return false
    else
      can_compute_map(original_complex(fac), 1, I) || return false
    end
  end

  if can_compute_index(original_complex(fac), I_m)
    if direction(c, 1) == :chain
      can_compute_map(original_complex(fac), 1, I) || return false
    else
      can_compute_map(original_complex(fac), 1, I_m) || return false
    end
  end
  return true
end

maps_to_original(fac::SimplifiedChainFactory) = fac.maps_to_original
maps_from_original(fac::SimplifiedChainFactory) = fac.maps_from_original
original_complex(fac::SimplifiedChainFactory) = fac.orig


### Simplified maps 

struct SimplifiedMapFactory{MorphismType} <: HyperComplexMapFactory{MorphismType}
  orig::AbsHyperComplex

  function SimplifiedMapFactory(orig::AbsHyperComplex{<:Any, MorphismType}) where {MorphismType}
    return new{MorphismType}(orig)
  end
end

function (fac::SimplifiedMapFactory)(c::AbsHyperComplex, p::Int, I::Tuple)
  # fill the cache
  c[I]
  i = first(I)
  next = i + (direction(c, 1) == :chain ? -1 : 1)
  c[next]
  chain_fac = chain_factory(c)
  dom_out = maps_to_original(chain_fac)[i]
  @assert domain(dom_out) === c[I]
  cod_in = maps_from_original(chain_fac)[next]
  return compose(compose(dom_out, map(fac.orig, 1, I)), cod_in)
end

function can_compute(fac::SimplifiedMapFactory, c::AbsHyperComplex, p::Int, I::Tuple)
  return can_compute_map(original_complex(chain_factory(c)), p, I)
end

### Factories for the base change morphisms from and to the original complex
struct BaseChangeToOriginalFactory{MorphismType} <: HyperComplexMorphismFactory{MorphismType} end

function (fac::BaseChangeToOriginalFactory)(phi::AbsHyperComplexMorphism, I::Tuple)
  d = domain(phi)
  chain_fac = chain_factory(underlying_complex(d))
  d[I] # fill the cache
  return maps_to_original(chain_fac)[first(I)]
end

function can_compute(fac::BaseChangeToOriginalFactory, phi::AbsHyperComplexMorphism, I::Tuple)
  d = domain(phi)
  return can_compute_index(d, I)
end

struct BaseChangeFromOriginalFactory{MorphismType} <: HyperComplexMorphismFactory{MorphismType} end

function (fac::BaseChangeFromOriginalFactory)(phi::AbsHyperComplexMorphism, I::Tuple)
  d = codomain(phi)
  chain_fac = chain_factory(underlying_complex(d))
  d[I] # fill the cache
  return maps_from_original(chain_fac)[first(I)]
end

function can_compute(fac::BaseChangeFromOriginalFactory, phi::AbsHyperComplexMorphism, I::Tuple)
  d = domain(phi)
  return can_compute_index(d, I)
end


function simplify(c::AbsHyperComplex{ChainType, MorphismType}) where {ChainType, MorphismType}
  @assert dim(c) == 1 "complex must be one-dimensional"
  chain_fac = SimplifiedChainFactory(c)
  mor_fac = SimplifiedMapFactory(c)
  upper_bounds = [has_upper_bound(c, 1) ? upper_bound(c, 1) : nothing]
  lower_bounds = [has_lower_bound(c, 1) ? lower_bound(c, 1) : nothing]
  directions = [direction(c, 1)]
  internal_complex = HyperComplex(1, chain_fac, mor_fac, 
                                  directions, upper_bounds=upper_bounds, 
                                  lower_bounds=lower_bounds
                                 )

  result = SimplifiedComplex(internal_complex, c)
  result.map_to_original = HyperComplexMorphism(result, c, 
                                                BaseChangeToOriginalFactory{MorphismType}()
                                               )
  result.map_from_original = HyperComplexMorphism(c, result,
                                                  BaseChangeFromOriginalFactory{MorphismType}()
                                                 )
  return result
end

### Helper functions
function _make_free_module(M::ModuleFP, g::Vector{T}) where {T<:ModuleFPElem}
  if is_graded(M)
    w = degree.(g)
    return graded_free_module(base_ring(M), w)
  else
    return FreeMod(base_ring(M), length(g))
  end
end

@doc raw""" 
    _simplify_matrix!(A::SMat; find_pivot=nothing)

For a matrix `A` this looks for units in `A` and uses row- and column-
operations to eliminate the entries in `A` above, below, to the left, 
and to the right of that unit. It returns a quintuple 
`(S, Sinv, T, Tinv, ind)` where `S` and `Sinv` (resp. `T` and `Tinv`) 
are mutual inverse matrices such that `Sinv*A*T` is the original matrix 
put in and `ind` is a `Vector` of pairs `(i, j)` which are the indices 
of the units in `A` which have been used for elimination.
The optional argument `find_pivot` must either be a function `f`, or an 
object with overloaded call syntax, which allows for the call 
`f(A)` to return a pair of indices `(i, j)` such that `A[i, j]` is a 
unit in the `base_ring` of `A` to be used as the next pivot element, 
or `nothing` if no suitable pivot was found.
"""
function _simplify_matrix!(A::SMat; find_pivot=nothing)
  R = base_ring(A)
  m = nrows(A)
  n = ncols(A)

  done = Vector{Tuple{Int, Int}}()
  done_rows = Vector{Int}()
  done_columns = Vector{Int}()

  # Initialize the base change matrices in domain and codomain
  S = sparse_matrix(R, m, m)
  for i in 1:m
    S[i] = sparse_row(R, [(i, one(R))])
  end
  Sinv_transp = deepcopy(S)

  T = sparse_matrix(R, n, n)
  for i in 1:n
    T[i] = sparse_row(R, [(i, one(R))])
  end
  Tinv_transp = deepcopy(T)

  br = 0
  found_unit = true
  i_start = 0
  while found_unit
    # Find any unit entry
    p = q = 0

    if find_pivot === nothing
      found_unit = false
      for i in i_start+1:nrows(A)
        i in done_rows && continue
        v = A[i]
        found_unit && break
        for (j, c) in v
          j in done_columns && continue
          if is_unit(c) 
            p = i
            q = j
            found_unit = true
            break
          end
        end
        found_unit && break
      end
      if !found_unit
        for i in 1:i_start
          i in done_rows && continue
          v = A[i]
          found_unit && break
          for (j, c) in v
            j in done_columns && continue
            if is_unit(c) 
              p = i
              q = j
              found_unit = true
              break
            end
          end
          found_unit && break
        end
      end
    else
      res = find_pivot(A)
      if res === nothing
        found_unit = false
      else
        (p, q) = res::Tuple{Int, Int}
        found_unit = true
      end
    end

    !found_unit && break

    i_start = p
    if i_start == nrows(A) 
      i_start = 0
    end

    push!(done, (p, q))
    push!(done_rows, p)
    push!(done_columns, q)
    u = _get_index_fast(A, p, q)
    #u = A[p, q]

    # We found a unit in the entry (p, q) of the matrix which has not yet 
    # been treated.

    a_row = deepcopy(A[p])
    a_row_del = a_row - sparse_row(R, [(q, u)])

    col_entries = Vector{Tuple{Int, elem_type(R)}}()
    for i in 1:m
      # c = A[i, q]
      # !iszero(c) && push!(col_entries, (i, c))
      success, c = _has_index(A, i, q)
      success && push!(col_entries, (i, c::elem_type(R)))
    end
    a_col = sparse_row(R, col_entries)
    a_col_del = a_col - sparse_row(R, [(p, u)])

    uinv = inv(u)

    # clear the q-th column
    for (i, b) in a_col_del
      #A[i] = A[i] - uinv*b*a_row # original operation, replaced by inplace arithmetic below
      Hecke.add_scaled_row!(a_row, A[i], -uinv*b)
    end

    A[p] = sparse_row(R, [(q, u)])

    # Adjust S
    v = S[p]
    for (i, b) in a_col_del
      #S[i] = S[i] - uinv*b*v # original operation, replaced by inplace arithmetic below
      Hecke.add_scaled_row!(v, S[i], -uinv*b)
    end

    # Adjust Sinv_transp
    inc_row = sum(a*Sinv_transp[i] for (i, a) in a_col_del; init=sparse_row(R))
    Hecke.add_scaled_row!(inc_row, Sinv_transp[p], uinv)

    # Adjust T
    #T[q] = T[q] + sum(uinv*a*T[i] for (i, a) in a_row_del; init=sparse_row(R))
    Hecke.add_scaled_row!(sum(a*T[i] for (i, a) in a_row_del if !iszero(T[i]); init=sparse_row(R)), T[q], uinv)

    # Adjust Tinv_transp
    v = deepcopy(Tinv_transp[q])
    for (j, b) in a_row_del
      Hecke.add_scaled_row!(v, Tinv_transp[j], -uinv*b)
    end

    # Kept here for debugging. These must be true:
    # @assert isone(matrix(T)*matrix(transpose(Tinv_transp)))
    # @assert isone(matrix(S)*matrix(transpose(Sinv_transp)))
    # @assert matrix(S)*matrix(B)*transpose(matrix(Tinv_transp)) == matrix(A)
  end
 
  return S, transpose(Sinv_transp), T, transpose(Tinv_transp), done
end


function sparse_matrix(phi::FreeModuleHom{FreeMod{T}, FreeMod{T}, Nothing}) where {T}
  V = domain(phi)
  W = codomain(phi)
  kk = base_ring(V)
  m = ngens(V)
  n = ngens(W)
  result = sparse_matrix(kk, m, n)
  for i in 1:m
    result[i] = phi(V[i]).coords
  end
  return result
end

@doc raw"""
    _has_index(a::SRow, i::Int)

Find out whether `a` has a non-zero entry in position `i` and return `true, a[i]` 
if this is the case and `false, nothing` otherwise.
"""
function _has_index(a::SRow, i::Int)
  isempty(a.pos) && return false, nothing
  k0 = 1
  k1 = length(a.pos)
  while k0 <= k1
    k = div(k0 + k1, 2)
    j = a.pos[k]
    if j > i
      k1 = k - 1
    elseif j < i
      k0 = k + 1
    else 
      return true, a.values[k]
    end
  end
  return false, nothing
end

# This addresses the TODO in Heckes src/Sparse/Row.jl:290
function _get_index_fast(a::SRow, i::Int)
  isempty(a.pos) && return zero(base_ring(a))
  k0 = 1
  k1 = length(a.pos)
  while k0 <= k1
    k = div(k0 + k1, 2)
    j = a.pos[k]
    if j > i
      k1 = k - 1
    elseif j < i
      k0 = k + 1
    else
      return a.values[k]
    end
  end
  return zero(base_ring(a))
end

function _get_index_fast(a::SMat, i::Int, j::Int)
  return _get_index_fast(a[i], j)
end

function _has_index(a::SMat, i::Int, j::Int)
  return _has_index(a[i], j)
end

### Taylormade functionality for modules
# Provided as a hotfix for issue #3108.
function _alt_simplify(M::SubquoModule)
  res, aug = free_resolution(SimpleFreeResolution, M)
  simp = simplify(res)
  simp_to_orig = map_to_original_complex(simp)
  orig_to_simp = map_from_original_complex(simp)
  result, Z0_to_result = homology(simp, 0)
  Z0, inc_Z0 = kernel(simp, 0)

  result_to_M = hom(result, M, 
                    elem_type(M)[aug[0](simp_to_orig[0](inc_Z0(preimage(Z0_to_result, x)))) for x in gens(result)])
  M_to_result = hom(M, result,
                    elem_type(result)[Z0_to_result(preimage(inc_Z0, orig_to_simp[0](preimage(aug[0], y)))) for y in gens(M)])
  return result, M_to_result, result_to_M
end

# Some special shortcuts
function boundary(c::SimplifiedComplex{ChainType}, p::Int, i::Tuple) where {ChainType<:ModuleFP}
  # try to find the boundary already computed
  und = underlying_complex(c)::HyperComplex
  if !isdefined(und, :boundary_cache) 
    und.boundary_cache = Dict{Tuple{Tuple, Int}, Map}()
  end
  if haskey(und.boundary_cache, (i, p))
    inc = und.boundary_cache[(i, p)]
    return domain(inc), inc
  end

  # compute the boundary from scratch
  orig_boundary, orig_inc = boundary(c.original_complex, p, i)
  from_orig = map_from_original_complex(c)[i]
  result, inc = sub(c[i], from_orig.(orig_inc.(gens(orig_boundary))))

  # Cache the result
  und.boundary_cache[(i, p)] = inc

  return result, inc
end

function kernel(simp::SimplifiedComplex{ChainType}, p::Int, i::Tuple) where {ChainType<:ModuleFP}
  c = underlying_complex(simp)::HyperComplex

  if !isdefined(c, :kernel_cache) 
    c.kernel_cache = Dict{Tuple{Tuple, Int}, Map}()
  end
  if haskey(c.kernel_cache, (i, p))
    inc = c.kernel_cache[(i, p)]
    return domain(inc), inc
  end

  if !can_compute_map(simp, p, i)
    M = simp[i]
    K, inc = sub(simp[i], gens(simp[i]))
    c.kernel_cache[(i, p)] = inc
    return K, inc
  end

  psi = map_to_original(simp)[i]
  phi = map(c, p, i)
  K, inc = kernel(compose(psi, phi))
  c.kernel_cache[(i, p)] = inc
  return K, inc
end

function homology(simp::SimplifiedComplex{ChainType}, p::Int, i::Tuple) where {ChainType <: ModuleFP}
  c = underlying_complex(simp)::HyperComplex
  
  if !isdefined(c, :homology_cache) 
    c.homology_cache = Dict{Tuple{Tuple, Int}, Map}()
  end
  if haskey(c.homology_cache, (i, p))
    pr = c.homology_cache[(i, p)]
    return codomain(pr), pr
  end

  H, pr = quo(kernel(simp, p, i)[1], boundary(simp, p, i)[1])
  c.homology_cache[(i, p)] = pr
  return H, pr
end

homology(simp::SimplifiedComplex{ChainType}, i::Int) where {ChainType <: ModuleFP} = homology(simp, 1, (i,))
