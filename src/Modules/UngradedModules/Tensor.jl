##################################################
# Tensor
##################################################

@doc raw"""
    tensor_product(F::FreeMod...; task::Symbol = :none)

Given a collection of free modules, say, $F_1, \dots, F_n$ over a ring $R$, return $F_1\otimes_R \cdots \otimes_R F_n$.


If `task = :map`, additionally return the map which sends a tuple $(f_1,\dots, f_n)$ of elements $f_i\in F_i$ to the pure tensor $f_1\otimes\dots\otimes f_n$.
"""
function tensor_product(G::FreeMod...; task::Symbol = :none)
  gs = [is_graded(g) for g in G]
  if !all(gs) && !all(!x for x in gs)
    error("All factors must either be graded or all must be ungraded.")
  end
  n = length(G)
  sizes = [ngens(g) for g in G]
  strides = Vector{Int}(undef, n)
  strides[n] = 1
  for k in (n - 1):-1:1
    strides[k] = strides[k + 1] * sizes[k + 1]
  end

  rankF = prod(sizes)
  F = FreeMod(G[1].R, rankF)
  F.S = function _get_tensor_symbols()
    return Symbol.([join(["$(symbol(G[k], i[length(G)-k+1]))" for k in 1:length(G)], " \\otimes ") for i in AbstractAlgebra.ProductIterator([1:ngens(g) for g in reverse(G)])])
  end

  set_attribute!(F, :show => Hecke.show_tensor_product, :tensor_product => G)

  function pure(g::FreeModElem...)
    @assert length(g) == n
    @assert all(i -> parent(g[i]) === G[i], 1:n)
    c = coordinates(g[1])
    if length(c.pos) == 0
      return zero(F)
    end
    indices = Vector{Int}(undef, length(c.pos))
    coeffs = Vector{elem_type(F.R)}(undef, length(c.pos))
    for a in 1:length(c.pos)
      indices[a] = 1 + (c.pos[a] - 1) * strides[1]
      coeffs[a] = c.values[a]
    end
    nz = length(indices)
    for k in 2:n
      ck = coordinates(g[k])
      if length(ck.pos) == 0 || nz == 0
        return zero(F)
      end
      new_indices = Vector{Int}(undef, nz * length(ck.pos))
      new_coeffs = Vector{elem_type(F.R)}(undef, nz * length(ck.pos))
      t = 0
      for a in 1:nz
        idx_base = indices[a]
        coeff_base = coeffs[a]
        for b in 1:length(ck.pos)
          prod = coeff_base * ck.values[b]
          iszero(prod) && continue
          t += 1
          new_indices[t] = idx_base + (ck.pos[b] - 1) * strides[k]
          new_coeffs[t] = prod
        end
      end
      if t == 0
        return zero(F)
      end
      resize!(new_indices, t)
      resize!(new_coeffs, t)
      indices = new_indices
      coeffs = new_coeffs
      nz = t
    end
    return FreeModElem(sparse_row(F.R, indices, coeffs), F)
  end
  function pure(T::Tuple)
    return pure(T...)
  end
  function inv_pure(e::FreeModElem)
    c = coordinates(e)
    if length(c.pos) == 0
      return Tuple(zero(g) for g = G)
    end
    @assert length(c.pos) == 1
    @assert isone(c.values[1])
    q = c.pos[1] - 1
    return Tuple(gen(G[i], (q ÷ strides[i]) % sizes[i] + 1) for i = 1:n)
  end

  set_attribute!(F, :tensor_pure_function => pure, :tensor_generator_decompose_function => inv_pure)

  if all(is_graded, G)
    GG = grading_group(F.R)
    tensor_degrees = Vector{elem_type(GG)}(undef, rankF)
    for q in 0:(rankF - 1)
      d = zero(GG)
      for i in 1:n
        d += G[i].d[(q ÷ strides[i]) % sizes[i] + 1]
      end
      tensor_degrees[q + 1] = d
    end
    F.d = tensor_degrees
  end

  @assert _is_tensor_product(F)[1]
  if task == :none
    return F
  end

  return F, MapFromFunc(Hecke.TupleParent(Tuple([zero(g) for g = G])), F, pure, inv_pure)
end

⊗(G::ModuleFP...) = tensor_product(G..., task = :none)

@doc raw"""
    tensor_product(M::ModuleFP...; task::Symbol = :none)

Given a collection of modules, say, $M_1, \dots, M_n$ over a ring $R$, return $M_1\otimes_R \cdots \otimes_R M_n$.

If `task = :map`, additionally return the map which sends a tuple $(m_1,\dots, m_n)$ of elements $m_i\in M_i$ to the pure tensor $m_1\otimes\dots\otimes m_n$.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> F = free_module(R, 1);

julia> A = R[x; y];

julia> B = R[x^2; y^3; z^4];

julia> M = SubquoModule(F, A, B);

julia> gens(M)
2-element Vector{SubquoModuleElem{QQMPolyRingElem}}:
 x*e[1]
 y*e[1]

julia> T, t = tensor_product(M, M; task = :map);

julia> gens(T)
4-element Vector{SubquoModuleElem{QQMPolyRingElem}}:
 (e[1] \otimes e[1])
 (e[1] \otimes e[2])
 (e[2] \otimes e[1])
 (e[2] \otimes e[2])

julia> domain(t)
parent of tuples of type Tuple{SubquoModuleElem{QQMPolyRingElem}, SubquoModuleElem{QQMPolyRingElem}}

julia> t((M[1], M[2]))
(e[1] \otimes e[2])
```
"""
function tensor_product(G::ModuleFP...; task::Symbol = :none, minimal::Bool = false)
  @assert length(G) > 0 "list of modules must not be empty"
  gs = [is_graded(g) for g in G]
  if !all(gs) && !all(!x for x in gs)
    error("All factors must either be graded or all must be ungraded.")
  end
  R = base_ring(G[1])
  @assert all(base_ring(M) === R for M in G) "modules must be defined over the same ring"

  n = length(G)
  cache = IdDict{Any, Tuple{Any, Any, Any}}()
  augs = Vector{Any}(undef, n)
  F0s = Vector{Any}(undef, n)
  mats = Vector{Any}(undef, n)
  for i in 1:n
    M = G[i]
    data = get!(cache, M) do
      p = presentation(M; minimal=minimal)
      d1 = map(p, 1)
      aug = map(p, 0)
      F0 = domain(aug)
      mat = _sparse_matrix_cached(d1)
      (aug, F0, mat)
    end
    augs[i] = data[1]
    F0s[i] = data[2]
    mats[i] = data[3]
  end
  if !minimal
    @assert all(ncols(mats[i]) == ngens(G[i]) for i in 1:n) "presentations do not implement a 1:1 correspondence for the generators"
  end

  F = tensor_product(F0s...; task=:none)

  sizes = [ngens(F0) for F0 in F0s]
  strides = Vector{Int}(undef, n)
  strides[n] = 1
  for k in (n - 1):-1:1
    strides[k] = strides[k + 1] * sizes[k + 1]
  end

  rankF = rank(F)
  I = elem_type(F)[]
  if rankF != 0
    sizehint!(I, sum(nrows(mats[i]) * div(rankF, sizes[i]) for i in 1:n))
    for i in 1:n
      _append_tensor_relations!(I, F, mats[i], sizes, strides, i)
    end
  end

  result = SubquoModule(F, gens(F), I)

  free_pure = tensor_pure_function(F)
  free_decomp = tensor_generator_decompose_function(F)

  function pure(tuple_elems::ModuleFPElem...)
    @assert length(tuple_elems) == n
    @assert all(i -> parent(tuple_elems[i]) === G[i], 1:n)
    rankF == 0 && return zero(result)
    pre = ntuple(i -> preimage(augs[i], tuple_elems[i]), n)
    ww = free_pure(pre...)
    return result(coordinates(ww))
  end
  pure(T::Tuple) = pure(T...)

  function decompose_generator(v::SubquoModuleElem)
    c = coordinates(v)
    if length(c.pos) == 0
      return Tuple(zero(M) for M in G)
    end
    @assert length(c.pos) == 1
    @assert isone(c.values[1])
    i = c.pos[1]
    dec = free_decomp(gen(F, i))
    return Tuple(augs[j](dec[j]) for j in 1:n)
  end

  set_attribute!(result, :tensor_pure_function => pure, :tensor_generator_decompose_function => decompose_generator)
  set_attribute!(result, :show => Hecke.show_tensor_product, :tensor_product => G)
  @assert _is_tensor_product(result)[1]

  if task == :none
    return result
  end

  return result, MapFromFunc(Hecke.TupleParent(Tuple([zero(g) for g = G])), result, pure, decompose_generator)
end

@attr Any function _sparse_matrix_cached(d::ModuleFPHom)
  return sparse_matrix(d)
end

function _append_tensor_relations!(rels::Vector, F::FreeMod, mat, sizes::Vector{Int}, strides::Vector{Int}, i::Int)
  n = length(sizes)
  nrows(mat) == 0 && return rels
  stride_i = strides[i]
  idx = ones(Int, n)
  base = 1
  while true
    for r in mat
      p = r.pos
      if length(p) == 0
        continue
      end
      pos = Vector{Int}(undef, length(p))
      for a in 1:length(p)
        pos[a] = base + (p[a] - 1) * stride_i
      end
      push!(rels, FreeModElem(sparse_row(base_ring(F), pos, copy(r.values); sort=false), F))
    end
    k = n
    while k >= 1
      if k == i
        k -= 1
        continue
      end
      idx[k] += 1
      base += strides[k]
      if idx[k] <= sizes[k]
        break
      end
      idx[k] = 1
      base -= sizes[k] * strides[k]
      k -= 1
    end
    k == 0 && break
  end
  return rels
end
