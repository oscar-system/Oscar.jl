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
  t = [[x] for x = 1:ngens(G[1])]
  for H = G[2:end]
    t = [push!(deepcopy(x), y) for x = t  for y = 1:ngens(H)]
  end

  F = FreeMod(G[1].R, prod([rank(g) for g in G]))
  F.S = function _get_tensor_symbols()
    return Symbol.([join(["$(symbol(G[k], i[length(G)-k+1]))" for k in 1:length(G)], " \\otimes ") for i in AbstractAlgebra.ProductIterator([1:ngens(g) for g in reverse(G)])])
  end

  set_attribute!(F, :show => Hecke.show_tensor_product, :tensor_product => G)

  function pure(g::FreeModElem...)
    @assert length(g) == length(G)
    @assert all(i -> parent(g[i]) === G[i], 1:length(G))
    z = [[x] for x = coordinates(g[1]).pos]
    zz = coordinates(g[1]).values
    for h = g[2:end]
      zzz = Vector{Int}[]
      zzzz = elem_type(F.R)[]
      for i = 1:length(z)
        for (p, v) in coordinates(h)
          push!(zzz, push!(deepcopy(z[i]), p))
          push!(zzzz, zz[i]*v)
        end
      end
      z = zzz
      zz = zzzz
    end
    indices = Vector{Int}([findfirst(==(y), t) for y = z])
    return FreeModElem(sparse_row(F.R, indices, zz), F)
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
    return Tuple(gen(G[i], t[c.pos[1]][i]) for i = 1:length(G))
  end

  set_attribute!(F, :tensor_pure_function => pure, :tensor_generator_decompose_function => inv_pure)

  if all(is_graded, G)
    tensor_degrees = [sum(G[i].d[tplidx[i]] for i in 1:length(G)) for tplidx in t] 
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

  pres = [presentation(M; minimal=minimal) for M in G]
  d1s = [map(p, 1) for p in pres]
  augs = [map(p, 0) for p in pres]
  F0s = [domain(a) for a in augs]
  mats = [sparse_matrix(d) for d in d1s]
  if !minimal
    @assert all(ncols(mats[i]) == ngens(G[i]) for i in 1:length(G)) "presentations do not implement a 1:1 correspondence for the generators"
  end

  F = tensor_product(F0s...; task=:none)

  ids = [_sparse_identity_matrix(R, ngens(F0)) for F0 in F0s]
  A = sparse_matrix(R, 0, ngens(F))
  for i in 1:length(G)
    saved = ids[i]
    ids[i] = mats[i]
    block = _kronecker_product(ids)
    vcat!(A, block)
    ids[i] = saved
  end

  I = [F(v) for v in A]
  result = SubquoModule(F, gens(F), I)

  free_pure = tensor_pure_function(F)
  free_decomp = tensor_generator_decompose_function(F)

  function pure(tuple_elems::ModuleFPElem...)
    @assert length(tuple_elems) == length(G)
    @assert all(i -> parent(tuple_elems[i]) === G[i], 1:length(G))
    pre = [preimage(augs[i], tuple_elems[i]) for i in 1:length(G)]
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
    return Tuple(augs[j](dec[j]) for j in 1:length(G))
  end

  set_attribute!(result, :tensor_pure_function => pure, :tensor_generator_decompose_function => decompose_generator)
  set_attribute!(result, :show => Hecke.show_tensor_product, :tensor_product => G)
  @assert _is_tensor_product(result)[1]

  if task == :none
    return result
  end

  return result, MapFromFunc(Hecke.TupleParent(Tuple([zero(g) for g = G])), result, pure, decompose_generator)
end

function _kronecker_product(mats::Vector)
  @assert !isempty(mats)
  return _kronecker_product(mats, 1, length(mats))
end

function _kronecker_product(mats::Vector, l::Int, r::Int)
  l == r && return mats[l]
  if l + 1 == r
    return _kronecker_product(mats[l], mats[r])
  end
  m = (l + r) >>> 1
  return _kronecker_product(_kronecker_product(mats, l, m), _kronecker_product(mats, m + 1, r))
end

function _kronecker_product(A, B)
  R = base_ring(A)
  m1, n1 = size(A)
  m2, n2 = size(B)
  C = sparse_matrix(R, 0, n1 * n2)
  T = elem_type(R)
  for v1 in A
    for v2 in B
      push!(C, sparse_row(R, [((j1-1)*n1+j2, c1*c2) for (j1, c1) in v1 for (j2, c2) in v2]))
    end
  end
  return C

  for i in 1:m1
    r1 = A[i]
    p1 = r1.pos
    v1 = r1.values
    nnz1 = length(p1)
    for j in 1:m2
      r2 = B[j]
      p2 = r2.pos
      v2 = r2.values
      nnz2 = length(p2)
      if nnz1 == 0 || nnz2 == 0
        push!(C, sparse_row(R))
        continue
      end
      pos = Vector{Int}(undef, nnz1 * nnz2)
      vals = Vector{T}(undef, nnz1 * nnz2)
      k = 0
      for a in 1:nnz1
        base = (p1[a] - 1) * n2
        c1 = v1[a]
        for b in 1:nnz2
          prod = c1 * v2[b]
          iszero(prod) && continue
          k += 1
          pos[k] = base + p2[b]
          vals[k] = prod
        end
      end
      if k == 0
        push!(C, sparse_row(R))
      else
        resize!(pos, k)
        resize!(vals, k)
        push!(C, sparse_row(R, pos, vals))
      end
    end
  end
  return C
end

function _sparse_identity_matrix(R, n::Int)
  A = sparse_matrix(R, 0, n)
  for i in 1:n
    push!(A, sparse_row(R, [i], [one(R)]))
  end
  return A
end
