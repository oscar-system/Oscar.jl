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
function tensor_product(G::ModuleFP...; task::Symbol = :none, minimal::Bool=false)
  is_empty(G) && error("list of modules must not be empty")
  R = base_ring(first(G))
  @assert all(base_ring(M) === R for M in G) "modules must be defined over the same ring"
  # The presentations need to map generators to generators! 
  # We can not match elements, unless this is the case.
  presentations = [presentation(M; minimal) for M in G]
  mats = [sparse_matrix(map(c, 1)) for c in presentations] 
  @assert all(ncols(B) == ngens(M) for (M, B) in zip(G, mats)) "presentations do not implement a 1:1 correspondence for the generators"
  multi_inds = collect(AbstractAlgebra.ProductIterator([1:ngens(M) for M in reverse(G)]))
  multi_inds = reverse.(multi_inds)
  ind_pos = Dict(i=>j for (j, i) in enumerate(multi_inds))
  ids = [_sparse_identity_matrix(R, ngens(M)) for M in G]
  A = sparse_matrix(R, 0, prod(ngens(M) for M in G; init=1))
  for i in 1:length(ids)
    id = ids[i]
    ids[i] = mats[i]
    B = _kronecker_product(ids)
    vcat!(A, B)
    ids[i] = id
  end
  F = FreeMod(R, ncols(A))
  @assert ngens(F) == length(multi_inds)

  # In the graded case we just swap over the degrees. 
  if all(is_graded, G)
    all_degs = [degrees_of_generators(M) for M in G]
    GG = grading_group(R)
    F.d = [sum(all_degs[k][j] for (k, j) in enumerate(ind); init=zero(GG)) for ind in multi_inds]
  end

  I = [F(v) for v in A]
  result = SubquoModule(F, gens(F), I)

  # Functions for tensor multiplication and decomposition
  function pure_helper(inds::Vector{Int}, rem::Vector)
    is_one(length(rem)) && return sum(c*F[ind_pos[vcat(inds, [i])]] for (i, c) in coordinates(only(rem)); init=zero(F))
    v = first(rem)
    return sum(c*pure_helper(vcat(inds, [i]), rem[2:end]) for (i, c) in coordinates(v); init=zero(F))
  end
  function pure(t::ModuleFPElem...)
    return result(pure_helper(Int[], collect(t)))
  end
  function pure(t::Tuple)
    return pure(t...)
  end

  function decomp(v::SubquoModuleElem)
    i, c = only(coordinates(v))
    is_one(c) || error("element must be a module generator")
    return Tuple([M[j] for (j, M) in zip(multi_inds[i], G)])
  end

  # set attributes as needed
  set_attribute!(result, :tensor_pure_function=>pure, :tensor_generator_decompose_function=>decomp)
  set_attribute!(result, :show => Hecke.show_tensor_product, :tensor_product => G)

  task == :none && return result
  return result, MapFromFunc(Hecke.TupleParent(Tuple([zero(g) for g = G])), result, pure, decomp)
end

# recursive method with bisection for many matrices
function _kronecker_product(mats::Vector)
  is_one(length(mats)) && return only(mats)
  length(mats) == 2 && return _kronecker_product(mats[1], mats[2])
  n = length(mats)
  h = div(n, 2)
  return _kronecker_product(_kronecker_product(mats[1:h]), _kronecker_product(mats[h+1:end]))
end

# Probably the same as in AA. However, we rely on a compatibility with the 
# other iterations here, so we need to make sure this does the right thing. 
function _kronecker_product(A::SMat, B::SMat)
  R = base_ring(A)
  m1, n1 = size(A)
  m2, n2 = size(B)
  C = sparse_matrix(R, 0, n1*n2)
  for v1 in A
    for v2 in B
      # expanded version of 
      # push!(C, sparse_row(R, [((j1-1)*n2+j2, c1*c2) for (j1, c1) in v1 for (j2, c2) in v2]; sort=false))
      # to avoid allocation of extra vectors
      new_pos = Int[(j1-1)*n2+j2 for j1 in v1.pos for j2 in v2.pos]
      new_vals = elem_type(R)[c1*c2 for c1 in v1.values for c2 in v2.values]
      push!(C, sparse_row(R, new_pos, new_vals; sort=false))
    end
  end
  return C
end

# Helper function to create the identity matrix. Does this exist somewhere already? 
function _sparse_identity_matrix(R::Ring, n::Int)
  A = sparse_matrix(R, 0, n)
  for i in 1:n
    push!(A, sparse_row(R, i, one(R); check=false))
  end
  return A
end

