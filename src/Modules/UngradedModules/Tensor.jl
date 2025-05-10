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
  s = G[1].S
  t = [[x] for x = 1:ngens(G[1])]
  for H = G[2:end]
    s = [Symbol("$x \\otimes $y") for x = s  for y = H.S]
    t = [push!(deepcopy(x), y) for x = t  for y = 1:ngens(H)]
  end

  F = FreeMod(G[1].R, prod([rank(g) for g in G]))
  F.S = s
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
 x^2*e[1] \otimes e[1]
 x*y*e[1] \otimes e[1]
 x*y*e[1] \otimes e[1]
 y^2*e[1] \otimes e[1]

julia> domain(t)
parent of tuples of type Tuple{SubquoModuleElem{QQMPolyRingElem}, SubquoModuleElem{QQMPolyRingElem}}

julia> t((M[1], M[2]))
x*y*e[1] \otimes e[1]
```
"""
function tensor_product(G::ModuleFP...; 
    ambient_tensor_product::FreeMod=tensor_product((ambient_free_module(x) for x in G)..., task = :none),
    task::Symbol = :none
  )
  F = ambient_tensor_product
  mF = tensor_pure_function(ambient_tensor_product)
  # We want to store a dict where the keys are tuples of indices and the values
  # are the corresponding pure vectors (i.e. a tuple (2,1,5) represents the 
  # 2nd, 1st and 5th generator of the 1st, 2nd and 3rd module, which we are 
  # tensoring. The corresponding value is then G[1][2] ⊗ G[2][1] ⊗ G[2][5]). 
  corresponding_tuples_as_indices = vec([x for x = Base.Iterators.ProductIterator(Tuple(1:ngens(x) for x = G))])
  # In corresponding_tuples we store tuples of the actual generators, so in 
  # the example above we would store (G[1][2], G[2][1], G[2][5]).
  corresponding_tuples = map(index_tuple -> Tuple(map(index -> G[index][index_tuple[index]],1:length(index_tuple))), corresponding_tuples_as_indices)

  generating_tensors::Vector{elem_type(F)} = map(mF, map(tuple -> map(x -> parent(x) isa FreeMod ? x : repres(x), tuple), corresponding_tuples))
  s, emb = sub(F, generating_tensors)
  #s, emb = sub(F, vec([mF(x) for x = Base.Iterators.ProductIterator(Tuple(gens(x, ambient_free_module(x)) for x = G))]), :with_morphism)
  q::Vector{elem_type(F)} = vcat([vec([mF(x) for x = Base.Iterators.ProductIterator(Tuple(i == j ? rels(G[i]) : gens(ambient_free_module(G[i])) for i=1:length(G)))]) for j=1:length(G)]...) 
  local projection_map
  if length(q) != 0
    s, projection_map = quo(s, q)
  end

  tuples_pure_tensors_dict = IdDict(zip(corresponding_tuples_as_indices, gens(s)))
  set_attribute!(s, :show => Hecke.show_tensor_product, :tensor_product => G)

  
  function pure(tuple_elems::Union{SubquoModuleElem,FreeModElem}...)
    coeffs_tuples = vec([x for x in Base.Iterators.ProductIterator(coordinates.(tuple_elems))])
    res = zero(s)
    for coeffs_tuple in coeffs_tuples
      indices = map(first, coeffs_tuple)
      coeff_for_pure = prod(map(x -> x[2], coeffs_tuple))
      res += coeff_for_pure*tuples_pure_tensors_dict[indices]
    end
    return res
  end
  function pure(T::Tuple)
    return pure(T...)
  end

  decompose_generator = function(v::SubquoModuleElem)
    i = index_of_gen(v)
    return corresponding_tuples[i]
  end

  set_attribute!(s, :tensor_pure_function => pure, :tensor_generator_decompose_function => decompose_generator)

  if task == :none
    return s
  end

  return s, MapFromFunc(Hecke.TupleParent(Tuple([zero(g) for g = G])), s, pure, decompose_generator)
end


