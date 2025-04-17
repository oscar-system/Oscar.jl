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
function tensor_product(G::ModuleFP...; task::Symbol = :none)
  resols = AbsHyperComplex[]
  augs = ModuleFPHom[]
  for M in G
    res, aug = free_resolution(SimpleFreeResolution, M)
    push!(resols, res)
    push!(augs, aug[0])
  end
  res_prod = tensor_product(resols)
  tot = total_complex(res_prod)
  pres = map(tot, 1)
  I, inc_I = image(pres)
  result, pr_res = quo(tot[0], I)
  
  # assemble the multiplication and decomposition functions
  z = Tuple([0 for _ in 1:length(G)])
  @assert _is_tensor_product(res_prod[z])[1]
  function pure(tuple_elems::Union{SubquoModuleElem,FreeModElem}...)
    w = [preimage(augs[i], x) for (i, x) in enumerate(tuple_elems)]
    free_pure = tensor_pure_function(res_prod[z])
    ww = free_pure(w...)
    return pr_res(canonical_injection(tot[0], 1)(ww))
  end
  function pure(T::Tuple)
    return pure(T...)
  end
  
  decompose_generator = function(v::SubquoModuleElem)
    w = canonical_projection(tot[0], 1)(preimage(pr_res, v))
    w_dec = tensor_generator_decompose_function(res_prod[z])(w)
    return Tuple([augs[i](x) for (i, x) in enumerate(w_dec)])
  end

  set_attribute!(result, :tensor_pure_function => pure, :tensor_generator_decompose_function => decompose_generator)
  set_attribute!(result, :show => Hecke.show_tensor_product, :tensor_product => G)
  @assert _is_tensor_product(result)[1]
  
  if task == :none
    return result
  end

  return result, MapFromFunc(Hecke.TupleParent(Tuple([zero(g) for g = G])), result, pure, decompose_generator)
end


