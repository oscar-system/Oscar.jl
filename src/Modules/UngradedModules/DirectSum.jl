##################################################
# direct product
##################################################
@doc raw"""
    direct_product(F::FreeMod{T}...; task::Symbol = :prod) where T

Given free modules $F_1\dots F_n$, say, return the direct product $\prod_{i=1}^n F_i$.

Additionally, return
- a vector containing the canonical projections  $\prod_{i=1}^n F_i\to F_i$ if `task = :prod` (default),
- a vector containing the canonical injections  $F_i\to\prod_{i=1}^n F_i$ if `task = :sum`,
- two vectors containing the canonical projections and injections, respectively, if `task = :both`,
- none of the above maps if `task = :none`.
"""
function direct_product(F::FreeMod{T}...; task::Symbol = :prod) where {T}
  R = base_ring(F[1])
  G = FreeMod(R, sum([rank(f) for f = F]))
  all_graded = all(is_graded, F)
  if all_graded
    G.d = vcat([f.d for f in F]...)
  end
  G.S = []
  for i = 1:length(F)
    s = "("
    for j=1:i-1
      s *= "0, "
    end
    e = ""
    if i<length(F)
      e*=", "
    end
    for j=i+1:length(F)-1
      e *= "0, "
    end
    if i<length(F)
      e *= "0"
    end
    e*=")"

    for t = F[i].S
      push!(G.S, Symbol(s*string(t)*e))
    end
  end
  set_attribute!(G, :show => Hecke.show_direct_product, :direct_product => F)
  emb = []
  pro = []
  projection_dictionary = IdDict{Int,ModuleFPHom}()
  injection_dictionary = IdDict{Int,ModuleFPHom}()
  set_attribute!(G, :projection_morphisms => projection_dictionary, :injection_morphisms => injection_dictionary)
  i = 0
  for f = F
    if task in [:sum, :both]
      push!(emb, hom(f, G, Vector{elem_type(G)}([gen(G, j+i) for j=1:ngens(f)]); check=false))
      injection_dictionary[length(emb)] = emb[length(emb)]
    end
    if task in [:prod, :both]
      push!(pro, hom(G, f, vcat(elem_type(f)[zero(f) for j=1:i], gens(f), elem_type(f)[zero(f) for j=i+ngens(f)+1:ngens(G)]); check=false))
      projection_dictionary[length(pro)] = pro[length(pro)]
    end
    i += ngens(f)
  end
  if task == :none
    return G
  elseif task == :sum
    return G, emb
  elseif task == :prod
    return G, pro
  elseif task == :both
    return G, pro, emb
  end
end

@doc raw"""
    direct_product(M::ModuleFP{T}...; task::Symbol = :prod) where T

Given modules $M_1\dots M_n$, say, return the direct product $\prod_{i=1}^n M_i$.

Additionally, return
- a vector containing the canonical projections  $\prod_{i=1}^n M_i\to M_i$ if `task = :prod` (default),
- a vector containing the canonical injections  $M_i\to\prod_{i=1}^n M_i$ if `task = :sum`,
- two vectors containing the canonical projections and injections, respectively, if `task = :both`,
- none of the above maps if `task = :none`.
"""
function direct_product(M::ModuleFP{T}...; task::Symbol = :prod) where T
  F, pro, mF = direct_product([ambient_free_module(x) for x = M]..., task = :both)
  s, emb_sF = sub(F, vcat([elem_type(F)[mF[i](y) for y = ambient_representatives_generators(M[i])] for i=1:length(M)]...))
  q::Vector{elem_type(F)} = vcat([elem_type(F)[mF[i](y) for y = rels(M[i])] for i=1:length(M)]...)
  pro_quo = nothing
  if length(q) != 0
    s, pro_quo = quo(s, q)
  end
  set_attribute!(s, :show => Hecke.show_direct_product, :direct_product => M)
  projection_dictionary = IdDict{Int,ModuleFPHom}()
  injection_dictionary = IdDict{Int,ModuleFPHom}()
  set_attribute!(s, :projection_morphisms => projection_dictionary, :injection_morphisms => injection_dictionary)
  if task == :none
    return s
  end
  if task == :prod || task != :sum
    if pro_quo === nothing
      for i=1:length(pro)
        pro[i] = hom(s, M[i], Vector{elem_type(M[i])}([M[i](pro[i](emb_sF(gen))) for gen in gens(s)]); check=false) # TODO distinction between pro on the left and pro on the right side!
        projection_dictionary[i] = pro[i]
      end
    else
      for i=1:length(pro)
        pro[i] = hom(s, M[i], Vector{elem_type(M[i])}([M[i](pro[i](emb_sF(preimage(pro_quo,gen)))) for gen in gens(s)]); check=false)
        projection_dictionary[i] = pro[i]
      end
    end
    if task == :prod
      return s, pro
    end
  end
  if task == :sum || task != :prod
    if pro_quo === nothing
      for i=1:length(mF)
        mF[i] = hom(M[i], s, Vector{elem_type(s)}([preimage(emb_sF, mF[i](repres(g))) for g in gens(M[i])]); check=false)
        injection_dictionary[i] = mF[i]
      end
    else
      for i=1:length(mF)
        mF[i] = hom(M[i], s, Vector{elem_type(s)}([pro_quo(preimage(emb_sF, mF[i](repres(g)))) for g in gens(M[i])]); check=false)
        injection_dictionary[i] = mF[i]
      end
    end
    if task == :sum
      return s, mF
    else
      return s, pro, mF
    end
  end
end
##################################################
# direct sum
##################################################
@doc raw"""
    direct_sum(M::ModuleFP{T}...; task::Symbol = :sum) where T

Given modules $M_1\dots M_n$, say, return the direct sum $\bigoplus_{i=1}^n M_i$.  
 
Additionally, return 
- a vector containing the canonical injections  $M_i\to\bigoplus_{i=1}^n M_i$ if `task = :sum` (default),
- a vector containing the canonical projections  $\bigoplus_{i=1}^n M_i\to M_i$ if `task = :prod`,
- two vectors containing the canonical injections and projections, respectively, if `task = :both`,
- none of the above maps if `task = :none`.
"""
function direct_sum(M::ModuleFP{T}...; task::Symbol = :sum) where {T}
  res = direct_product(M...; task)
  if task == :sum || task == :prod
    ds, f = res
    set_attribute!(ds, :show => Hecke.show_direct_sum, :direct_sum => M)
    return ds, f
  elseif task == :both
    ds, p, i = res
    set_attribute!(ds, :show => Hecke.show_direct_sum, :direct_sum => M)
    return ds, i, p
  else
    set_attribute!(res, :show => Hecke.show_direct_sum, :direct_sum => M)
    return res
  end
end

⊕(M::ModuleFP...) = direct_sum(M..., task = :none)

@doc raw"""
    canonical_injections(G::ModuleFP)

Return the canonical injections from all components into $G$
where $G = G_1 \oplus \cdot \oplus G_n$.
"""
function canonical_injections(G::ModuleFP)
  H = get_attribute(G, :direct_product)
  @req H !== nothing "module not a direct product"
  return [canonical_injection(G, i) for i in 1:length(H)]
end

@doc raw"""
    canonical_injection(G::ModuleFP, i::Int)

Return the canonical injection $G_i \to G$ where $G = G_1 \oplus \cdot \oplus G_n$.
"""
function canonical_injection(G::ModuleFP, i::Int)
  H = get_attribute(G, :direct_product)
  @req H !== nothing "module not a direct product"
  injection_dictionary = get_attribute(G, :injection_morphisms)
  if haskey(injection_dictionary, i)
    return injection_dictionary[i]
  end
  @req 0 < i <= length(H) "index out of bound"
  j = sum(ngens(H[l]) for l in 1:i-1; init=0)
  emb = hom(H[i], G, Vector{elem_type(G)}([G[l+j] for l in 1:ngens(H[i])]); check=false)
  injection_dictionary[i] = emb
  return emb
end

@doc raw"""
    canonical_projections(G::ModuleFP)

Return the canonical projections from $G$ to all components
where $G = G_1 \oplus \cdot \oplus G_n$.
"""
function canonical_projections(G::ModuleFP)
  H = get_attribute(G, :direct_product)
  @req H !== nothing "module not a direct product"
  return [canonical_projection(G, i) for i in 1:length(H)]
end

@doc raw"""
    canonical_projection(G::ModuleFP, i::Int)

Return the canonical projection $G \to G_i$ where $G = G_1 \oplus \cdot \oplus G_n$.
"""
function canonical_projection(G::ModuleFP, i::Int)
  H = get_attribute(G, :direct_product)
  @req H !== nothing "module not a direct product"
  projection_dictionary = get_attribute(G, :projection_morphisms)
  if haskey(projection_dictionary, i)
    return projection_dictionary[i]
  end
  @req 0 < i <= length(H) "index out of bound"
  j = sum(ngens(H[l]) for l in 1:i-1; init=0) 
  pro = hom(G, H[i], Vector{elem_type(H[i])}(vcat([zero(H[i]) for l in 1:j], gens(H[i]), [zero(H[i]) for l in 1+j+ngens(H[i]):ngens(G)])); check=false)
  projection_dictionary[i] = pro
  return pro
end
    

