export Covering, patches, npatches, glueings, add_glueing!, standard_covering, glueing_graph, update_glueing_graph, transition_graph, edge_dict, disjoint_union
export fill_transitions!

import Oscar.Graphs: Graph, Directed, Undirected, add_edge!, vertices, edges, all_neighbors, neighbors, add_vertex!, nv, ne, has_edge


@Markdown.doc """
    Covering{SpecType<:Spec, GlueingType<:Glueing}

A covering of a scheme ``X`` by affine patches ``Uᵢ`` which are glued 
along isomorphisms ``gᵢⱼ : Uᵢ⊃ Vᵢⱼ →  Vⱼᵢ ⊂ Uⱼ``.

 * `SpecType` is the type of the affine patches;
 * `GlueingType` is the type of the glueings.

**Note:** The distinction between the different affine patches of the scheme 
is made from their hashes. Thus, an affine scheme must not appear more than once 
in any covering!
"""
mutable struct Covering{SpecType<:Spec, GlueingType<:Glueing}
  patches::Vector{SpecType}
  glueings::Dict{Tuple{SpecType, SpecType}, GlueingType}

  # fields for caching
  glueing_graph::Graph{Undirected}
  transition_graph::Graph{Undirected}
  edge_dict::Dict{Tuple{Int, Int}, Int}

  function Covering(
      patches::Vector{SpecType},
      glueings::Dict{Tuple{SpecType, SpecType}, GlueingType}
    ) where {SpecType<:Spec, GlueingType<:Glueing}
    n = length(patches)
    n > 0 || error("can not glue the empty scheme")
    kk = coefficient_ring(base_ring(OO(patches[1])))
    for i in 2:n
      kk == coefficient_ring(base_ring(OO(patches[i]))) || error("schemes are not defined over the same base ring")
    end
    # Check that no patch appears twice
    for i in 1:n-1
      for j in i+1:n
        patches[i] === patches[j] && error("affine schemes must not appear twice among the patches")
      end
    end
    for (X, Y) in keys(glueings)
      X in patches || error("glueings are not compatible with the patches")
      Y in patches || error("glueings are not compatible with the patches")
      if haskey(glueings, (Y, X))
	inverse(glueings[(X, Y)]) == glueings[(Y, X)] || error("glueings are not inverse of each other")
      else
	glueings[(Y, X)] = inverse(glueings[(X, Y)])
      end
    end
    return new{SpecType, GlueingType}(patches, glueings)
  end
end

patches(C::Covering) = C.patches
npatches(C::Covering) = length(C.patches)
glueings(C::Covering) = C.glueings
getindex(C::Covering, i::Int) = C.patches[i]
getindex(C::Covering, i::Int, j::Int) = glueings(C)[(patches(C)[i], patches(C)[j])]
getindex(C::Covering, X::SpecType, Y::SpecType) where {SpecType<:Spec} = glueings(C)[(X, Y)]
edge_dict(C::Covering) = C.edge_dict

# TODO: For some reason, the `indexin` method won't work. In the long 
# run, one should probably find out why and fix it. 
function getindex(C::Covering, X::SpecType) where {SpecType<:Spec} 
  for i in 1:length(patches(C))
    X === patches(C)[i] && return i
  end
  error("affine scheme could not be found among the patches")
end

function add_glueing!(C::Covering, G::Glueing)
  (X, Y) = patches(G)
  C.glueings[(X, Y)] = G
  C.glueings[(Y, X)] = inverse(G)
  return C
end

### The default constructor
# Returns a scheme in which every affine patch is only 
# glued to itself via the identity.
function Covering(patches::Vector{SpecType}) where {SpecType<:Spec}
  g = Dict{Tuple{SpecType, SpecType}, glueing_type(SpecType)}()
  for X in patches
    f = maximal_extension(X, X, fraction.(gens(OO(X))))
    g[X,X] = Glueing(X, X, f, f)
  end
  return Covering(patches, g)
end


### standard constructors 

function standard_covering(X::ProjectiveScheme{CRT}) where {CRT<:AbstractAlgebra.Field}
  CX = affine_cone(X)
  kk = base_ring(X)
  S = homogeneous_coordinate_ring(X)
  r = fiber_dimension(X)
  U = Vector{affine_patch_type(X)}()
  # TODO: Check that all weights are equal to one. Otherwise the routine is not implemented.
  s = symbols(S)
  for i in 0:r
    R, x = PolynomialRing(kk, [Symbol("("*String(s[k+1])*"//"*String(s[i+1])*")") for k in 0:r if k != i])
    push!(U, Spec(R))
  end
  result = Covering(U)
  for i in 2:r+1
    x = gens(base_ring(OO(U[1])))
    y = gens(base_ring(OO(U[i])))
    f = maximal_extension(U[1], U[i], vcat([1//x[i-1]], [x[k]//x[i-1] for k in 1:i-2], [x[k]//x[i-1] for k in i:r]))
    g = maximal_extension(U[i], U[1], vcat([y[k]//y[1] for k in 2:i-1], [1//y[1]], [y[k]//y[1] for k in i:r]))
    add_glueing!(result, Glueing(U[1], U[i], restriction(f, domain(f), domain(g)), restriction(g, domain(g), domain(f))))
  end
  fill_transitions!(result)
  return result
end

function standard_covering(X::ProjectiveScheme{CRT}) where {CRT<:MPolyQuoLocalizedRing}
  CX = affine_cone(X)
  Y = base_scheme(X)
  L = OO(Y)
  W = localized_ring(L)
  R = base_ring(L)
  kk = coefficient_ring(R)
  S = homogeneous_coordinate_ring(X)
  r = fiber_dimension(X)
  U = Vector{affine_patch_type(X)}()
  # TODO: Check that all weights are equal to one. Otherwise the routine is not implemented.
  s = symbols(S)
  for i in 0:r
    R, x = PolynomialRing(kk, [Symbol("("*String(s[k+1])*"//"*String(s[i+1])*")") for k in 0:r if k != i])
    F = Spec(R)
    patch, _, _ = product(F, Y)
    # TODO: Harvest the projection maps
    push!(U, patch)
  end
  result = Covering(U)
  for i in 2:r+1
    x = gens(base_ring(OO(U[1])))
    y = gens(base_ring(OO(U[i])))
    f = maximal_extension(U[1], U[i], vcat([1//x[i-1]], [x[k]//x[i-1] for k in 1:i-2], [x[k]//x[i-1] for k in i:r], x[r+1:end]))
    g = maximal_extension(U[i], U[1], vcat([y[k]//y[1] for k in 2:i-1], [1//y[1]], [y[k]//y[1] for k in i:r], y[r+1:end]))
    add_glueing!(result, Glueing(U[1], U[i], restriction(f, domain(f), domain(g)), restriction(g, domain(g), domain(f))))
  end
  fill_transitions!(result)
  return result
end

function glueing_graph(C::Covering) 
  if !isdefined(C, :glueing_graph)
    update_glueing_graph(C)
  end
  return C.glueing_graph
end

function update_glueing_graph(C::Covering)
  n = npatches(C)
  gg = Graph{Undirected}(n)
  for (X, Y) in keys(glueings(C))
    (U, V) = glueing_domains(C[X,Y])
    is_dense(U) && add_edge!(gg, C[X], C[Y])
    is_dense(V) && add_edge!(gg, C[Y], C[X])
  end
  C.glueing_graph = gg
  return gg
end

function transition_graph(C::Covering)
  if !isdefined(C, :transition_graph)
    p = length(keys(glueings(C)))
    edge_dict = Dict{Tuple{Int, Int}, Int}()
    edge_count = 1::Int
    C.transition_graph = Graph{Undirected}(0)
    for v in 1:nv(glueing_graph(C))
      W = neighbors(glueing_graph(C), v)
      for i in 1:length(W)-1
        for j in i+1:length(W)
          if is_dense(intersect(glueing_domains(C[W[i],v])[2], glueing_domains(C[v,W[j]])[1]))
            if !haskey(edge_dict, (W[i],v))
              edge_dict[(W[i],v)] = edge_dict[(v, W[i])] = edge_count
              edge_count+=1
              add_vertex!(C.transition_graph)
            end
            if !haskey(edge_dict, (v, W[j]))
              edge_dict[(v, W[j])] = edge_dict[(W[j], v)] = edge_count
              edge_count+=1
              add_vertex!(C.transition_graph)
            end
            add_edge!(C.transition_graph, edge_dict[(W[i],v)], edge_dict[(v,W[j])])
          end
        end
      end
    end
  end
  C.edge_dict = edge_dict
  return C.transition_graph
end

### fill transitions
# Whenever three schemes X ↩ U ↪ Y ↩ V ↪ Z are glued with 
# U ∩ V dense in both U and V, one can infer a glueing of 
# X and Z. This is done by crawling through the glueing graph, 
# updating the glueings, and proceeding until nothing more 
# can be done.
function fill_transitions!(C::Covering)
  gg = glueing_graph(C)
  dirty = true
  while dirty
    dirty = false
    for v in 1:nv(gg)
      W = neighbors(gg, v)
      for i in 1:length(W)-1
        for j in i+1:length(W)
          # TODO: replace the `is_dense` command by one that really checks that the 
          # intersection of U and V is dense in both U and V and not in their ambient 
          # variety. This implementation certainly works, but it is not as general 
          # as it could be, yet. 
          if !has_edge(gg, W[i], W[j]) && is_dense(intersect(glueing_domains(C[W[i],v])[2], glueing_domains(C[v,W[j]])[1]))
            new_glueing = maximal_extension(compose(C[W[i], v], C[v, W[j]]))
            add_glueing!(C, new_glueing)
            add_edge!(gg, W[i], W[j])
            dirty = true
          end
        end
      end
    end
  end
  return C
end

function disjoint_union(C1::T, C2::T) where {T<:Covering}
  C = Covering(vcat(patches(C1), patches(C2)))
  for (X, Y) in keys(glueings(C1))
    add_glueing!(C, C1[X, Y])
  end
  for (X, Y) in keys(glueings(C2))
    add_glueing!(C, C2[X, Y])
  end
  return C
end

mutable struct CoveringMorphism{SpecType<:Spec, CoveringType<:Covering, SpecMorType<:SpecMor}
  domain::CoveringType
  codomain::CoveringType
  morphisms::Dict{SpecType, SpecMorType} # on a patch X of the domain covering, this 
                                         # returns the morphism φ : X → Y to the corresponding 
                                         # patch Y of the codomain covering. 

  function CoveringMorphism(
      domain::CoveringType, 
      codomain::CoveringType, 
      morphisms::Dict{SpecType, SpecMorType}; 
      check::Bool=true
    ) where {
             SpecType<:Spec, 
             CoveringType<:Covering, 
             SpecMorType<:SpecMor
            }
    # TODO: check domain/codomain compatibility
    # TODO: if check is true, check that all morphisms glue.
    return new{SpecType, CoveringType, SpecMorType}(domain, codomain, morphisms)
  end
end


mutable struct CoveredScheme{CoveringType<:Covering, CoveringMorphismType<:CoveringMorphism}
  coverings::Vector{CoveringType}
  refinements::Dict{Tuple{CoveringType, CoveringType}, CoveringMorphismType}
  refinement_graph::Graph{Directed}

  function CoveredScheme(coverings::Vector{CoveringType}, refinements::Dict{Tuple{CoveringType, CoveringType}, CoveringMorphismType}) where {CoveringType<:Covering, CoveringMorphismType<:CoveringMorphism}
    # TODO: Check whether the refinements form a connected graph.
    return new{CoveringType, CoveringMorphismType}(coverings, refinements)
  end
end

