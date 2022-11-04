export neighbor_patches, update_glueing_graph, transition_graph, fill_transitions!, add_glueing!, all_patches
########################################################################
# Finding patches                                                      #
########################################################################

function getindex(C::Covering, X::AbsSpec)
  for i in 1:length(patches(C))
    X === patches(C)[i] && return i
  end
  error("affine scheme could not be found among the patches")
end

function Base.in(U::AbsSpec, C::Covering)
  for i in 1:npatches(C)
    U === C[i] && return true
  end
  return false
  ### Affine refinements not implemented at the moment!
#  for i in 1:npatches(C)
#    if haskey(affine_refinements(C), C[i])
#      V = affine_refinements(C)[C[i]]
#      for (V, a) in affine_refinements(C)[C[i]]
#        U in affine_patches(V) && return true
#      end
#    end
#  end
#  return false
end

function Base.indexin(U::AbsSpec, C::Covering)
  for i in 1:npatches(C)
    U === C[i] && return (i, 0, 0)
  end
  return (0,0,0)
  ### Affine refinements not implemented at the moment!
#  for i in 1:npatches(C)
#    if haskey(affine_refinements(C), C[i])
#      V = affine_refinements(C)[C[i]]
#      for j in 1:length(V)
#        (W, a) = V[j]
#        for k in 1:length(gens(W))
#          U === W[k] && return (i, j, k)
#        end
#      end
#    end
#  end
#  return (0,0,0)
end

########################################################################
# Functionality for lazy glueing                                       #
########################################################################
function neighbor_patches(C::Covering, U::AbsSpec)
  gg = glueing_graph(C)
  n = neighbors(gg, C[U])
  return [C[i] for i in n]
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
    C.edge_dict = edge_dict
  end
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

########################################################################
# Iteration over patches                                               #
########################################################################
function Base.length(C::Covering)
  result = 0
  for U in patches(C)
    result = result + 1
    if haskey(affine_refinements(C), U)
      for W in affine_refinements(C)[U]
        result = result + length(gens(W))
      end
    end
  end
  return result
end

function all_patches(C::Covering)
  result = Vector{AbsSpec}()
  for U in patches(C)
    push!(result, U)
    if haskey(affine_refinements(C), U)
      for (W, a) in affine_refinements(C)[U]
        result = vcat(result, affine_patches(W))
      end
    end
  end
  return result
end

function Base.iterate(C::Covering, s::Int=1)
  U = all_patches(C)
  s > length(U) && return nothing
  return U[s], s+1
end

Base.eltype(C::Covering) = AbsSpec

########################################################################
# Building a Covering                                                  #
########################################################################
function add_glueing!(C::Covering, G::AbsGlueing)
  (X, Y) = patches(G)
  C.glueings[(X, Y)] = G
  C.glueings[(Y, X)] = inverse(G)
  return C
end

########################################################################
# Printing                                                             #
########################################################################
function Base.show(io::IO, C::Covering) 
  print(io, 
          "Covering with $(npatches(C)) patch" * 
          (npatches(C) == 1 ? "" : "es"))
#           " and glueing graph")
#   print(io, glueing_graph(C))
end

