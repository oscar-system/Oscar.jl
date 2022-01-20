export CoveredScheme

import Oscar.Graphs: Graph, Directed


mutable struct CoveredScheme{SpecType<:Spec, GlueingType<:Glueing, BRT, BRET} <:Scheme{BRT, BRET}
  patches::Vector{SpecType}
  glueings::Dict{Tuple{SpecType, SpecType}, GlueingType}

  # fields for caching
  glueing_graph::Graph{Directed}

  function CoveredScheme(
      patches::Vector{SpecType},
      glueings::Dict{Tuple{SpecType}, GlueingType}
    ) where {SpecType<:Spec, GlueingType<:Glueing}
    n = length(patches)
    n > 0 || error("can not glue the empty scheme")
    kk = coefficient_ring(base_ring(OO(patches[1])))
    for i in 2:n
      kk == coefficient_ring(base_ring(OO(patches[i]))) || error("schemes are not defined over the same base ring")
    end
    for (X, Y) in keys(glueings)
      X in patches || error("glueings are not compatible with the patches")
      Y in patches || error("glueings are not compatible with the patches")
      if haskey(g, (Y, X))
	inverse(g[(X, Y)]) == g[(Y, X)] || error("glueings are not inverse of each other")
      else
	g[(Y, X)] = inverse(g[(X, Y)])
      end
    end
    X = patches[1]
    return new{SpecType, GlueingType, typeof(base_ring(X)), elem_type(base_ring(X)) }(patches, glueings)
  end
end

patches(X::CoveredScheme) = X.patches
npatches(X::CoveredScheme) = length(X.patches)
glueings(X::CoveredScheme) = X.glueings
getindex(X::CoveredScheme, i::Int) = X.patches[i]

function glueing_graph(X::CoveredScheme) 
  if !isdefined(X, :glueing_graph)
    n = npatches(X)
    gg = Graph{Directed}(n)
    for (i,j) in keys(glueings(X))
      g = X[(i,j)] 
      Y = X[i]
      Z = X[j]
      (U, V) = glueing_domains(g)
      is_dense(U, Y) && add_edge!(gg, i, j)
      is_dense(V, Z) && add_edge!(gg, j, i)
    end
    X.glueing_graph = gg
  end
  return X.glueing_graph
end

function are_glued(X::CoveredScheme, U::T, V::T) where {T<:Spec}
  return exists_path(glueing_graph(X), U, V)
end

function getindex(X::CoveredScheme, a::Tuple{T, T}) where {T<:Spec}
  if haskey(glueings(X), a)
    return glueings(X)[a]
  end
  if are_glued(X, a[1], a[2])
    p = get_path(glueing_graph(X), a[1], a[2])
    # construct the glueing along p
    glueings(X)[a] = g
    return g
  end
  glueings(X)[a] = trivial_glueing(X[a[1]], X[a[2]])
  glueings(X)[(a[2], a[1])] = inverse(glueings(X)[a])
  return glueings(X)[a]
end
    
