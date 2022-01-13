export CoveredScheme

import Oscar.Graphs: Graph, Directed


mutable struct CoveredScheme{BRT, BRET, RT, RET} <:Scheme{BRT, BRET}
  patches::Vector{Spec{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}}
  glueings::Dict{Tuple{Spec{BRT, BRET, RT, RET}, Spec{BRT, BRET, RT, RET}}, Glueing{BRT, BRET, RT, RET}}

  # fields for caching
  glueing_graph::Graph{Directed}

  function CoveredScheme(
      patches::Vector{Spec{BRT, BRET, RT, RET, MPolyPowersOfElement{BRT, BRET, RT, RET}}},
      glueings::Dict{Tuple{Spec{BRT, BRET, RT, RET}}, Glueing{BRT, BRET, RT, RET}}
    ) where {BRT, BRET, RT, RET}
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
    return new{BRT, BRET, RT, RET}(X, glueings)
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

function are_glued(X::CoveredScheme, i::Int64, j::Int64)
  return exists_path(glueing_graph(X), i, j)
end

function getindex(X::CoveredScheme, a::Tuple{Int64, Int64})
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
    
