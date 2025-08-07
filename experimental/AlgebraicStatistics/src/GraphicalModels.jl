import Base.show

###################################################################################
#
#       Additional functions for matrices
#
###################################################################################

function adjugate(A::Generic.MatSpaceElem)
  if nrows(A) != ncols(A)
    error("adjugate can only be computed for a square matrix")
  end
  adj, _ = pseudo_inv(A)
  # pseudo_inv guarantees only that A*adj == d*Id but that allows adj and d
  # to both differ from the adjugate and the determinant, respectively, by
  # a sign. We compute one cofactor by hand to make sure we get sign right.
  n = nrows(A)
  c = det(A[1:n-1, 1:n-1])
  if c == adj[n,n]
    return adj
  elseif c == -adj[n,n]
    return -adj
  else # unlikely, but fall back to computing many determinants
    return matrix(base_ring(A), [[(-1)^(i+j)*det(A[setdiff(1:n, j), setdiff(1:n, i)]) for j in 1:n] for i in 1:n])
  end
end

###################################################################################
#
#       Generic Graphical Models
#
###################################################################################
const GraphTypes = Union{Directed, Undirected, Mixed}

abstract type GraphicalModel{T <: GraphTypes, L <: Union{NamedTuple, Nothing}} end

graph(M::GraphicalModel) = M.graph

@doc raw"""
    model_ring(GM::GraphicalModel)

Returns a polynomial ring where the indeterminants correspond to the defining parameters of the model.
```
"""
model_ring(M::T) where T <: GraphicalModel = error("Please implement the method model_ring for $T")

@doc raw"""
    parameter_ring(GM::GraphicalModel)

Returns a polynomial ring where the indeterminants parametrize the variety of possible defining parameter values for the model.
"""
parameter_ring(M::T) where T <: GraphicalModel = error("Please implement the method parameter_ring for $T")

@doc raw"""
    parametrization(GM::GraphicalModel)

Returns a map from the model ring to the parameter ring.
"""
parametrization(M::T) where T <: GraphicalModel = error("Please implement the method parametrization for $T")

#TODO update docs with an example?
@doc raw"""
    vanishing_ideal(M::GraphicalModel)

Compute the vanishing ideal for a graphical model `M`.
This is done by computing the kernel of the parametrization.
"""
function vanishing_ideal(M::GraphicalModel)
  kernel(parametrization(M))
end
