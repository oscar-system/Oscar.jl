import Base: show, getindex

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

abstract type GraphicalModel{T <: AbstractGraph, L <: Union{NamedTuple, Nothing}} end

graph(M::GraphicalModel) = M.graph


@doc raw"""
    model_ring(GM::GraphicalModel)

Returns a polynomial ring where the indeterminants correspond to the defining parameters of the model.
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

# Below we make the ideal of the graph of a map. These methods allow us
# to work transparently with a polynomial and a rational map.
numerator(f::MPolyRingElem) = f
denominator(f::MPolyRingElem) = parent(f)(1)

#TODO update docs with an example?
@doc raw"""
    vanishing_ideal(M::GraphicalModel; algorithm::Symbol = :eliminate)

Compute the vanishing ideal for a graphical model `M`.
This is done by computing the kernel of the parametrization.
"""
function vanishing_ideal(GM::GraphicalModel; algorithm::Symbol = :eliminate)
  phi = parametrization(GM)
  if algorithm == :kernel
    invariants = kernel(phi)
  else
    S = domain(phi)
    R = codomain(phi)

    if R isa MPolyLocRing # Rational parametrization
      U = inverted_set(R)
      R = base_ring(R)
    else # Polynomial parametrization
      U = nothing
    end
    @req coefficient_ring(R) == coefficient_ring(S) "coefficient fields of parameter and model ring must be the same"
    elim_ring, elim_gens = polynomial_ring(coefficient_ring(R), vcat(symbols(R), symbols(S)))
    inject_R = hom(R, elim_ring, elim_gens[1:ngens(R)])
    inject_S = hom(S, elim_ring, elim_gens[ngens(R) + 1:ngens(elim_ring)])
    I = ideal([inject_S(s)*inject_R(denominator(phi(s))) - inject_R(numerator(phi(s))) for s in gens(S)])
    if U != nothing
      RU, iota = localization(elim_ring, inject_R(U))
      I = saturated_ideal(iota(I))
    end

    if algorithm == :eliminate
      invariants = eliminate(I, elim_gens[1:ngens(R)])
    elseif algorithm == :f4
      invariants = ideal(groebner_basis_f4(I, eliminate=ngens(R)))
      return invariants
    elseif algorithm == :markov
      o = lex(elim_ring)
      groebner_basis(I; ordering=o, algorithm=:markov)
      invariants = eliminate(I, o, ngens(R))
    else
      error("Didn't recognize  algorithm $algorithm")
    end
    project_S = hom(elim_ring, S, z -> z, [repeat([1], ngens(R)); gens(S)])
    invariants = project_S(invariants)
  end
  invariants
end
