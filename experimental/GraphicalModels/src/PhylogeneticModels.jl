struct PhylogeneticModel
  graph::Graph{Directed}
  n_states::Int
  prob_ring::MPolyRing{QQFieldElem}
  fourier_ring::MPolyRing{QQFieldElem}
  trans_matrices::Dict{Edge, MatElem{QQMPolyRingElem}}
  fourier_params::Dict{Edge, Vector{QQMPolyRingElem}}
  group::Vector{Vector{Int64}}
end

# calling elements of the structure
@doc raw"""
    graph(pm::PhylogeneticModel)

Return the graph of a PhylogeneticModel `pm`.
"""
graph(pm::PhylogeneticModel) = pm.graph

@doc raw"""
    transition_matrices(pm::PhylogeneticModel)

Returns a dictionary between the edges of the tree specifying a PhylogeneticModel `pm` and attached transition matrices.
"""
transition_matrices(pm::PhylogeneticModel) = pm.trans_matrices

@doc raw"""
    number_states(pm::PhylogeneticModel)

Return the number of states of the PhylogeneticModel `pm`.
"""
number_states(pm::PhylogeneticModel) = pm.n_states

@doc raw"""
    probabilities_ring(pm::PhylogeneticModel)

Returns the ring of probability coordinates of the PhylogeneticModel `pm`.
"""
probabilities_ring(pm::PhylogeneticModel) = pm.prob_ring

@doc raw"""
    fourier_ring(pm::PhylogeneticModel)

    Returns the ring of Fourier coordinates of the PhylogeneticModel `pm`.
"""
fourier_ring(pm::PhylogeneticModel) = pm.fourier_ring

@doc raw"""
    fourier_param(pm::PhylogeneticModel)

Return the Fourier parameters of the PhylogeneticModel `pm` as a vector of eigenvalues of the transition matrices.
"""
fourier_parameters(pm::PhylogeneticModel) = pm.fourier_params

group_model(pm::PhylogeneticModel) = pm.group


#define group-based models
@doc raw"""
    jukes_cantor_model(graph::Graph{Directed})

Creates a `PhylogeneticModel` based on `graph` whose transition matrices are Jukes Cantor matrices. 

# Examples
```jldoctest
julia> pm = jukes_cantor_model(graph_from_edges(Directed,[[4,1],[4,2],[4,3]]));
```
We can recover the information of `pm`
```jldoctest
julia> graph(pm)
Directed graph with 4 nodes and the following edges:
(4, 1)(4, 2)(4, 3)

julia> number_states(pm)
4

julia> prob_ring(pm)
fourier_ring(pm)
trans_matrices(pm)
fourier_params(pm)
group(pm)
```
"""
function jukes_cantor_model(graph::Graph{Directed})
    ns = 4
    ne = n_edges(graph)
    R, list_a, list_b = polynomial_ring(QQ, :a => 1:ne, :b => 1:ne)
    
    matrices = Dict{Edge, MatElem}(e => matrix(R, [
      a b b b
      b a b b
      b b a b
      b b b a]) for (a,b,e) in zip(list_a, list_b, edges(graph))
    )

    S, list_x = polynomial_ring(QQ, :x => (1:ne, 1:2))
    fourier_param = Dict{Edge, Vector{QQMPolyRingElem}}(e => 
            [list_x[i,1], list_x[i,2], list_x[i,2], list_x[i,2]] for (i, e) in zip(1:ne, edges(graph)))
    
    group = [[0,0], [0,1], [1,0], [1,1]]

    return PhylogeneticModel(ns, graph, R , S, matrices, fourier_param, group)
end

# Kimura 2
function kimura2_model(graph::Graph{Directed})
    ns = 4
    ne = n_edges(graph)
    R, list_a, list_b, list_c = polynomial_ring(QQ, :a => 1:ne, :b => 1:ne, :c => 1:ne)
    
    matrices = Dict{Edge, MatElem}(e => matrix(R, [
      a b c b
      b a b c
      c b a b
      b c b a]) for (a,b,c,e) in zip(list_a, list_b, list_c, edges(graph))
    )

    S, list_x = polynomial_ring(QQ, :x => (1:ne, 1:3))
    fourier_param = Dict{Edge, Vector{QQMPolyRingElem}}(e => 
            [list_x[i,1], list_x[i,2], list_x[i,3], list_x[i,2]] for (i, e) in zip(1:ne, edges(graph)))
    
    group = [[0,0], [0,1], [1,0], [1,1]]

    return PhylogeneticModel(ns, graph, R, S, matrices, fourier_param, group)
end

# Kimura 3
function kimura3_model(graph::Graph{Directed})
    ns = 4
    ne = n_edges(graph)
    R, list_a, list_b , list_c, list_d= polynomial_ring(QQ, :a => 1:ne, :b => 1:ne, :c => 1:ne, :d => 1:ne)
    
    matrices = Dict{Edge, MatElem}(e => matrix(R, [
      a b c d
      b a d c
      c d a b
      d c b a]) for (a,b,c,d,e) in zip(list_a, list_b, list_c, list_d, edges(graph))
    )

    S, list_x = polynomial_ring(QQ, :x => (1:ne, 1:4))
    fourier_param = Dict{Edge, Vector{QQMPolyRingElem}}(e => 
            [list_x[i,1], list_x[i,2], list_x[i,3], list_x[i,4]] for (i, e) in zip(1:ne, edges(graph)))
    
    group = [[0,0], [0,1], [1,0], [1,1]]

    return PhylogeneticModel(ns, graph, R, S, matrices, fourier_param, group)
end

# general Markov model
# work in progress
# not sure whether to use this https://docs.sciml.ai/Symbolics/stable/manual/arrays/
# or array of polynomial ring elements in OSCAR
function generalmarkov_model(graph::Graph{Directed})
    ns = 4
    ne = n_edges(graph)
    R, list_m = polynomial_ring(QQ, :m => 1:(ne^2))
    
    matrices = Dict{Edge, MatElem}(e => matrix(R, [
      m11 m12 m13 m14
      m11 m12 m13 m14
      b b a b
      b b b a]) for (a,b,e) in zip(list_a, list_b, edges(graph))
    )
    return PhylogeneticModel(ns, graph, R, matrices)
end

# more models to come, Ch Mar 21


#### SPECIALIZED FOURIER TRANSFORM MATRIX ####

function hadamardmatrix()
  H = [1 1 1 1
       1 1 -1 -1 
       1 -1 1 -1
       1 -1 -1 1]
  return H
end

# TODO: line 138 mutates f_equivclasses, i.e. an object from outside the function.
# This makes the function a !-function. We do not want f_equivclasses to be changed, 
# so I add the missing equivalence class again in line 175. 
# We should find out whether there is a better way to do this.
function specialized_fourier_transform(pm::PhylogeneticModel, p_equivclasses::Dict{Vector{Vector{Int64}}, QQMPolyRingElem},f_equivclasses::Dict{Vector{Vector{Int64}}, QQMPolyRingElem})
  R = probabilities_ring(pm)

  class0 = findall(x -> x ==0, f_equivclasses)[1]
  delete!(f_equivclasses, class0)

  np = length(p_equivclasses)
  nq = length(f_equivclasses)

  ## We need to sort the equivalence classes: both inside each class as well as the collection of classes. 
  p_equivclasses_sorted = collect(keys(p_equivclasses))
  for p_eqclass in p_equivclasses_sorted
      sort!(p_eqclass)
  end
  sort!(p_equivclasses_sorted)

  f_equivclasses_sorted = collect(keys(f_equivclasses))
  for f_eqclass in f_equivclasses_sorted
      sort!(f_eqclass)
  end
  sort!(f_equivclasses_sorted)

  H = hadamardmatrix()

  specialized_ft_matrix = R.(Int.(zeros(nq, np)))
  for i in 1:nq
      current_fourier_classes = f_equivclasses_sorted[i]
      for j in 1:np
          current_prob_classes = p_equivclasses_sorted[j]
          current_entriesin_M = [prod([H[y,x] for (x,y) in zip(p,q)]) for p in current_prob_classes, q in current_fourier_classes]
          specialized_ft_matrix[i,j] = R.(1//(length(current_prob_classes)*length(current_fourier_classes))*sum(current_entriesin_M))
      end
  end
  get!(f_equivclasses, class0, R(0))
  return specialized_ft_matrix
end

function inverse_specialized_fourier_transform(pm::PhylogeneticModel, p_equivclasses::Dict{Vector{Vector{Int64}}, QQMPolyRingElem},f_equivclasses::Dict{Vector{Vector{Int64}}, QQMPolyRingElem})
  R = probabilities_ring(pm)

  class0 = findall(x -> x ==0, f_equivclasses)[1]
  delete!(f_equivclasses, class0)

  np = length(p_equivclasses)
  nq = length(f_equivclasses)

  ## We need to sort the equivalence classes: both inside each class as well as the collection of classes. 
  p_equivclasses_sorted = collect(keys(p_equivclasses))
  for p_eqclass in p_equivclasses_sorted
      sort!(p_eqclass)
  end
  sort!(p_equivclasses_sorted)

  f_equivclasses_sorted = collect(keys(f_equivclasses))
  for f_eqclass in f_equivclasses_sorted
      sort!(f_eqclass)
  end
  sort!(f_equivclasses_sorted)

  H = hadamardmatrix()
  Hinv = 1//4 * H 

  inverse_spec_ft_matrix = R.(Int.(zeros(np, nq)))
  for i in 1:np
      current_prob_class = p_equivclasses_sorted[i]
      for j in 1:nq
          current_fourier_class = f_equivclasses_sorted[j]
          current_entriesin_Minv = [prod([Hinv[x,y] for (x,y) in zip(p,q)]) for p in current_prob_class, q in current_fourier_class] 
          inverse_spec_ft_matrix[i,j] = R.(sum(current_entriesin_Minv))
      end
  end
  get!(f_equivclasses, class0, R(0))
  return inverse_spec_ft_matrix
end

-
# @doc raw"""
#     my_access_func(S::ExampleStruct)

# This is a dummy sample function to teach the use of docstrings.
# """
