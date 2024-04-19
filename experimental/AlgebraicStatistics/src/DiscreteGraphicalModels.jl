# creates matrix of parameters which correspond to the directed edges of a digraph
# creates the parameterization of the model as a ring map
function markov_ring(G::Graph{Undirected})

end


# creates a ring with variables corresponding to the joint probabilities of discrete random variables
# associated to the vertices of G which have states (d1, ... , dn)
function markov_ring(G::Graph{Directed})

end



# creates the parameterization of the model as a ring map
function discrete_parameterization(R::MarkovRing)

end



# computes the vanishing ideal of a discrete graphical model
function discrete_vanishing_ideal(R::MarkovRing)

end