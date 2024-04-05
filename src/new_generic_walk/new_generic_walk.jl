include("new_generic_step.jl")
include("new_next_gamma.jl")

function new_generic_walk(G::Oscar.IdealGens, start::MonomialOrdering, target::MonomialOrdering)
    Lm = leading_term.(G, ordering = start)
    v = new_next_gamma(G, Lm, ZZ.([0]), start, target)
  
    while !isempty(v)
      G, Lm = new_generic_step(G, Lm, v, target)
  
      G = Oscar.IdealGens(G)
      v = new_next_gamma(G, Lm, v, start, target)
    end
    return G
  end
  

