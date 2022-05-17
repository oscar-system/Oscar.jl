
"""
An ascending chain of minimal supergroups linking `U` and `G`.
Not a good algorithm, but Max' version `RefinedChain` does not produce
minimal supergroups, ie. it fails in ``C_4``.
"""
function maximal_subgroup_chain(G::PermGroup, U::PermGroup)
  l = [G]
  while order(l[end]) > order(U)
    m = maximal_subgroups(l[end])
    push!(l, m[findfirst(x -> issubgroup(x, U)[1], m)])
  end
  return reverse(l)

  l = GAP.Globals.AscendingChain(G.X,U.X)
  map(GAP.Globals.MaximalSubgroupClassReps, l)
  ll = GAP.Globals.RefinedChain(G.X,l)
  return [Oscar._as_subgroup(G, x)[1] for x = ll]
end

function block_system(G::PermGroup, B::Vector{Int})
  return collect(orbit(G, on_sets, B))
end

# given a perm group G and a block B, compute a homomorphism into Sym(B^G)
function action_on_blocks(G::PermGroup, B::Vector{Int})
  Omega = gset(G, on_sets, [B])
  return action_homomorphism(Omega)
end

function action_on_block_system(G::PermGroup, B::Vector{Vector{Int}})
  Omega = gset(G, on_sets, B)
  set_attribute!(Omega, :elements => Omega.seeds)
  return action_homomorphism(Omega)
end
