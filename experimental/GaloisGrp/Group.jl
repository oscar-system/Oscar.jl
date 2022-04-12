"""
Tests if `U` is a maximal subgroup of `G`. Suboptimal algorithm...
  Missing in GAP
"""
function ismaximal(G::Oscar.PermGroup, U::Oscar.PermGroup)
  if !issubgroup(G, U)[1]
    return false
  end
  if order(G)//order(U) < 100
    t = right_transversal(G, U)[2:end] #drop the identity
    if any(x->order(sub(G, vcat(gens(U), [x]))[1]) < order(G), t)
      return false
    end
    return true
  end
  S = maximal_subgroup_reps(G) 
  error("not done yet")
  if any(x->x == U, S)
    return true
  end
  return false
end

"""
Tests if a conjugate of `V` by some element in `G` is a subgroup of `U`.
  Missing from GAP
"""
function isconjugate_subgroup(G::Oscar.PermGroup, U::Oscar.PermGroup, V::Oscar.PermGroup)
  if order(V) == 1
    return true, one(U)
  end
  local sigma
  while true
    sigma = rand(V)
    if order(sigma) > 1
      break
    end
  end
  s = short_right_transversal(G, U, sigma)
  for t = s
    if issubgroup(U, V^inv(t))[1]
      return true, inv(t)
    end
  end
  return false, one(U)
end

export maximal_subgroup_reps
function maximal_subgroup_reps(G::Oscar.GAPGroup)
  return Oscar._as_subgroups(G, GAP.Globals.MaximalSubgroupClassReps(G.X))
end

function low_index_subgroups(G::PermGroup, n::Int)
  ll = GAP.Globals.LowIndexSubgroups(G.X, n)
  return [Oscar._as_subgroup(G, x)[1] for x = ll]
end

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

# TODO: add a GSet Julia type which does something similar Magma's,
# or also to GAP's ExternalSet (but don't use ExternalSet, to avoid the overhead)

# TODO: add type BlockSystem

# TODO: add lots of more orbit related stuff

#provided by Thomas Breuer:

Hecke.orbit(G::PermGroup, i::Int) = Vector{Int}(GAP.Globals.Orbit(G.X, GAP.Obj(i)))
orbits(G::PermGroup) = Vector{Vector{Int}}(GAP.Globals.Orbits(G.X, GAP.GapObj(1:degree(G))))

function action_homomorphism(G::PermGroup, omega::AbstractVector{Int})
  mp = GAP.Globals.ActionHomomorphism(G.X, GAP.GapObj(omega))
  if mp == GAP.Globals.fail throw(ArgumentError("Invalid input")) end
  H = PermGroup(GAP.Globals.Range(mp))
  return GAPGroupHomomorphism(G, H, mp)
end

function block_system(G::PermGroup, B::Vector{Int})
  orb = GAP.Globals.Orbit(G.X, GAP.GapObj(B), GAP.Globals.OnSets)
  return Vector{Vector{Int}}(orb)
end

# given a perm group G and a block B, compute a homomorphism into Sym(B^G)
function action_on_blocks(G::PermGroup, B::Vector{Int})
  orb = GAP.Globals.Orbit(G.X, GAP.GapObj(B), GAP.Globals.OnSets)
  act = GAP.Globals.ActionHomomorphism(G.X, orb, GAP.Globals.OnSets)
  H = GAPWrap.Image(act)
  T = Oscar._get_type(H)
  H = T(H)
  return Oscar.GAPGroupHomomorphism(G, H, act)
end

function action_on_block_system(G::PermGroup, B::Vector{Vector{Int}})
  orb = GAP.GapObj(B, recursive = true)
  act = GAP.Globals.ActionHomomorphism(G.X, orb, GAP.Globals.OnSets)
  H = GAPWrap.Image(act)
  T = Oscar._get_type(H)
  H = T(H)
  return Oscar.GAPGroupHomomorphism(G, H, act)
end


@doc Markdown.doc"""
    short_right_transversal(G::PermGroup, H::PermGroup, s::PermGroupElem) ->

Determines representatives `g` for all right-cosets of `G` modulo `H`
  such that `H^g` contains the element `s`.
"""
function short_right_transversal(G::PermGroup, H::PermGroup, s::PermGroupElem)
  C = conjugacy_classes(H)
  cs = cycle_structure(s)
  can = PermGroupElem[]
  for c in C
    r = representative(c)
    if cs == cycle_structure(r)
      push!(can, r)
    end
  end

  R = PermGroupElem[]
  for c in can
    success, d = representative_action(G, c, s)
    if success
      push!(R, d)
      @assert c^R[end] == s
    end
  end

  S = PermGroupElem[]
  C = centralizer(G, s)[1]
  for r in R
    CH = centralizer(H^r, s)[1]
    for t = right_transversal(C, CH)
      push!(S, r*t)
    end
  end

  return S
end

