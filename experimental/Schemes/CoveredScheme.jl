export Covering, patches, npatches, glueings, add_glueing!, standard_covering, glueing_graph, update_glueing_graph, transition_graph, edge_dict, disjoint_union, neighbor_patches, affine_refinements, add_affine_refinement!, all_patches

export fill_transitions!

export affine_patch_type, glueing_type

export CoveringMorphism
export morphism_type, morphisms

export CoveredScheme
export empty_covered_scheme
export coverings, refinements, default_covering, set_name!, name_of, has_name, dim
export covering_type, covering_morphism_type, affine_patch_type, covered_scheme_type

import Oscar.Graphs: Graph, Directed, Undirected, add_edge!, edges, all_neighbors, neighbors, add_vertex!, nv, ne, has_edge

export CoveredSchemeMorphism, domain, codomain, covering_morphism

export simplify

@Markdown.doc """
    Covering{SpecType<:Spec, GlueingType<:Glueing, SpecOpenType<:SpecOpen}

A covering of a scheme ``X`` by affine patches ``Uᵢ`` which are glued 
along isomorphisms ``gᵢⱼ : Uᵢ⊃ Vᵢⱼ →  Vⱼᵢ ⊂ Uⱼ``.

 * `SpecType` is the type of the affine patches;
 * `GlueingType` is the type of the glueings;
 * `SpecOpenType` is the type of the affine refinements of the ``Uᵢ``.

**Note:** The distinction between the different affine patches of the scheme 
is made from their hashes. Thus, an affine scheme must not appear more than once 
in any covering!
"""
mutable struct Covering{SpecType<:Spec, GlueingType<:Glueing, SpecOpenType<:SpecOpen, RingElemType<:RingElem}
  patches::Vector{SpecType} # the basic affine patches of X
  glueings::Dict{Tuple{SpecType, SpecType}, GlueingType} # the glueings of the basic affine patches
  affine_refinements::Dict{SpecType, Vector{Tuple{SpecOpenType, Vector{RingElemType}}}} # optional lists of refinements 
      # of the basic affine patches.
      # These are stored as pairs (U, a) where U is a 'trivial' SpecOpen, 
      # meaning that its list of hypersurface equation (f₁,…,fᵣ) has empty 
      # intersection in the basic affine patch X and hence satisfies 
      # some equality 1 ≡ a₁⋅f₁ + a₂⋅f₂ + … + aᵣ⋅fᵣ on X. 
      # Since the coefficients aᵢ of this equality are crucial for computations, 
      # we store them in an extra tuple. 

  # fields for caching
  glueing_graph::Graph{Undirected}
  transition_graph::Graph{Undirected}
  edge_dict::Dict{Tuple{Int, Int}, Int}

  function Covering(
      patches::Vector{SpecType},
      glueings::Dict{Tuple{SpecType, SpecType}, GlueingType};
      affine_refinements::Dict{
          SpecType, 
          Vector{Tuple{SpecOpenType, Vector{RingElemType}}}
         }=Dict{SpecType, Vector{Tuple{open_subset_type(SpecType), Vector{elem_type(ring_type(SpecType))}}}}(),
      check::Bool=true
    ) where {SpecType<:Spec, GlueingType<:Glueing, SpecOpenType<:SpecOpen, RingElemType<:RingElem}
    n = length(patches)
    n > 0 || error("can not glue the empty scheme")
    kk = coefficient_ring(base_ring(OO(patches[1])))
    for i in 2:n
      kk == coefficient_ring(base_ring(OO(patches[i]))) || error("schemes are not defined over the same base ring")
    end
    # Check that no patch appears twice
    for i in 1:n-1
      for j in i+1:n
        patches[i] === patches[j] && error("affine schemes must not appear twice among the patches")
      end
    end
    for (X, Y) in keys(glueings)
      X in patches || error("glueings are not compatible with the patches")
      Y in patches || error("glueings are not compatible with the patches")
      if haskey(glueings, (Y, X))
        if check
          inverse(glueings[(X, Y)]) == glueings[(Y, X)] || error("glueings are not inverse of each other")
        end
      else
	glueings[(Y, X)] = inverse(glueings[(X, Y)])
      end
    end

    # check the affine refinements
    for U in keys(affine_refinements)
      for (V, a) in affine_refinements[U]
        ambient(V) == U && error("the ambient scheme of the refinement of X must be X")
        U in patches && error("the ambient scheme of the refinement can not be found in the affine patches")
        if check
          isone(OO(U)(sum([c*g for (c, g) in zip(a, gens(U))]))) || error("the patch $V does not cover $U")
        end
      end
    end
    return new{SpecType, GlueingType, SpecOpenType, RingElemType}(patches, glueings, affine_refinements)
  end
end

### type constructors
covering_type(::Type{T}) where {T<:Spec} = Covering{T, glueing_type(T)}
covering_type(X::Spec) = covering_type(typeof(X))

patches(C::Covering) = C.patches
basic_patches(C::Covering) = C.patches
npatches(C::Covering) = length(C.patches)
glueings(C::Covering) = C.glueings
getindex(C::Covering, i::Int) = C.patches[i]
getindex(C::Covering, i::Int, j::Int) = glueings(C)[(patches(C)[i], patches(C)[j])]
getindex(C::Covering, X::SpecType, Y::SpecType) where {SpecType<:Spec} = glueings(C)[(X, Y)]
edge_dict(C::Covering) = C.edge_dict

affine_refinements(C::Covering) = C.affine_refinements

function add_affine_refinement!(
    C::Covering, U::SpecOpen; 
    a::Vector{RingElemType}=Vector{poly_type(ambient_type(U))}(), 
    check::Bool=true
  ) where {RingElemType<:RingElem}
  X = ambient(U)
  X in patches(C) || error("ambient scheme not found in the basic patches of the covering")
  if check
    if all(x->iszero(x), a) # in case of default argument, do the computations now
      g = gens(U)
      R = base_ring(OO(X))
      F = FreeMod(OO(X), 1)
      A = MatrixSpace(OO(X), length(g), 1)(g)
      M, inc_M = sub(F, A)
      represents_element(F[1], M) || error("patches of $U do not cover $X")
      a = as_vector(coordinates(F[1], M), length(g))
    end
    isone(OO(X)(sum([c*g for (c, g) in zip(a, gens(U))]))) || error("patches of $U do not cover $X")
  end
  if !haskey(affine_refinements(C), X)
    affine_refinements(C)[X] = [(U, a)]
  else
    push!(affine_refinements(C)[X], (U, a))
  end
  return C
end

# functions for handling sets in coverings

using Infiltrator
function Base.intersect(U::SpecType, V::SpecType, C::Covering) where {SpecType<:Spec}
  U in C || error("first patch not found in covering")
  V in C || error("second patch not found in covering")
  (i, j, k) = indexin(U, C)
  (l, m, n) = indexin(V, C)
  if i == l # U and V are affine opens of the same patch X in C
    X = C[i]
    if X == U == V 
      iso = identity_map(SpecOpen(X))
      return iso, iso, iso, iso
    end
    f = one(base_ring(OO(X))) # For the case where U is already a basic patch
    if j != 0 # In this case, U appears in some refinement
      f = gens(affine_refinements(C)[C[i]][j][1])[k]
    end
    g = one(base_ring(OO(X)))
    if m != 0
      g = gens(affine_refinements(C)[C[l]][m][1])[n]
    end
    W = SpecOpen(X, [f*g])
    isoW = identity_map(W)
    incWtoU = inclusion_morphism(W, SpecOpen(U), check=false)
    incWtoV = inclusion_morphism(W, SpecOpen(V), check=false)
    return isoW, isoW, incWtoU, incWtoV
  else
    G = C[i, l]
    (f, g) = glueing_morphisms(G)
    preimV = preimage(f, V)
    preimU = preimage(g, U)
    WU = intersect(preimV, U)
    WV = intersect(preimU, V)
    isoWUtoWV = restriction(f, WU, WV, check=false)
    isoWVtoWU = restriction(g, WV, WU, check=false)
    incWUtoU = inclusion_morphism(WU, SpecOpen(U), check=false)
    incWVtoV = inclusion_morphism(WV, SpecOpen(V), check=false)
    return isoWUtoWV, isoWVtoWU, incWUtoU, incWVtoV
  end
end

function neighbor_patches(C::Covering, U::Spec)
  gg = glueing_graph(C)
  n = neighbors(gg, C[U])
  return [C[i] for i in n]
end

affine_patch_type(C::Covering) = affine_patch_type(typeof(C))
glueing_type(C::Covering{SpecType, GlueingType, SpecOpenType}) where {SpecType<:Spec, GlueingType<:Glueing, SpecOpenType<:SpecOpen} = GlueingType
affine_patch_type(::Type{Covering{SpecType, GlueingType, SpecOpenType, RingElemType}}) where {SpecType<:Spec, GlueingType<:Glueing, SpecOpenType<:SpecOpen, RingElemType<:RingElem} = SpecType
glueing_type(::Type{Covering{SpecType, GlueingType, SpecOpenType}}) where {SpecType<:Spec, GlueingType<:Glueing, SpecOpenType<:SpecOpen} = GlueingType
open_subset_type(::Type{Covering{R, S, T}}) where {R, S, T} = T
open_subset_type(C::Covering) = open_subset_type(typeof(C))

# TODO: For some reason, the `indexin` method won't work. In the long 
# run, one should probably find out why and fix it. 
function getindex(C::Covering, X::SpecType) where {SpecType<:Spec} 
  for i in 1:length(patches(C))
    X === patches(C)[i] && return i
  end
  error("affine scheme could not be found among the patches")
end

function add_glueing!(C::Covering, G::Glueing)
  (X, Y) = patches(G)
  C.glueings[(X, Y)] = G
  C.glueings[(Y, X)] = inverse(G)
  return C
end

### The default constructor
# Returns a scheme in which every affine patch is only 
# glued to itself via the identity.
function Covering(patches::Vector{SpecType}) where {SpecType<:Spec}
  g = Dict{Tuple{SpecType, SpecType}, glueing_type(SpecType)}()
  for X in patches
    U = SpecOpen(X)
    f = identity_map(U)
    g[X,X] = Glueing(X, X, f, f, check=false)
  end
  return Covering(patches, g, check=false)
end

Covering(X::SpecType) where {SpecType<:Spec} = Covering([X])

function Base.show(io::IO, C::Covering) 
  print(io, 
          "Covering with $(npatches(C)) patch" * 
          (npatches(C) == 1 ? "" : "es"))
#           " and glueing graph")
#   print(io, glueing_graph(C))
end

function Base.in(U::Spec, C::Covering)
  for i in 1:npatches(C)
    U == C[i] && return true
  end
  for i in 1:npatches(C)
    if haskey(affine_refinements(C), C[i])
      V = affine_refinements(C)[C[i]]
      for (V, a) in affine_refinements(C)[C[i]]
        U in affine_patches(V) && return true
      end
    end
  end
  return false
end

function Base.indexin(U::Spec, C::Covering)
  for i in 1:npatches(C)
    U == C[i] && return (i, 0, 0)
  end
  for i in 1:npatches(C)
    if haskey(affine_refinements(C), C[i])
      V = affine_refinements(C)[C[i]]
      for j in 1:length(V)
        (W, a) = V[j]
        for k in 1:length(gens(W))
          U == W[k] && return (i, j, k)
        end
      end
    end
  end
  return (0,0,0)
end


### standard constructors 

function standard_covering(X::ProjectiveScheme{CRT}) where {CRT<:AbstractAlgebra.Ring}
  if has_attribute(X, :standard_covering) 
    return get_attribute(X, :standard_covering)
  end
  CX = affine_cone(X)
  kk = base_ring(X)
  S = homog_poly_ring(X)
  r = fiber_dimension(X)
  U = Vector{affine_patch_type(X)}()
  # TODO: Check that all weights are equal to one. Otherwise the routine is not implemented.
  s = symbols(S)
  for i in 0:r
    R, x = PolynomialRing(kk, [Symbol("("*String(s[k+1])*"//"*String(s[i+1])*")") for k in 0:r if k != i])
    phi = hom(S, R, vcat(gens(R)[1:i], [one(R)], gens(R)[i+1:r]))
    I = ideal(R, phi.(gens(defining_ideal(X))))
    push!(U, Spec(R, I))
  end
  result = Covering(U)
  for i in 1:r
    for j in i+1:r+1
      x = gens(base_ring(OO(U[i])))
      y = gens(base_ring(OO(U[j])))
      f = SpecOpenMor(U[i], x[j-1], 
                      U[j], y[i],
                      vcat([x[k]//x[j-1] for k in 1:i-1],
                           [1//x[j-1]],
                           [x[k-1]//x[j-1] for k in i+1:j-1],
                           [x[k]//x[j-1] for k in j:r],
                           x[r+1:end]),
                      check=false
                     )
      g = SpecOpenMor(U[j], y[i],
                      U[i], x[j-1],
                      vcat([y[k]//y[i] for k in 1:i-1],
                           [y[k+1]//y[i] for k in i:j-2],
                           [1//y[i]],
                           [y[k]//y[i] for k in j:r],
                           y[r+1:end]),
                      check=false
                     )
      add_glueing!(result, Glueing(U[i], U[j], f, g, check=false))
    end
  end
  set_attribute!(X, :standard_covering, result)
  return result
end

function standard_covering(X::ProjectiveScheme{CRT}) where {CRT<:MPolyQuoLocalizedRing}
  if has_attribute(X, :standard_covering) 
    return get_attribute(X, :standard_covering)
  end
  CX = affine_cone(X)
  Y = base_scheme(X)
  L = OO(Y)
  W = localized_ring(L)
  R = base_ring(L)
  kk = coefficient_ring(R)
  S = homog_poly_ring(X)
  r = fiber_dimension(X)
  U = Vector{affine_patch_type(X)}()
  pU = Dict{affine_patch_type(X), morphism_type(affine_patch_type(X), typeof(Y))}()
  # TODO: Check that all weights are equal to one. Otherwise the routine is not implemented.
  s = symbols(S)
  # for each homogeneous variable, set up the chart 
  for i in 0:r
    R_fiber, x = PolynomialRing(kk, [Symbol("("*String(s[k+1])*"//"*String(s[i+1])*")") for k in 0:r if k != i])
    F = Spec(R_fiber)
    ambient_space, pF, pY = product(F, Y)
    fiber_vars = pullback(pF).(gens(R_fiber))
    mapped_polys = [map_coefficients(pullback(pY), f) for f in gens(defining_ideal(X))]
    patch = subscheme(ambient_space, [evaluate(f, vcat(fiber_vars[1:i], [one(OO(ambient_space))], fiber_vars[i+1:end])) for f in mapped_polys])
    push!(U, patch)
    pU[patch] = restrict(pY, patch, Y, check=false)
  end
  result = Covering(U)
  for i in 1:r
    for j in i+1:r+1
      x = gens(base_ring(OO(U[i])))
      y = gens(base_ring(OO(U[j])))
      f = SpecOpenMor(U[i], x[j-1], 
                      U[j], y[i],
                      vcat([x[k]//x[j-1] for k in 1:i-1],
                           [1//x[j-1]],
                           [x[k-1]//x[j-1] for k in i+1:j-1],
                           [x[k]//x[j-1] for k in j:r],
                           x[r+1:end]),
                      check=false
                     )
      g = SpecOpenMor(U[j], y[i],
                      U[i], x[j-1],
                      vcat([y[k]//y[i] for k in 1:i-1],
                           [y[k+1]//y[i] for k in i:j-2],
                           [1//y[i]],
                           [y[k]//y[i] for k in j:r],
                           y[r+1:end]),
                      check=false
                     )
      add_glueing!(result, Glueing(U[i], U[j], f, g, check=false))
    end
  end
  covered_projection = CoveringMorphism(result, Covering(Y), pU)
  set_attribute!(X, :standard_covering, result)
  set_attribute!(X, :covered_projection_to_base, covered_projection)
  return result
end

function glueing_graph(C::Covering) 
  if !isdefined(C, :glueing_graph)
    update_glueing_graph(C)
  end
  return C.glueing_graph
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
  end
  C.edge_dict = edge_dict
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

function disjoint_union(C1::T, C2::T) where {T<:Covering}
  C = Covering(vcat(patches(C1), patches(C2)))
  for (X, Y) in keys(glueings(C1))
    add_glueing!(C, C1[X, Y])
  end
  for (X, Y) in keys(glueings(C2))
    add_glueing!(C, C2[X, Y])
  end
  return C
end

@Markdown.doc """
    CoveringMorphism{SpecType<:Spec, CoveringType<:Covering, SpecMorType<:SpecMor}

A morphism ``f : C → D`` of two coverings. For every patch ``U`` of ``C`` this 
provides a map `f[U']` of type `SpecMorType` from ``U' ⊂ U`` to 
some patch `codomain(f[U])` in `D` for some affine patches ``U'`` covering ``U``.

**Note:** For two affine patches ``U₁, U₂ ⊂ U`` the codomains of `f[U₁]` and `f[U₂]`
do not need to coincide! However, given the glueings in `C` and `D`, all affine maps 
have to coincide on their overlaps.
"""
mutable struct CoveringMorphism{SpecType<:Spec, CoveringType<:Covering, SpecMorType<:SpecMor}
  domain::CoveringType
  codomain::CoveringType
  morphisms::Dict{SpecType, SpecMorType} # on a patch X of the domain covering, this 
                                         # returns the morphism φ : X → Y to the corresponding 
                                         # patch Y of the codomain covering. 

  function CoveringMorphism(
      dom::CoveringType, 
      cod::CoveringType, 
      mor::Dict{SpecType, SpecMorType}; 
      check::Bool=true
    ) where {
             SpecType<:Spec, 
             CoveringType<:Covering, 
             SpecMorType<:SpecMor
            }
    # TODO: check domain/codomain compatibility
    # TODO: if check is true, check that all morphisms glue and that the domain patches 
    # cover the basic patches of `dom`.
    for U in keys(mor)
      U in dom || error("patch $U of the map not found in domain")
      codomain(mor[U]) in cod || error("codomain patch not found")
    end
    # check that the whole domain is covered
    for U in basic_patches(dom)
      if !haskey(mor, U)
        !haskey(affine_refinements(dom), U) || error("patch $U of the domain not covered")
        found = false
        for (V, a) in affine_refinements(dom)[U] 
          all(x->(haskey(mor, x)), affine_patches(V)) && (found = true)
        end
        !found && error("patch $U of the domain not covered")
      end
    end
    return new{SpecType, CoveringType, SpecMorType}(dom, cod, mor)
  end
end

### type getters
covering_type(C::CoveringMorphism{R, S, T}) where {R, S, T} = S
covering_type(::Type{CoveringMorphism{R, S, T}}) where {R, S, T} = S
affine_patch_type(C::CoveringMorphism{R, S, T}) where {R, S, T} = R
affine_patch_type(::Type{CoveringMorphism{R, S, T}}) where {R, S, T} = R

morphism_type(C::Covering{SpecType, GlueingType, SpecOpenType}) where {SpecType<:Spec, GlueingType<:Glueing, SpecOpenType<:SpecOpen} = CoveringMorphism{SpecType, Covering{SpecType, GlueingType, SpecOpenType}, morphism_type(SpecType, SpecType)}
morphism_type(::Type{Covering{SpecType, GlueingType, SpecOpenType}}) where {SpecType<:Spec, GlueingType<:Glueing, SpecOpenType<:SpecOpen} = CoveringMorphism{SpecType, Covering{SpecType, GlueingType, SpecOpenType}, morphism_type(SpecType, SpecType)}

domain(f::CoveringMorphism) = f.domain
codomain(f::CoveringMorphism) = f.codomain
getindex(f::CoveringMorphism, U::Spec) = f.morphisms[U]
morphisms(f::CoveringMorphism) = f.morphisms

function compose(f::T, g::T) where {T<:CoveringMorphism}
  domain(g) == codomain(f) || error("morphisms can not be composed")
  morphism_dict = Dict{affine_patch_type(f), morphism_type(affine_patch_type(f), affine_patch_type(f))}()
  for U in patches(domain(f))
    morphism_dict[U] = compose(f[U], g[codomain(f[U])])
  end
  return CoveringMorphism(domain(f), codomain(g), morphism_dict)
end

@Markdown.doc """
    mutable struct CoveredScheme{
      CoveringType<:Covering, 
      CoveringMorphismType<:CoveringMorphism
    }

A covered scheme ``X`` given by means of at least one covering 
of type `CoveringType`. 

A scheme may posess several coverings which are partially ordered 
by refinement. Such refinements are special instances of `CoveringMorphism`

    ρ : C1 → C2

where for each patch ``U`` in `C1` the inclusion map ``ρ[U] : U → V`` 
into the corresponding patch ``V`` of `C2` is an open embedding for which 
both ``𝒪(U)`` and ``𝒪(V)`` have the same `base_ring` (so that they can be 
canonically compared). 
"""
@attributes mutable struct CoveredScheme{CoveringType<:Covering, CoveringMorphismType<:CoveringMorphism}
  coverings::Vector{CoveringType}
  refinements::Dict{Tuple{CoveringType, CoveringType}, CoveringMorphismType}
  refinement_graph::Graph{Directed}

  default_covering::CoveringType

  function CoveredScheme(coverings::Vector{CoveringType}, refinements::Dict{Tuple{CoveringType, CoveringType}, CoveringMorphismType}) where {CoveringType<:Covering, CoveringMorphismType<:CoveringMorphism}
    # TODO: Check whether the refinements form a connected graph.
    X = new{CoveringType, CoveringMorphismType}(coverings, refinements)
    X.default_covering = X.coverings[1]
    return X
  end
end

### type getters
covering_type(X::CoveredScheme{S, T}) where {S, T} = S
covering_type(::Type{CoveredScheme{S, T}}) where {S, T} = S
covering_morphism_type(X::CoveredScheme{S, T}) where {S, T} = T
covering_morphism_type(::Type{CoveredScheme{S, T}}) where {S, T} = T
affine_patch_type(X::CoveredSchemeType) where {CoveredSchemeType<:CoveredScheme} = affine_patch_type(covering_type(CoveredSchemeType))
affine_patch_type(::Type{CoveredSchemeType}) where {CoveredSchemeType<:CoveredScheme} = affine_patch_type(covering_type(CoveredSchemeType))

### type constructors
covered_scheme_type(::Type{T}) where {T<:Spec} = CoveredScheme{covering_type(T), morphism_type(covering_type(T))}
covered_scheme_type(X::Spec) = covered_scheme_type(typeof(X))

covered_scheme_type(::Type{T}) where {T<:ProjectiveScheme} = covered_scheme_type(affine_patch_type(P))
covered_scheme_type(P::ProjectiveScheme) = covered_scheme_type(typeof(P))

### getter methods
coverings(X::CoveredScheme) = X.coverings
refinements(X::CoveredScheme) = X.refinements
getindex(X::CoveredScheme, C::CoveringType, D::CoveringType) where {CoveringType<:Covering} = X.refinements[(C, D)]
setindex(X::CoveredScheme, f::CoveringMorphismType, C::CoveringType, D::CoveringType) where {CoveringMorphismType<:CoveringMorphism, CoveringType<:Covering} = X.refinements[(C, D)]
default_covering(X::CoveredScheme) = X.default_covering
getindex(X::CoveredScheme, i::Int) = coverings(X)[i]

patches(X::CoveredScheme) = patches(default_covering(X))
glueings(X::CoveredScheme) = glueings(default_covering(X))

function getindex(X::CoveredScheme, C::Covering)
  for i in 1:length(coverings(X))
    C == coverings(X)[i] && return i
  end
  error("covering not listed")
end

getindex(X::CoveredScheme, i::Int, j::Int) = X[X[i], X[j]]


set_name!(X::CoveredScheme, name::String) = set_attribute!(X, :name, name)
name_of(X::CoveredScheme) = get_attribute(X, :name)::String
has_name(X::CoveredScheme) = has_attribute(X, :name)

function dim(X::CoveredScheme) 
  if !has_attribute(X, :dim)
    d = -1
    is_equidimensional=true
    for U in patches(default_covering(X))
      e = dim(U)
      if e > d
        d == -1 || (is_equidimensional=false)
        d = e
      end
    end
    set_attribute!(X, :dim, d)
    if !is_equidimensional
      # the above is not an honest check for equidimensionality,
      # because in each chart the output of `dim` is only the 
      # supremum of all components. Thus we can only infer 
      # non-equidimensionality in case this is already visible
      # from comparing the diffent charts
      set_attribute(X, :is_equidimensional, false)
    end
  end
  return get_attribute(X, :dim)::Int
end

function set_default_covering!(X::CoveredScheme, C::Covering) 
  C in coverings(X) || error("covering is not listed")
  X.default_covering = C
  return X
end

### constructors 

function CoveredScheme(C::Covering)
  refinements = Dict{Tuple{typeof(C), typeof(C)}, morphism_type(C)}()
  X = CoveredScheme([C], refinements)
  set_attribute!(X, :seed_covering, C)
  return X
end

CoveredScheme(X::Spec) = CoveredScheme(Covering(X))

# construct the empty covered scheme over the ring R
function empty_covered_scheme(R::RT) where {RT<:AbstractAlgebra.Ring}
  return CoveredScheme(empty_spec(R))
end


function Base.show(io::IO, X::CoveredScheme)
  if has_name(X)
    print(io, name_of(X))
    return
  end
  print(io, "covered scheme with $(npatches(default_covering(X))) affine patches in its default covering")
end

_compose_along_path(X::CoveredScheme, p::Vector{Int}) = _compose_along_path(X, [X[i] for i in p])

function _compose_along_path(X::CoveredScheme, p::Vector{CoveringType}) where {CoveringType<:Covering}
  root = pop!(p)
  next = pop!(p)
  mor = X[next, root]
  while length(p) > 0
    leaf = pop!(p)
    mor = compose(X[leaf, next], mor)
    next = leaf
  end
  X[leaf, root] = mor
  add_edge!(refinement_graph(X), X[leaf], X[root])
  return mor
end

# TODO: Replace by the polymake routines, once provided!
function find_common_root(G::Graph{Directed}, i::Int, j::Int)
  p = [i]
  Ni = neighbors(G, i)
  while length(Ni) > 0
    push!(p, Ni[1])
    Ni = neigbors(G, Ni[1])
  end
  q = [j] 
  Nj = neighbors(G, j)
  while length(Nj) > 0
    push!(p, Nj[1])
    Nj = neigbors(G, Nj[1])
  end
  last(p) == last(q) || error("no common root found")
  return last(p), p, q
end

@Markdown.doc """
    common_refinement(X::CoveredScheme, C1::T, C2::T) where {T<:Covering}

Given two coverings of ``X``, return a triple `(C_new, f, g)` consisting 
of a common refinement `C_new` of `C1` and `C2` and the refinement morphisms 
`f : C_new → C1` and `g : C_new → C2`.
"""
function common_refinement(X::CoveredScheme, C1::T, C2::T) where {T<:Covering}
  # shortcut for the trivial cases
  C1 == C2 && return (C1, identity_map(C1), identity_map(C1))

  # find the minimal common root using the refinement graph
  r, p1, p2 = find_common_root(refinement_graph(X), X[C1], X[C2])

  # if one covering sits strictly on top of the other, take the shortcut
  if length(p1) == 0
    return (C2, identity_map(C1), _compose_along_path(X, p2))
  end
  if length(p2) == 0
    return (C1, _compose_along_path(X, p1), identity_map(C2))
  end
  
  # now we may assume that neither one of the coverings is contained in the other
  C0 = X[r]
  f = _compose_along_path(X, p1)
  g = _compose_along_path(X, p2)

  # prepare for the common refinement
  new_patches = Vector{affine_patch_type(X)}()
  inc1 = Dict{affine_patch_type(X), morphism_type(affine_patch_type(X))}()
  inc2 = Dict{affine_patch_type(X), morphism_type(affine_patch_type(X))}()
  inc0 = Dict{affine_patch_type(X), morphism_type(affine_patch_type(X))}()
  for U in patches(C1)
    W = codomain(f[U])
    V_candidates = [V for V in patches(C2) if codomain(g[V]) == W]

    # first try to find a patch in C2 which fully includes U
    patch_found = false
    while length(V_candidates) > 0
      V = pop!(V_candidates) 
      if issubset(U, V) 
        inc1[U] = identity_map(U)
        inc2[U] = inclusion_map(U, V)
        inc0[U] = f[U]
        push!(new_patches, U)
        patch_found = true
        break
      end
    end
    patch_found && break

    # this is the worst case where there is no patch in C2 containing U.
    for V in V_candidates
      UV = intersect(U, V)
      inc1[UV] = inclusion_map(UV, U)
      inc2[UV] = inclusion_map(UV, V)
      inc0[UV] = inclusion_map(UV, W)
      push!(new_patches, UV)
    end
  end
  
  # cook up the glueings for the new patches from those in the common root.
  new_glueings = Dict{Tuple{affine_patch_type(X), affine_patch_type(X)}, glueing_type(affine_patch_type(X))}()
  for (W1, W2) in keys(glueings(C0))
    U_patches = [U for U in new_patches if codomain(inc0[U]) == W1]
    V_patches = [V for V in new_patches if codomain(inc0[V]) == W2]
    for U in U_patches
      for V in V_patches
        new_glueings[(U, V)] = restriction(C0[W1, W2], U, V)
      end
    end
  end

  C_new = Covering(new_patches, new_glueings)
  f = CoveringMorphism(C_new, C1, inc1, check=true) # set to false after debugging
  g = CoveringMorphism(C_new, C2, inc2, check=true)
  h = CoveringMorphism(C_new, C0, inc0, check=true)
  X[C_new, C1] = f
  X[C_new, C2] = g
  X[C_new, C0] = h
  add_edge!(refinement_graph(X), X[C_new], X[C1])
  add_edge!(refinement_graph(X), X[C_new], X[C2])
  add_edge!(refinement_graph(X), X[C_new], X[C0])
  return (C_new, f, g)
end


@attributes mutable struct CoveredSchemeMorphism{
    DomainType<:CoveredScheme, 
    CodomainType<:CoveredScheme, 
    CoveringMorphismType<:CoveringMorphism
  }
  X::DomainType
  Y::CodomainType
  f::CoveringMorphismType

  function CoveredSchemeMorphism(
      X::DomainType, 
      Y::CodomainType, 
      f::CoveringMorphismType;
      check::Bool=true
    ) where {
             DomainType<:CoveredScheme, 
             CodomainType<:CoveredScheme, 
             CoveringMorphismType<:CoveringMorphism
            }
    domain(f) in coverings(X) || error("covering not found in domain")
    codomain(f) in coverings(Y) || error("covering not found in codomain")
    return new{DomainType, CodomainType, CoveringMorphismType}(X, Y, f)
  end
end

domain(f::CoveredSchemeMorphism) = f.X
codomain(f::CoveredSchemeMorphism) = f.Y
covering_morphism(f::CoveredSchemeMorphism) = f.f
domain_covering(f::CoveredSchemeMorphism) = domain(covering_morphism(f))
codomain_covering(f::CoveredSchemeMorphism) = codomain(covering_morphism(f))
getindex(f::CoveredSchemeMorphism, U::Spec) = covering_morphism(f)[U]

function simplify(C::Covering)
  n = npatches(C)
  SpecType = affine_patch_type(C)
  new_patches = [simplify(X) for X in patches(C)]
  GD = glueings(C)
  new_glueings = Dict{Tuple{SpecType, SpecType}, glueing_type(SpecType)}()
  for (X, Y) in keys(GD)
    Xsimp, iX, jX = new_patches[C[X]]
    Ysimp, iY, jY = new_patches[C[Y]]
    G = GD[(X, Y)]
    U, V = glueing_domains(G)
    f, g = glueing_morphisms(G)
    Usimp = SpecOpen(Xsimp, lifted_numerator.(pullback(iX).(gens(U))), check=false)
    Vsimp = SpecOpen(Ysimp, lifted_numerator.(pullback(iY).(gens(V))), check=false)

    fsimp = SpecOpenMor(Usimp, Vsimp, 
                        [
                         compose(
                                 compose(
                                         restrict(iX, Usimp[k], U[k], check=false), 
                                         f[k]),
                                 jY)
                         for k in 1:npatches(Usimp)],
                        check=false
                       )
    gsimp = SpecOpenMor(Vsimp, Usimp, 
                        [
                         compose(
                                 compose(
                                         restrict(iY, Vsimp[k], V[k], check=false), 
                                         g[k]),
                                 jX)
                         for k in 1:npatches(Vsimp)],
                        check=false
                       )
    new_glueings[(Xsimp, Ysimp)] = Glueing(Xsimp, Ysimp, fsimp, gsimp, check=false)
  end
  iDict = Dict{SpecType, morphism_type(SpecType, SpecType)}()
  jDict = Dict{SpecType, morphism_type(SpecType, SpecType)}()
  for i in 1:length(new_patches)
    iDict[new_patches[i][1]] = new_patches[i][2]
    jDict[C[i]] = new_patches[i][3]
  end
  Cnew = Covering([ U for (U, _, _) in new_patches], new_glueings, check=false)
  i_cov_mor = CoveringMorphism(Cnew, C, iDict)
  j_cov_mor = CoveringMorphism(C, Cnew, jDict)
  return Cnew, i_cov_mor, j_cov_mor
end



function simplify!(X::CoveredScheme)
  C = default_covering(X)
  Csimp, i, j = simplify(C)
  push!(coverings(X), Csimp)
  refinements(X)[(C, Csimp)] = j
  refinements(X)[(Csimp, C)] = i
  return (X, Csimp)
end

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
  result = Vector{affine_patch_type(C)}()
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

Base.eltype(C::Covering) = affine_patch_type(C)

### Miscellaneous helper routines
function as_vector(v::SRow{T}, n::Int) where {T<:RingElem}
  R = base_ring(v)
  result = elem_type(R)[zero(R) for i in 1:n]
  for (i, a) in v
    result[i] = a
  end
  return result
end
