export patches, npatches, glueings, add_glueing!, standard_covering, glueing_graph, update_glueing_graph, transition_graph, edge_dict, disjoint_union, neighbor_patches, affine_refinements, add_affine_refinement!, all_patches

export fill_transitions!

export affine_patch_type, glueing_type

export morphism_type, morphisms

export empty_covered_scheme
export coverings, refinements, default_covering, set_name!, name_of, has_name, dim
export covering_type, covering_morphism_type, affine_patch_type, covered_scheme_type

export domain, codomain, covering_morphism

export simplify

########################################################################
# Methods for Covering                                                 #
########################################################################
### type getters 
base_ring_type(::Type{Covering{T}}) where {T} = T
base_ring_type(C::Covering) = base_ring_type(typeof(C))

### type constructors
#covering_type(::Type{T}) where {T<:Spec} = Covering{T, glueing_type(T)}
#covering_type(X::Spec) = covering_type(typeof(X))

### essential getters
patches(C::Covering) = C.patches
basic_patches(C::Covering) = C.patches
npatches(C::Covering) = length(C.patches)
glueings(C::Covering) = C.glueings
getindex(C::Covering, i::Int) = C.patches[i]
getindex(C::Covering, i::Int, j::Int) = glueings(C)[(patches(C)[i], patches(C)[j])]
getindex(C::Covering, X::AbsSpec, Y::AbsSpec) = glueings(C)[(X, Y)]
#edge_dict(C::Covering) = C.edge_dict

affine_refinements(C::Covering) = C.affine_refinements

#function add_affine_refinement!(
#    C::Covering, U::SpecOpen; 
#    a::Vector{RingElemType}=as_vector(coordinates(one(OO(ambient(U))), 
#                                                  ideal(OO(ambient(U)), 
#                                                        OO(ambient(U)).(gens(U)))), 
#                                      length(gens(U))),
#    check::Bool=true
#  ) where {RingElemType<:RingElem}
#  X = ambient(U)
#  @show typeof(OO(X))
#  X in patches(C) || error("ambient scheme not found in the basic patches of the covering")
#  if check
#    isone(OO(X)(sum([c*g for (c, g) in zip(a, gens(U))]))) || error("patches of $U do not cover $X")
#  end
#  if !haskey(affine_refinements(C), X)
#    affine_refinements(C)[X] = [(U, a)]
#  else
#    push!(affine_refinements(C)[X], (U, a))
#  end
#  return C
#end

# functions for handling sets in coverings

function intersect_in_covering(U::AbsSpec, V::AbsSpec, C::Covering)
  U in C || error("first patch not found in covering")
  V in C || error("second patch not found in covering")
  (i, j, k) = indexin(U, C)
  (l, m, n) = indexin(V, C)
  if i == l # U and V are affine opens of the same patch X in C
    error("affine refinements of coverings not implemented at the moment")
#    X = C[i]
#    if X === U === V
#      iso = identity_map(X)
#      return iso, iso, iso, iso
#    end
#    f = one(ambient_ring(X)) # For the case where U is already a basic patch
#    if j != 0 # In this case, U appears in some refinement
#      f = gens(affine_refinements(C)[X][j][1])[k]
#    end
#    g = one(ambient_ring(X))
#    if m != 0
#      g = gens(affine_refinements(C)[X][m][1])[n]
#    end
#    W = PrincipalOpenSubset(X, [f*g])
#    isoW = identity_map(W)
#    incWtoU = inclusion_morphism(W, U, check=false)
#    incWtoV = inclusion_morphism(W, V, check=false)
#    return isoW, isoW, incWtoU, incWtoV
  else
    G = C[i, l]
    (f, g) = glueing_morphisms(G)
    preimV = preimage(f, V)
    preimU = preimage(g, U)
    WU = intersect(preimV, U)
    WV = intersect(preimU, V)
    isoWUtoWV = restrict(f, WU, WV, check=false)
    isoWVtoWU = restrict(g, WV, WU, check=false)
    incWUtoU = inclusion_morphism(WU, U, check=false)
    incWVtoV = inclusion_morphism(WV, V, check=false)
    return isoWUtoWV, isoWVtoWU, incWUtoU, incWVtoV
  end
end

function neighbor_patches(C::Covering, U::AbsSpec)
  gg = glueing_graph(C)
  n = neighbors(gg, C[U])
  return [C[i] for i in n]
end

#affine_patch_type(C::Covering) = affine_patch_type(typeof(C))
#glueing_type(C::Covering{SpecType, GlueingType, SpecOpenType}) where {SpecType<:Spec, GlueingType<:Glueing, SpecOpenType<:SpecOpen} = GlueingType
#affine_patch_type(::Type{Covering{SpecType, GlueingType, SpecOpenType, RingElemType}}) where {SpecType<:Spec, GlueingType<:Glueing, SpecOpenType<:SpecOpen, RingElemType<:RingElem} = SpecType
#glueing_type(::Type{Covering{SpecType, GlueingType, SpecOpenType}}) where {SpecType<:Spec, GlueingType<:Glueing, SpecOpenType<:SpecOpen} = GlueingType
#open_subset_type(::Type{Covering{R, S, T}}) where {R, S, T} = T
#open_subset_type(C::Covering) = open_subset_type(typeof(C))

# TODO: For some reason, the `indexin` method won't work. In the long 
# run, one should probably find out why and fix it. 
function getindex(C::Covering, X::AbsSpec)
  for i in 1:length(patches(C))
    X === patches(C)[i] && return i
  end
  error("affine scheme could not be found among the patches")
end

function add_glueing!(C::Covering, G::AbsGlueing)
  (X, Y) = patches(G)
  C.glueings[(X, Y)] = G
  C.glueings[(Y, X)] = inverse(G)
  return C
end

### The default constructor
# Returns a scheme in which every affine patch is only 
# glued to itself via the identity.
function Covering(patches::Vector{<:AbsSpec})
  g = IdDict{Tuple{AbsSpec, AbsSpec}, AbsGlueing}()
  for X in patches
    U = PrincipalOpenSubset(X)
    f = identity_map(U)
    g[X,X] = SimpleGlueing(X, X, f, f, check=false)
  end
  return Covering(patches, g, check=false)
end

Covering(X::AbsSpec) = Covering([X])

empty_covering(kk::Ring) = Covering(kk)

function Base.show(io::IO, C::Covering) 
  print(io, 
          "Covering with $(npatches(C)) patch" * 
          (npatches(C) == 1 ? "" : "es"))
#           " and glueing graph")
#   print(io, glueing_graph(C))
end

function Base.in(U::AbsSpec, C::Covering)
  for i in 1:npatches(C)
    U === C[i] && return true
  end
  return false
  ### Affine refinements not implemented at the moment!
#  for i in 1:npatches(C)
#    if haskey(affine_refinements(C), C[i])
#      V = affine_refinements(C)[C[i]]
#      for (V, a) in affine_refinements(C)[C[i]]
#        U in affine_patches(V) && return true
#      end
#    end
#  end
#  return false
end

function Base.indexin(U::AbsSpec, C::Covering)
  for i in 1:npatches(C)
    U === C[i] && return (i, 0, 0)
  end
  return (0,0,0)
  ### Affine refinements not implemented at the moment!
#  for i in 1:npatches(C)
#    if haskey(affine_refinements(C), C[i])
#      V = affine_refinements(C)[C[i]]
#      for j in 1:length(V)
#        (W, a) = V[j]
#        for k in 1:length(gens(W))
#          U === W[k] && return (i, j, k)
#        end
#      end
#    end
#  end
#  return (0,0,0)
end


### standard constructors 

@attr function standard_covering(X::ProjectiveScheme{CRT}) where {CRT<:AbstractAlgebra.Ring}
  CX = affine_cone(X)
  kk = base_ring(X)
  S = ambient_ring(X)
  r = fiber_dimension(X)
  U = Vector{AbsSpec}()
  # TODO: Check that all weights are equal to one. Otherwise the routine is not implemented.
  s = symbols(S)
  for i in 0:r
    R, x = PolynomialRing(kk, [Symbol("("*String(s[k+1])*"//"*String(s[i+1])*")") for k in 0:r if k != i])
    phi = hom(S, R, vcat(gens(R)[1:i], [one(R)], gens(R)[i+1:r]))
    I = ideal(R, phi.(gens(defining_ideal(X))))
    push!(U, Spec(quo(R, I)[1]))
  end
  result = Covering(U)
  for i in 1:r
    for j in i+1:r+1
      x = gens(base_ring(OO(U[i])))
      y = gens(base_ring(OO(U[j])))
      Ui = PrincipalOpenSubset(U[i], OO(U[i])(x[j-1]))
      Uj = PrincipalOpenSubset(U[j], OO(U[j])(y[i]))
      f = SpecMor(Ui, Uj,
                      vcat([x[k]//x[j-1] for k in 1:i-1],
                           [1//x[j-1]],
                           [x[k-1]//x[j-1] for k in i+1:j-1],
                           [x[k]//x[j-1] for k in j:r],
                           x[r+1:end]),
                      check=false
                     )
      g = SpecMor(Uj, Ui,
                      vcat([y[k]//y[i] for k in 1:i-1],
                           [y[k+1]//y[i] for k in i:j-2],
                           [1//y[i]],
                           [y[k]//y[i] for k in j:r],
                           y[r+1:end]),
                      check=false
                     )
      add_glueing!(result, SimpleGlueing(U[i], U[j], f, g, check=false))
    end
  end
  return result
end

@attr function standard_covering(X::ProjectiveScheme{CRT}) where {CRT<:Union{<:MPolyQuoLocalizedRing, <:MPolyLocalizedRing, <:MPolyRing, <:MPolyQuo}}
  CX = affine_cone(X)
  Y = base_scheme(X)
  R = ambient_ring(Y)
  kk = coefficient_ring(R)
  S = ambient_ring(X)
  r = fiber_dimension(X)
  U = Vector{AbsSpec}()
  pU = IdDict{AbsSpec, AbsSpecMor}()
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
    C.edge_dict = edge_dict
  end
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

### type getters
# Not required at the moment; should eventually be deleted.
#base_morphism_type(::Type{T}) where {DT, CT, BMT, T<:CoveringMorphism{DT, CT, BMT}} = BMT
#base_morphism_type(C::Covering) = base_morphism_type(typeof(C))

#domain_type(::Type{T}) where {DT, CT, BMT, T<:CoveringMorphism{DT, CT, BMT}} = DT
#domain_type(C::Covering) = domain_type(typeof(C))

#codomain_type(::Type{T}) where {DT, CT, BMT, T<:CoveringMorphism{DT, CT, BMT}} = CT
#codomain_type(C::Covering) = codomain_type(typeof(C))

########################################################################
# Methods for CoveringMorphism                                         #
########################################################################

#covering_type(C::CoveringMorphism{R, S, T}) where {R, S, T} = S
#covering_type(::Type{CoveringMorphism{R, S, T}}) where {R, S, T} = S
#affine_patch_type(C::CoveringMorphism{R, S, T}) where {R, S, T} = R
#affine_patch_type(::Type{CoveringMorphism{R, S, T}}) where {R, S, T} = R

#morphism_type(C::Covering{SpecType, GlueingType, SpecOpenType}) where {SpecType<:Spec, GlueingType<:Glueing, SpecOpenType<:SpecOpen} = CoveringMorphism{SpecType, Covering{SpecType, GlueingType, SpecOpenType}, morphism_type(SpecType, SpecType)}
#morphism_type(::Type{Covering{SpecType, GlueingType, SpecOpenType}}) where {SpecType<:Spec, GlueingType<:Glueing, SpecOpenType<:SpecOpen} = CoveringMorphism{SpecType, Covering{SpecType, GlueingType, SpecOpenType}, morphism_type(SpecType, SpecType)}

domain(f::CoveringMorphism) = f.domain
codomain(f::CoveringMorphism) = f.codomain
getindex(f::CoveringMorphism, U::AbsSpec) = f.morphisms[U]
morphisms(f::CoveringMorphism) = f.morphisms

function compose(f::CoveringMorphism, g::CoveringMorphism)
  domain(g) === codomain(f) || error("morphisms can not be composed")
  morphism_dict = IdDict{AbsSpec, AbsSpecMor}()
  for U in patches(domain(f))
    morphism_dict[U] = compose(f[U], g[codomain(f[U])])
  end
  return CoveringMorphism(domain(f), codomain(g), morphism_dict)
end

########################################################################
# Interface for AbsCoveredScheme                                       #
########################################################################
### type getters
base_ring_type(::Type{T}) where {BRT, T<:AbsCoveredScheme{BRT}} = BRT
base_ring_type(X::AbsCoveredScheme) = base_ring_type(typeof(X))

### essential getters for the interface
base_ring(X::AbsCoveredScheme) = base_ring(underlying_scheme(X))

@Markdown.doc """
    coverings(X::AbsCoveredScheme)

Returns the available coverings for ``X``.
"""
function coverings(X::AbsCoveredScheme) ::Vector{<:Covering}
  return coverings(underlying_scheme(X))
end

@Markdown.doc """
    default_covering(X::AbsCoveredScheme)::Covering

Returns the default covering for ``X``.
"""
function default_covering(X::AbsCoveredScheme)
  return default_covering(underlying_scheme(X))::Covering
end

@Markdown.doc """
    patches(X::AbsCoveredScheme) = patches(default_covering(X))

Returns the affine patches in the `default_covering` of ``X``.
"""
patches(X::AbsCoveredScheme) = patches(default_covering(X))

### Forwarding of the non-documented getters
Base.in(U::AbsSpec, X::AbsCoveredScheme) = (U in underlying_scheme(X))::Bool
getindex(X::AbsCoveredScheme, i::Int) = coverings(underlying_scheme(X))[i]::Covering
getindex(X::AbsCoveredScheme, C::Covering, D::Covering) = getindex(underlying_scheme(X), C, D)::CoveringMorphism
#setindex(X::AbsCoveredScheme, f::CoveringMorphism, C::Covering, D::Covering) = setindex(underlying_scheme(X), f, C, D)::CoveringMorphism
refinements(X::AbsCoveredScheme) = refinements(underlying_scheme(X))::Dict{<:Tuple{<:Covering, <:Covering}, <:CoveringMorphism}
glueings(X::AbsCoveredScheme) = glueings(underlying_scheme(X))::IdDict{<:Tuple{<:AbsSpec, <:AbsSpec}, <:AbsGlueing}


base_ring(X::CoveredScheme) = X.kk

########################################################################
# Methods for CoveredScheme                                            #
########################################################################
### type getters
#covering_type(X::CoveredScheme{S, T}) where {S, T} = S
#covering_type(::Type{CoveredScheme{S, T}}) where {S, T} = S
#covering_morphism_type(X::CoveredScheme{S, T}) where {S, T} = T
#covering_morphism_type(::Type{CoveredScheme{S, T}}) where {S, T} = T
#affine_patch_type(X::CoveredSchemeType) where {CoveredSchemeType<:CoveredScheme} = affine_patch_type(covering_type(CoveredSchemeType))
#affine_patch_type(::Type{CoveredSchemeType}) where {CoveredSchemeType<:CoveredScheme} = affine_patch_type(covering_type(CoveredSchemeType))

### type constructors
#covered_scheme_type(::Type{T}) where {T<:Spec} = CoveredScheme{covering_type(T), morphism_type(covering_type(T))}
#covered_scheme_type(X::Spec) = covered_scheme_type(typeof(X))
#
#covered_scheme_type(::Type{T}) where {T<:ProjectiveScheme} = covered_scheme_type(affine_patch_type(P))
#covered_scheme_type(P::ProjectiveScheme) = covered_scheme_type(typeof(P))

### getter methods
coverings(X::CoveredScheme) = X.coverings
refinements(X::CoveredScheme) = X.refinements
getindex(X::CoveredScheme, C::CoveringType, D::CoveringType) where {CoveringType<:Covering} = X.refinements[(C, D)]
setindex(X::CoveredScheme, f::CoveringMorphismType, C::CoveringType, D::CoveringType) where {CoveringMorphismType<:CoveringMorphism, CoveringType<:Covering} = X.refinements[(C, D)]
default_covering(X::CoveredScheme) = X.default_covering
getindex(X::CoveredScheme, i::Int) = coverings(X)[i]

patches(X::CoveredScheme) = patches(default_covering(X))
glueings(X::CoveredScheme) = glueings(default_covering(X))

function getindex(X::AbsCoveredScheme, C::Covering)
  for i in 1:length(coverings(X))
    C == coverings(X)[i] && return i
  end
  error("covering not listed")
end

getindex(X::AbsCoveredScheme, i::Int, j::Int) = X[X[i], X[j]]

Base.in(U::AbsSpec, X::CoveredScheme) = any(C->(U in C), coverings(X))

set_name!(X::AbsCoveredScheme, name::String) = set_attribute!(X, :name, name)
name_of(X::AbsCoveredScheme) = get_attribute(X, :name)::String
has_name(X::AbsCoveredScheme) = has_attribute(X, :name)

function dim(X::AbsCoveredScheme) 
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

#function set_default_covering!(X::CoveredScheme, C::Covering) 
#  C in coverings(X) || error("covering is not listed")
#  X.default_covering = C
#  return X
#end

### constructors 

function CoveredScheme(C::Covering)
  refinements = Dict{Tuple{Covering, Covering}, CoveringMorphism}()
  X = CoveredScheme([C], refinements)
  set_attribute!(X, :seed_covering, C)
  return X
end

CoveredScheme(X::AbsSpec) = CoveredScheme(Covering(X))

# construct the empty covered scheme over the ring R
function empty_covered_scheme(R::RT) where {RT<:AbstractAlgebra.Ring}
  return CoveredScheme(R)
end

function is_empty(X::AbsCoveredScheme)
  if !isdefined(X, :coverings) 
    return true
  end
  return all(x->isempty(x), all_patches(default_covering(X)))
end



function Base.show(io::IO, X::AbsCoveredScheme)
  if has_name(X)
    print(io, name_of(X))
    return
  end
  print(io, "covered scheme with $(npatches(default_covering(X))) affine patches in its default covering")
end

_compose_along_path(X::CoveredScheme, p::Vector{Int}) = _compose_along_path(X, [X[i] for i in p])

#function _compose_along_path(X::CoveredScheme, p::Vector{CoveringType}) where {CoveringType<:Covering}
#  root = pop!(p)
#  next = pop!(p)
#  mor = X[next, root]
#  while length(p) > 0
#    leaf = pop!(p)
#    mor = compose(X[leaf, next], mor)
#    next = leaf
#  end
#  X[leaf, root] = mor
#  add_edge!(refinement_graph(X), X[leaf], X[root])
#  return mor
#end
#
## TODO: Replace by the polymake routines, once provided!
#function find_common_root(G::Graph{Directed}, i::Int, j::Int)
#  p = [i]
#  Ni = neighbors(G, i)
#  while length(Ni) > 0
#    push!(p, Ni[1])
#    Ni = neigbors(G, Ni[1])
#  end
#  q = [j] 
#  Nj = neighbors(G, j)
#  while length(Nj) > 0
#    push!(p, Nj[1])
#    Nj = neigbors(G, Nj[1])
#  end
#  last(p) == last(q) || error("no common root found")
#  return last(p), p, q
#end
#
#@Markdown.doc """
#    common_refinement(X::CoveredScheme, C1::T, C2::T) where {T<:Covering}
#
#Given two coverings of ``X``, return a triple `(C_new, f, g)` consisting 
#of a common refinement `C_new` of `C1` and `C2` and the refinement morphisms 
#`f : C_new → C1` and `g : C_new → C2`.
#"""
#function common_refinement(X::CoveredScheme, C1::T, C2::T) where {T<:Covering}
#  # shortcut for the trivial cases
#  C1 == C2 && return (C1, identity_map(C1), identity_map(C1))
#
#  # find the minimal common root using the refinement graph
#  r, p1, p2 = find_common_root(refinement_graph(X), X[C1], X[C2])
#
#  # if one covering sits strictly on top of the other, take the shortcut
#  if length(p1) == 0
#    return (C2, identity_map(C1), _compose_along_path(X, p2))
#  end
#  if length(p2) == 0
#    return (C1, _compose_along_path(X, p1), identity_map(C2))
#  end
#  
#  # now we may assume that neither one of the coverings is contained in the other
#  C0 = X[r]
#  f = _compose_along_path(X, p1)
#  g = _compose_along_path(X, p2)
#
#  # prepare for the common refinement
#  new_patches = Vector{affine_patch_type(X)}()
#  inc1 = IdDict{affine_patch_type(X), morphism_type(affine_patch_type(X))}()
#  inc2 = IdDict{affine_patch_type(X), morphism_type(affine_patch_type(X))}()
#  inc0 = IdDict{affine_patch_type(X), morphism_type(affine_patch_type(X))}()
#  for U in patches(C1)
#    W = codomain(f[U])
#    V_candidates = [V for V in patches(C2) if codomain(g[V]) === W]
#
#    # first try to find a patch in C2 which fully includes U
#    patch_found = false
#    while length(V_candidates) > 0
#      V = pop!(V_candidates) 
#      if issubset(U, V) 
#        inc1[U] = identity_map(U)
#        inc2[U] = inclusion_morphism(U, V)
#        inc0[U] = f[U]
#        push!(new_patches, U)
#        patch_found = true
#        break
#      end
#    end
#    patch_found && break
#
#    # this is the worst case where there is no patch in C2 containing U.
#    for V in V_candidates
#      UV = intersect(U, V)
#      inc1[UV] = inclusion_morphism(UV, U)
#      inc2[UV] = inclusion_morphism(UV, V)
#      inc0[UV] = inclusion_morphism(UV, W)
#      push!(new_patches, UV)
#    end
#  end
#  
#  # cook up the glueings for the new patches from those in the common root.
#  new_glueings = IdDict{Tuple{affine_patch_type(X), affine_patch_type(X)}, glueing_type(affine_patch_type(X))}()
#  for (W1, W2) in keys(glueings(C0))
#    U_patches = [U for U in new_patches if codomain(inc0[U]) === W1]
#    V_patches = [V for V in new_patches if codomain(inc0[V]) === W2]
#    for U in U_patches
#      for V in V_patches
#        new_glueings[(U, V)] = restrict(C0[W1, W2], U, V)
#      end
#    end
#  end
#
#  C_new = Covering(new_patches, new_glueings)
#  f = CoveringMorphism(C_new, C1, inc1, check=true) # set to false after debugging
#  g = CoveringMorphism(C_new, C2, inc2, check=true)
#  h = CoveringMorphism(C_new, C0, inc0, check=true)
#  X[C_new, C1] = f
#  X[C_new, C2] = g
#  X[C_new, C0] = h
#  add_edge!(refinement_graph(X), X[C_new], X[C1])
#  add_edge!(refinement_graph(X), X[C_new], X[C2])
#  add_edge!(refinement_graph(X), X[C_new], X[C0])
#  return (C_new, f, g)
#end

########################################################################
# Interface for AbsCoveredSchemeMorphism                               #
########################################################################
### essential getters 
function domain(f::AbsCoveredSchemeMorphism)::CoveredScheme
  return domain(underlying_morphism(f))
end

function codomain(f::AbsCoveredSchemeMorphism)::CoveredScheme
  return codomain(underlying_morphism(f))
end

function covering_morphism(f::AbsCoveredSchemeMorphism)::CoveringMorphism
  return covering_morphism(underlying_morphism(f))
end

### generically derived getters
domain_covering(f::AbsCoveredSchemeMorphism) = domain(covering_morphism(f))
codomain_covering(f::AbsCoveredSchemeMorphism) = codomain(covering_morphism(f))
getindex(f::AbsCoveredSchemeMorphism, U::Spec) = covering_morphism(f)[U]

########################################################################
# Methods for CoveredSchemeMorphism                                    #
########################################################################
domain(f::CoveredSchemeMorphism) = f.X
codomain(f::CoveredSchemeMorphism) = f.Y
covering_morphism(f::CoveredSchemeMorphism) = f.f

lifted_numerator(f::MPolyElem) = f
lifted_numerator(f::MPolyQuoElem) = lift(f)

@Markdown.doc """
    simplify(C::Covering)

Given a covering ``C`` apply `simplify` to all basic affine patches 
in ``C`` and return a triple ``(C', f, g)`` consisting of the 
resulting covering ``C'`` and the identifying isomorphism 
``f : C' ↔ C``.
"""
function simplify(C::Covering)
  n = npatches(C)
  new_patches = [simplify(X) for X in patches(C)]
  GD = glueings(C)
  new_glueings = IdDict{Tuple{AbsSpec, AbsSpec}, AbsGlueing}()
  for (X, Y) in keys(GD)
    Xsimp, iX, jX = new_patches[C[X]]
    Ysimp, iY, jY = new_patches[C[Y]]
    G = GD[(X, Y)]
    new_glueings[(Xsimp, Ysimp)] = restrict(G, jX, jY, check=false)
  end
  iDict = IdDict{AbsSpec, AbsSpecMor}()
  jDict = IdDict{AbsSpec, AbsSpecMor}()
  for i in 1:length(new_patches)
    iDict[new_patches[i][1]] = new_patches[i][2]
    jDict[C[i]] = new_patches[i][3]
  end
  Cnew = Covering([ U for (U, _, _) in new_patches], new_glueings, check=false)
  i_cov_mor = CoveringMorphism(Cnew, C, iDict)
  j_cov_mor = CoveringMorphism(C, Cnew, jDict)
  return Cnew, i_cov_mor, j_cov_mor
end



@Markdown.doc """
    simplify!(X::AbsCoveredScheme)

Apply `simplify` to the `default_covering` of `X` and store the 
resulting `Covering` ``C'`` and its identification with the default 
covering in `X`; return a tuple ``(X, C')``.
"""
function simplify!(X::AbsCoveredScheme)
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
  result = Vector{AbsSpec}()
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

Base.eltype(C::Covering) = AbsSpec

### Miscellaneous helper routines
#function as_vector(v::SRow{T}, n::Int) where {T<:RingElem}
#  R = base_ring(v)
#  result = elem_type(R)[zero(R) for i in 1:n]
#  for (i, a) in v
#    result[i] = a
#  end
#  return result
#end
