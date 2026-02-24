# A morphism from a polynomial ring phi : R[x_1,...,x_n] -> S into any
# (noncommutative) ring is uniquely determined by the restriction of phi to K
# and the images of the x_i under phi.
#
# Thus, the datum for the morphism from a polynomial ring into a ring consists of
#
# - the morphism restricted to the coefficient ring (coeff_map)
# - the images of the generators (img_gens)
#
# To ease the usage we also add the domain and the codomain.
#
# In case there is a canonical embedding R -> S (which happens for example
# if R === S, or R === ZZ and S arbitrary, we store the coefficient ring map
# as nothing.
#
# Evaluating a multivariate polynomial f  under the map is thus either
#
# evaluate((coeff_map(f), img_gens)
#
# or just
#
# evaluate(f, img_gens)
# 
# Note that the evaluate function is designed to work whenever there are
# canonical embeddings in place (aka S(...))
#
# Finally note that the same strategy works for
# - MPolyRing
# - MPolyQuoRing
# - MPolyDecRing

### See `Types.jl` for the declaration of the type.

################################################################################
#
#  Field access
#
################################################################################

domain(f::MPolyAnyMap) = f.domain
codomain(f::MPolyAnyMap) = f.codomain


# Not sure if we want to expose the following function to the user.
# It might be `nothing`. We could return `identity in the `nothing` case.
coefficient_map(f::MPolyAnyMap) = f.coeff_map

_images(f::MPolyAnyMap) = f.img_gens

@doc raw"""
    _maps_variables_to_variables(f::MPolyAnyMap)

A cheap check to see whether `f` is taking the variables of its domain into 
pairwise different variables of its codomain. In the affirmative case 
return `(true, ind)` where `ind` is a `Vector` of the indices so that 
`gens(codomain(f))[ind]` equals the images of the generators. 
Otherwise, return `(false, garbage)`.

Note: For rings different from plain polynomial rings in the codomain 
this is not a rigorous check (which would probably be very expensive), 
but based only a quick look on the representatives!

Do not mutate the second return value.
"""
function _maps_variables_to_variables(f::MPolyAnyMap)
  if ngens(domain(f)) == 0
    return true, f.variable_indices
  end
  return !isempty(f.variable_indices), f.variable_indices
end

# This can be overwritten in order to avoid making the above check a bottleneck.
_cmp_reps(a) = ==(a)

function _assert_has_maps_variables_to_variables!(f::MPolyAnyMap)
  if !isdefined(f, :variable_indices)
    f.variable_indices = __maps_variables_to_variables(_images(f), codomain(f))
  end
end

function __maps_variables_to_variables(img_gens::Vector, C)
  # C is codomain

  # this is a dirty implicit check, that the codomain is of a useful type for
  # everything to make sense
  # proper check would be to see if C isa MPolyRing etc
  if !all(_is_gen, img_gens)
    return Int[]
  end

  if !_allunique(img_gens)
    return Int[]
  end

  l = length(img_gens)
  r = Vector{Int}(undef, l)
  Cgens = gens(C)

  for i in 1:length(img_gens)
    j = findfirst(_cmp_reps(img_gens[i]), Cgens)
    if j isa Nothing
      # image not a variable
      return Int[]
    end
    r[i] = j::Int
  end

  return r
end

################################################################################
#
#  String I/O
#
################################################################################
#
function Base.show(io::IO, ::MIME"text/plain", f::MPolyAnyMap)
  io = pretty(io)
  println(terse(io), f)
  print(io, Indent())
  println(io, "from ", Lowercase(), domain(f))
  println(io, "to ", Lowercase(), codomain(f))
  println(io, Dedent(), "defined by", Indent())
  R = domain(f)
  g = gens(R)
  for i in 1:(ngens(R)-1)
    println(io, g[i], " -> ", f(g[i]))
  end
  print(io, g[end], " -> ", f(g[end]), Dedent())
  # the last print statement must not add a new line
  phi = coefficient_map(f)
  if !(phi isa Nothing)
    println(io)
    println(io, "with map on coefficients")
    print(io, Indent(), phi, Dedent())
  end
end

function Base.show(io::IO, f::MPolyAnyMap)
  io = pretty(io)
  if is_terse(io)
    print(io, "Ring homomorphism")
  else
    print(io, "Hom: ")
    print(terse(io), Lowercase(), domain(f), " -> ")
    print(terse(io), Lowercase(), codomain(f))
  end
end

################################################################################
#
#  Helper
#
################################################################################

# # Since we want to allow specifying images in a "subring", we need to coerce
# if necessary. For example, hom(Qx, Qx, [1, 1]), should work, although
# 1 is not an element of the codomain.
function _coerce(S, img_gens)
  if eltype(img_gens) === elem_type(S)
    return img_gens::Vector{elem_type(S)}
  else
    _img_gens = S.(img_gens)
    @req eltype(_img_gens) === elem_type(S) "Elements cannot be coerced into the codomain"
    return _img_gens::Vector{elem_type(S)}
  end
end

# The following is used when the domain is graded

# if the codomain is graded, the images must be homogeneous?!
# Grading helpers for polynomial ring maps
_isgraded(::NCRing) = false
_isgraded(R::MPolyDecRing) = is_graded(R)
_isgraded(R::MPolyQuoRing) = _isgraded(base_ring(R))

function _check_homo(S::NCRing, images)
  if _isgraded(S)
    for i in images
      @req is_homogeneous(i) "Images must be homogeneous"
    end
  end
end

# Kind of a map, determined by types

_kind_for(::NCRing, ::NCRing) = HomUngraded

_kind_for(::MPolyDecRing, ::NCRing) = HomGradedToUngraded
_kind_for(::NCRing, ::MPolyDecRing) = HomUngradedToGraded
_kind_for(::MPolyDecRing, ::MPolyDecRing) = HomGraded

_kind_for(R::MPolyQuoRing, S::NCRing) = _kind_for(base_ring(R), S)
_kind_for(R::NCRing, S::MPolyQuoRing) = _kind_for(R, base_ring(S))
_kind_for(R::MPolyQuoRing, S::MPolyQuoRing) = _kind_for(base_ring(R), base_ring(S))

_ungraded_kind_for(::NCRing, ::NCRing) = HomUngraded
_ungraded_kind_for(::MPolyDecRing, ::NCRing) = HomGradedToUngraded
_ungraded_kind_for(::NCRing, ::MPolyDecRing) = HomUngradedToGraded
_ungraded_kind_for(::MPolyDecRing, ::MPolyDecRing) = HomUngraded

_ungraded_kind_for(R::MPolyQuoRing, S::NCRing) = _ungraded_kind_for(base_ring(R), S)
_ungraded_kind_for(R::NCRing, S::MPolyQuoRing) = _ungraded_kind_for(R, base_ring(S))
_ungraded_kind_for(R::MPolyQuoRing, S::MPolyQuoRing) = _ungraded_kind_for(base_ring(R), base_ring(S))

hom_kind(f::MPolyAnyMap{D, C, U, V, K}) where {D, C, U, V, K} = K

# Compute degree shift for graded maps.
#
# Condition:
#   degree(img_i) == phi(degree(gen(R,i))) + deg_shift  for all i.
#
# If grading groups differ, phi must be supplied
# except if both groups are the infinite cyclic group
# If groups are equal and phi is omitted, we use hom(GR, GS, gens(GS))

function _is_Z(G)
  return (G isa FinGenAbGroup) && !isfinite(G) &&number_of_generators(G) == 1
end

function _graded_data(R::NCRing, S::NCRing, imgs;
                      grading_group_hom = nothing,
                      degree_shift = nothing)
  GR = grading_group(R)
  GS = grading_group(S)

  phi = grading_group_hom
  if phi === nothing
    if GR == GS
      phi = hom(GR, GS, gens(GS))
    elseif _is_Z(GR) && _is_Z(GS)
      phi = hom(GR, GS, gens(GS))
    else
      error("Source and target have different grading groups; please supply grading_group_hom.")
    end
  else
    @req grading_group_hom isa Map "grading_group_hom must be a Map"
    @req domain(grading_group_hom) === GR "grading_group_hom has wrong domain"
    @req codomain(grading_group_hom) === GS "grading_group_hom has wrong codomain"
    phi = grading_group_hom
  end
  deg_shift = degree_shift
  if deg_shift === nothing
    n = ngens(R)
    if n == 0
      deg_shift = zero(GS)
    else
      deg_shift = degree(imgs[1]) - phi(degree(gen(R, 1)))
    end
  end

  return phi, deg_shift
end



# When evaluating the map F at a polynomial f, we first construct the polynomial
# map_coefficients(coefficient_map(F), f), which is a polynomial over
# codomain(coefficient_map(F)).

_nvars(R::MPolyRing) = nvars(R)

_nvars(R::MPolyQuoRing) = nvars(base_ring(R))

function temp_ring(f::MPolyAnyMap{<:Any, <: Any, <: Map})
  if isdefined(f, :temp_ring)
    return f.temp_ring::mpoly_ring_type(codomain(coefficient_map(f)))
  end

  S, = polynomial_ring(codomain(coefficient_map(f)), _nvars(domain(f)); cached = false)
  f.temp_ring = S
  return S
end

# If the coefficient_map is nothing, we can just evaluate the polynomial directly
function temp_ring(f::MPolyAnyMap{<:Any, <: Any, Nothing})
  return nothing
end

# If the coefficient_map is e.g. a julia function, there is not too much we can
# do, so we do the defensive thing
function temp_ring(f::MPolyAnyMap{<:Any, <: Any})
  return nothing
end

################################################################################
#
#  Kernel
#
################################################################################

function kernel(f::MPolyAnyMap)
  error("Cannot compute this kernel!")
end

################################################################################
#
#  Injectivity
#
################################################################################

function is_injective(F::MPolyAnyMap)
  error("Cannot decide injectivity!")
end

################################################################################
#
#  Surjectivity
#
################################################################################

function is_surjective(F::MPolyAnyMap)
  error("Cannot decide surjectivity!")
end

################################################################################
#
#  Bijectivity
#
################################################################################

function is_bijective(F::MPolyAnyMap)
  error("Cannot decide bijectivity!")
end

################################################################################
#
#  Composition
#
################################################################################

# This is getting difficult, because Map{C, D} does not yield information
# on the type of the domain, codomain

# First consider the case where both coefficient maps are maps in the Map
# sense
function compose(F::MPolyAnyMap{D, C, S}, G::MPolyAnyMap{C, E, U}) where {D, C, E, S <: Map, U <: Map}
  @req codomain(F) === domain(G) "Incompatible (co)domain in composition"
  f = coefficient_map(F)
  g = coefficient_map(G)
  if typeof(codomain(f)) === typeof(domain(g))
    newcoeffmap = compose(f, g)
    return hom(domain(F), codomain(G), newcoeffmap, G.(_images(F)), check=false)
  else
    return Generic.CompositeMap(F, G)
  end
end

# No coefficient maps in the second argument
function compose(F::MPolyAnyMap{D, C, S}, G::MPolyAnyMap{C, E, Nothing}) where {D, C, E, S <: Map}
  @req codomain(F) === domain(G) "Incompatible (co)domain in composition"
  f = coefficient_map(F)
  if typeof(codomain(f)) === typeof(coefficient_ring(domain(G)))
    return hom(domain(F), codomain(G), f, G.(_images(F)), check=false)
  elseif typeof(codomain(f)) === typeof(domain(G))
    new_coeff_map = compose(f, G)
    return hom(domain(F), codomain(G), new_coeff_map, G.(_images(F)), check=false)
  else
    return Generic.CompositeMap(F, G)
  end
end

# No coefficient maps in the first argument
function compose(F::MPolyAnyMap{D, C, Nothing}, G::MPolyAnyMap{C, E, S}) where {D, C, E, S <: Map}
  @req codomain(F) === domain(G) "Incompatible (co)domain in composition"
  g = coefficient_map(G)
  if domain(g) === coefficient_ring(domain(F))
    return hom(domain(F), codomain(G), g, G.(_images(F)), check=false)
  else
    new_coeff_map = MapFromFunc(coefficient_ring(domain(F)), codomain(g), x->g(domain(g)(x)))
    return hom(domain(F), codomain(G), new_coeff_map, G.(_images(F)), check=false)
  end
end

# No coefficient maps in both maps
#function compose(F::MPolyAnyMap{D, C, Nothing}, G::MPolyAnyMap{C, E, Nothing}) where {D, C, E}
#  @req codomain(F) === domain(G) "Incompatible (co)domain in composition"
#  return hom(domain(F), codomain(G), G.(_images(F)), check=false)
#end

# No coefficient maps in both maps; now respect grading kinds.
function compose(F::MPolyAnyMap{D, C, Nothing, V1, K1},
                 G::MPolyAnyMap{C, E, Nothing, V2, K2}) where {D, C, E, V1, V2,
                                                              K1 <: MPolyHomKind,
                                                              K2 <: MPolyHomKind}
  @req codomain(F) === domain(G) "Incompatible (co)domain in composition"

  R = domain(F)
  T = codomain(G)

  imgs = G.(_images(F))

  if K1 === HomGraded && K2 === HomGraded
    h = MPolyAnyMap(HomGraded, R, T, nothing, imgs)

    phi_F = get_attribute(F, :grading_group_hom, nothing)
    dsh_F = get_attribute(F, :degree_shift, nothing)
    phi_G = get_attribute(G, :grading_group_hom, nothing)
    dsh_G = get_attribute(G, :degree_shift, nothing)

    @req phi_F !== nothing && dsh_F !== nothing &&
         phi_G !== nothing && dsh_G !== nothing "Missing grading data on graded maps."

    phi_H = compose(phi_G, phi_F)
    dsh_H = phi_G(dsh_F) + dsh_G

    set_attribute!(h, :grading_group_hom => phi_H)
    set_attribute!(h, :degree_shift => dsh_H)

    return h
  end

  Kout = _ungraded_kind_for(R, T)
  return MPolyAnyMap(Kout, R, T, nothing, imgs)
end


# Julia functions in both maps
function compose(F::MPolyAnyMap{D, C, <: Function}, G::MPolyAnyMap{C, E, <: Function}) where {D, C, E}
  @req codomain(F) === domain(G) "Incompatible (co)domain in composition"
  b = coefficient_map(F)(one(coefficient_ring(domain(F))))
  if parent(b) === domain(G)
    return hom(domain(F), codomain(G), x -> G(coefficient_map(F)(x)), G.(_images(F)), check=false)
  elseif parent(b) === coefficient_ring(domain(G))
    return hom(domain(F), codomain(G), x -> coefficient_map(G)(coefficient_map(F)(x)), G.(_images(F)), check=false)
  else
    error("coefficient map is not admissible")
  end
end

# Now compose with arbitrary maps

# I technically cannot do the Nothing version

# I can only do the Map version of the coefficient map has codomain C
function compose(F::MPolyAnyMap{D, C, <: Map, <: Any}, G::S) where {D, C, S <: Map{C, <: Any}}
  @req codomain(F) === domain(G) "Incompatible (co)domain in composition"
  f = coefficient_map(F)
  if typeof(codomain(f)) === C
    newcoeffmap = compose(f, G)
    return hom(domain(F), codomain(G), newcoeffmap, G.(_images(F)), check=false)
  else
    return Generic.CompositeMap(F, G)
  end
end

function compose(F::MPolyAnyMap{D, C, <: Map, <: Any}, G::S) where {D, C, S <: Generic.IdentityMap{C}}
  @req codomain(F) === domain(G) "Incompatible (co)domain in composition"
  f = coefficient_map(F)
  if typeof(codomain(f)) === C
    newcoeffmap = compose(f, G)
    return hom(domain(F), codomain(G), newcoeffmap, G.(_images(F)), check=false)
  else
    return Generic.CompositeMap(F, G)
  end
end

################################################################################
#
#  Types computers
#
################################################################################

function morphism_type(::Type{D}, ::Type{C}) where {D <: Union{MPolyRing, MPolyQuoRing}, C <: NCRing}
  return MPolyAnyMap{D, C, Nothing, elem_type(C), HomUngraded}
end

morphism_type(::D, ::C) where {D <: Union{MPolyRing, MPolyQuoRing}, C <: NCRing} =
  morphism_type(D, C)

function morphism_type(::Type{D}, ::Type{C}, f::Type{F}) where {D <: Union{MPolyRing, MPolyQuoRing}, C <: NCRing, F}
  return MPolyAnyMap{D, C, F, elem_type(C), HomUngraded}
end

morphism_type(::D, ::C, ::F) where {D <: Union{MPolyRing, MPolyQuoRing}, C <: NCRing, F} =
  morphism_type(D, C, F)

function (f::MPolyAnyMap{<:MPolyRing, <:AbstractAlgebra.NCRing})(I::MPolyIdeal)
  return ideal(codomain(f), [f(g) for g in gens(I)])
end

images_of_generators(phi::MPolyAnyMap) = _images(phi)

identity_map(R::MPolyRing) = AbstractAlgebra.identity_map(R)

identity_map(R::MPolyQuoRing) = AbstractAlgebra.identity_map(R)

identity_map(Z::ZZRing) = AbstractAlgebra.identity_map(Z)
