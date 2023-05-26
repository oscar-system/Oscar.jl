export CoveredProjectiveScheme
export ProjectiveGlueing
export base_covering
export base_glueing
export base_scheme
export blow_up
export controlled_transform
export covered_scheme
export empty_covered_projective_scheme
export glueing_type
export projective_patches
export strict_transform
export weak_transform

abstract type AbsProjectiveGlueing{
                                   GlueingType<:AbsGlueing,
                                  }
end

### getters for the essential functionality
base_glueing(PG::AbsProjectiveGlueing) = base_glueing(underlying_glueing(PG))
inclusion_maps(PG::AbsProjectiveGlueing) = inclusion_maps(underlying_glueing(PG))
glueing_domains(PG::AbsProjectiveGlueing) = glueing_domains(underlying_glueing(PG))
patches(PG::AbsProjectiveGlueing) = patches(underlying_glueing(PG))
glueing_morphisms(PG::AbsProjectiveGlueing) = glueing_morphisms(underlying_glueing(PG))


@doc raw"""
  LazyProjectiveGlueing(
      X::AbsProjectiveScheme,
      Y::AbsProjectiveScheme,
      BG::AbsGlueing,
      compute_function::Function, 
      glueing_data
    )

Produce a container `pg` to host a non-computed `ProjectiveGlueing` of ``X`` with ``Y``.

The arguments consist of 

 * the patches ``X`` and ``Y`` to be glued;
 * a glueing `BG` of the `base_scheme`s of ``X`` and ``Y`` over which the 
   `ProjectiveGlueing` to be computed sits;
 * a function `compute_function` which takes a single argument `glueing_data` 
   of arbitrary type and actually carries out the computation; 
 * an arbitrary struct `glueing_data` that the user can fill with whatever 
   information is needed to properly feed their `compute_function`. 

The container `pg` can then be stored as the glueing of ``X`` and ``Y``. As soon 
as it is asked about any data on the glueing beyond its two `patches`, it will 
invoke the internally stored `compute_function` to actually carry out the computation 
of the glueing and then serve the incoming request on the basis of that result. 
The latter actual `ProjectiveGlueing` will then be cached. 
"""
mutable struct LazyProjectiveGlueing{
                                     GlueingType<:AbsGlueing,
                                     GlueingDataType
                                    } <: AbsProjectiveGlueing{GlueingType}
  base_glueing::GlueingType
  patches::Tuple{AbsProjectiveScheme, AbsProjectiveScheme}
  compute_function::Function
  glueing_data::GlueingDataType
  underlying_glueing::AbsProjectiveGlueing

  function LazyProjectiveGlueing(
      X::AbsProjectiveScheme,
      Y::AbsProjectiveScheme,
      BG::AbsGlueing,
      compute_function::Function, 
      glueing_data
    )
    (base_scheme(X), base_scheme(Y)) == patches(BG) || error("glueing is incompatible with provided patches")
    return new{typeof(BG), typeof(glueing_data)}(BG, (X, Y), compute_function, glueing_data)
  end
end

# Essential getters
patches(G::LazyProjectiveGlueing) = G.patches
base_glueing(G::LazyProjectiveGlueing) = G.base_glueing

# Everything else will trigger the computation
function underlying_glueing(G::LazyProjectiveGlueing)
  if !isdefined(G, :underlying_glueing)
    G.underlying_glueing = G.compute_function(G.glueing_data)
  end
  return G.underlying_glueing
end

@doc raw"""
    ProjectiveGlueing(
        G::GlueingType, 
        incP::IncType, incQ::IncType,
        f::IsoType, g::IsoType;
        check::Bool=true
      ) where {GlueingType<:AbsGlueing, IncType<:ProjectiveSchemeMor, IsoType<:ProjectiveSchemeMor}

The `AbsProjectiveSchemeMorphism`s `incP` and `incQ` are open embeddings over open 
embeddings of their respective `base_scheme`s. 

        PX ↩ PU ≅ QV ↪ QY
      π ↓    ↓    ↓    ↓ π
    G : X  ↩ U  ≅ V  ↪ Y 

This creates a glueing of the projective schemes `codomain(incP)` and `codomain(incQ)` 
over a glueing `G` of their `base_scheme`s along the morphisms of `AbsProjectiveScheme`s 
`f` and `g`, identifying `domain(incP)` and `domain(incQ)`, respectively.
"""
mutable struct ProjectiveGlueing{
                                 GlueingType<:AbsGlueing,
                                 IsoType<:ProjectiveSchemeMor,
                                 IncType<:ProjectiveSchemeMor
                                } <: AbsProjectiveGlueing{GlueingType}
  G::GlueingType # the underlying glueing of the base schemes
  inc_to_P::IncType
  inc_to_Q::IncType
  f::IsoType
  g::IsoType

  ### 
  # Given two relative projective schemes and a glueing 
  #
  #       PX ↩ PU ≅ QV ↪ QY
  #     π ↓    ↓    ↓    ↓ π
  #   G : X  ↩ U  ≅ V  ↪ Y 
  #
  # this constructs the glueing of PX and QY along 
  # their open subsets PU and QV, given the two inclusions 
  # and isomorphisms over the glueing G in the base schemes.
  function ProjectiveGlueing(
      G::GlueingType, 
      incP::IncType, incQ::IncType,
      f::IsoType, g::IsoType;
      check::Bool=true
    ) where {GlueingType<:AbsGlueing, IncType<:ProjectiveSchemeMor, IsoType<:ProjectiveSchemeMor}
    (X, Y) = patches(G)
    (U, V) = glueing_domains(G)
    (fb, gb) = glueing_morphisms(G)
    (PX, QY) = (codomain(incP), codomain(incQ))
    (PU, QV) = (domain(incP), domain(incQ))
    (base_scheme(PX) == X && base_scheme(QY) == Y) || error("base glueing is incompatible with the projective schemes")
    domain(f) == codomain(g) == PU && domain(g) == codomain(f) == QV || error("maps are not compatible")
    SPU = homogeneous_coordinate_ring(domain(f))
    SQV = homogeneous_coordinate_ring(codomain(f))
    if check
      # check the commutativity of the pullbacks
      all(y->(pullback(f)(SQV(OO(V)(y))) == SPU(pullback(fb)(OO(V)(y)))), gens(base_ring(OO(Y)))) || error("maps do not commute")
      all(x->(pullback(g)(SPU(OO(U)(x))) == SQV(pullback(gb)(OO(U)(x)))), gens(base_ring(OO(X)))) || error("maps do not commute")
      fc = map_on_affine_cones(f, check=false)
      gc = map_on_affine_cones(g, check=false)
      idCPU = compose(fc, gc)
      idCPU == identity_map(domain(fc)) || error("composition of maps is not the identity")
      idCQV = compose(gc, fc)
      idCQV == identity_map(domain(gc)) || error("composition of maps is not the identity")
      # idPU = compose(f, g)
      # all(t->(pullback(idPU)(t) == t), gens(SPU)) || error("composition of maps is not the identity")
      # idQV = compose(g, f)
      # all(t->(pullback(idQV)(t) == t), gens(SQV)) || error("composition of maps is not the identity")
    end
    return new{GlueingType, IsoType, IncType}(G, incP, incQ, f, g)
  end
end

### type getters

glueing_type(P::T) where {T<:ProjectiveScheme} = ProjectiveGlueing{glueing_type(base_scheme_type(T)), T, morphism_type(T)}
glueing_type(::Type{T}) where {T<:ProjectiveScheme} = ProjectiveGlueing{glueing_type(base_scheme_type(T)), T, morphism_type(T)}

### essential getters

base_glueing(PG::ProjectiveGlueing) = PG.G
inclusion_maps(PG::ProjectiveGlueing) = (PG.inc_to_P, PG.inc_to_Q)
glueing_domains(PG::ProjectiveGlueing) = (domain(PG.f), domain(PG.g))
patches(PG::ProjectiveGlueing) = (codomain(PG.inc_to_P), codomain(PG.inc_to_Q))
glueing_morphisms(PG::ProjectiveGlueing) = (PG.f, PG.g)

### Proper schemes π : Z → X over a covered base scheme X
# 
# When {Uᵢ} is an affine covering of X, the datum stored 
# consists of a list of projective schemes 
#
#   Zᵢ ⊂ ℙʳ⁽ⁱ⁾(𝒪(Uᵢ)) → Uᵢ
#
# with varying ambient spaces ℙʳ⁽ⁱ⁾(𝒪(Uᵢ)) and a list of 
# identifications (transitions) 
#
#   Zᵢ ∩ π⁻¹(Uⱼ) ≅ Zⱼ ∩ π⁻¹(Uᵢ)
#
# of projective schemes over Uᵢ∩ Uⱼ for all pairs (i,j).
#
# These structs are designed to accommodate blowups of 
# covered schemes along arbitrary centers, as well as 
# projective bundles. 

@attributes mutable struct CoveredProjectiveScheme{BRT} <: Scheme{BRT}
  Y::AbsCoveredScheme # the base scheme
  BC::Covering # the reference covering of the base scheme
  patches::IdDict{AbsSpec, AbsProjectiveScheme} # the projective spaces over the affine patches in the base covering
  glueings::IdDict{Tuple{AbsSpec, AbsSpec}, AbsProjectiveGlueing} # the transitions sitting over the affine patches in the glueing domains of the base scheme

  function CoveredProjectiveScheme(
      Y::AbsCoveredScheme,
      C::Covering,
      projective_patches::IdDict{AbsSpec, AbsProjectiveScheme},
      projective_glueings::IdDict{Tuple{AbsSpec, AbsSpec}, AbsProjectiveGlueing};
      check::Bool=true
    )
    C in coverings(Y) || error("covering not listed")
    for P in values(projective_patches)
      base_scheme(P) in patches(C) || error("base scheme not found in covering")
    end
    for (U, V) in keys(glueings(C))
      (U, V) in keys(projective_glueings) || error("not all projective glueings were provided")
    end
    return new{base_ring_type(Y)}(Y, C, projective_patches, projective_glueings)
  end
end

base_scheme(P::CoveredProjectiveScheme) = P.Y
base_covering(P::CoveredProjectiveScheme) = P.BC
base_patches(P::CoveredProjectiveScheme) = patches(P.BC)
projective_patches(P::CoveredProjectiveScheme) = values(P.patches)
getindex(P::CoveredProjectiveScheme, U::AbsSpec) = (P.patches)[U]
getindex(P::CoveredProjectiveScheme, U::AbsSpec, V::AbsSpec) = (P.glueings)[(U, V)]

function empty_covered_projective_scheme(R::T) where {T<:AbstractAlgebra.Ring}
  Y = empty_covered_scheme(R)
  C = default_covering(Y)
  #U = C[1]
  #ST = affine_patch_type(Y)
  pp = Dict{AbsSpec, AbsProjectiveScheme}()
  #P = projective_space(U, 0)
  #pp[U] = P
  tr = Dict{Tuple{AbsSpec, AbsSpec}, AbsProjectiveGlueing}()
  #W = SpecOpen(U)
  #PW, inc = fiber_product(restriction_map(U, W), P)
  #tr[(U, U)] = ProjectiveGlueing(Glueing(U, U, identity_map(W), identity_map(W)), 
                                 #inc, inc, identity_map(PW), identity_map(PW))
  return CoveredProjectiveScheme(Y, C, pp, tr)
end

@doc raw"""
    blow_up_chart(W::AbsSpec, I::Ideal)

Return the blowup of ``W`` at the ideal ``I``; this is a `ProjectiveScheme`
with `base_scheme` ``W``.

!!! note
    blow_up relies on this internal method for computing the blow ups of all chartsand appropriately assembles the returned projective schemes to a single coverec scheme.
"""

function blow_up_chart(W::AbsSpec, I::Ideal; var_name::VarName = :s)
  error("method `blow_up_chart` not implemented for arguments of type $(typeof(W)) and $(typeof(I))")
end

########################################################################
# Blowups of affine schemes                                            #
########################################################################

function blow_up_chart(W::AbsSpec{<:Field, <:MPolyRing}, I::MPolyIdeal;
    var_name::VarName = :s
  )
  base_ring(I) === OO(W) || error("ideal does not belong to the correct ring")
#  if one(OO(W)) in I 
#    error("blowing up along the unit ideal; this case should be caught earlier")
#    # construct a relative ℙ⁰ over W with its identifying projection?
#  end
  r = ngens(I) - 1
  g = gens(I)
  IPW = projective_space(W, r, var_name=var_name)
  S = homogeneous_coordinate_ring(IPW)
  t = gens(S)
  if is_regular_sequence(gens(I))
    # construct the blowup manually
    J = ideal(S, [t[i]*g[j] - t[j]*g[i] for i in 1:r, j in i+1:r+1])
    IPY = subscheme(IPW, J)
    # Compute the IdealSheaf for the exceptional divisor
    ID = IdDict{AbsSpec, RingElem}()
    Y = covered_scheme(IPY)
    p = covered_projection_to_base(IPY)
    p_cov = covering_morphism(p)
    for i in 1:ngens(I)
      U = affine_charts(Y)[i]
      p_res = p_cov[U]
      W === codomain(p_res) || error("codomain not correct")
      ID[affine_charts(Y)[i]] = pullback(p_res)(gen(I, i))
    end
    E = oscar.EffectiveCartierDivisor(Y, ID, trivializing_covering=domain(p_cov), check=false)
    set_attribute!(Y, :exceptional_divisor, E)
    set_attribute!(IPY, :exceptional_divisor, E)
    return IPY
  else
    # construct the blowup by elimination.
    CIPW, pullback_to_cone = affine_cone(IPW)
    R = base_ring(I)
    x = gens(R)
    kk = coefficient_ring(R)
    A, x_ext = polynomial_ring(kk, vcat(symbols(R), [:t]))
    t = last(x_ext) 
    inc = hom(R, A, x_ext[1:end-1])
    phi = hom(OO(CIPW), A, vcat([inc(g[i])*t for i in 1:r+1], x_ext[1:end-1], )) # the homogeneous variables come first
    J = kernel(phi)
    pb = inverse(pullback_to_cone)
    Jh = ideal(homogeneous_coordinate_ring(IPW), pb.(lifted_numerator.(gens(J))))
    IPY = subscheme(IPW, Jh)
    # Compute the IdealSheaf for the exceptional divisor
    ID = IdDict{AbsSpec, RingElem}()
    Y = covered_scheme(IPY)
    p = covered_projection_to_base(IPY)
    p_cov = covering_morphism(p)
    for i in 1:ngens(I)
      U = affine_charts(Y)[i]
      p_res = p_cov[U]
      W === codomain(p_res) || error("codomain not correct")
      ID[affine_charts(Y)[i]] = pullback(p_res)(gen(I, i))
    end
    E = oscar.EffectiveCartierDivisor(Y, ID, trivializing_covering=domain(p_cov), check=false) 
    set_attribute!(Y, :exceptional_divisor, E)
    set_attribute!(IPY, :exceptional_divisor, E)
    return IPY
  end
end

function blow_up_chart(W::AbsSpec{<:Field, <:RingType}, I::Ideal;
    var_name::VarName = :s
  ) where {RingType<:Union{MPolyQuoRing, MPolyLocRing, MPolyQuoLocRing}}
  base_ring(I) === OO(W) || error("ideal does not belong to the correct ring")

  # It follows the generic Proj construction
  R = OO(W)
  T, (t,) = polynomial_ring(R, ["t"])
  S, s = grade(polynomial_ring(R, [Symbol(var_name, i-1) for i in 1:ngens(I)])[1])
  phi = hom(S, T, [t*g for g in gens(I)])
  K = kernel(phi)
  K = ideal(S, [g for g in gens(K) if !iszero(g)]) # clean up superfluous generators
  Bl_W = ProjectiveScheme(S, K)
  set_base_scheme!(Bl_W, W)
  # Compute the IdealSheaf for the exceptional divisor
  ID = IdDict{AbsSpec, RingElem}()
  Y = covered_scheme(Bl_W)
  p = covered_projection_to_base(Bl_W)
  p_cov = covering_morphism(p)
  for i in 1:ngens(I)
    U = affine_charts(Y)[i]
    p_res = p_cov[U]
    W === codomain(p_res) || error("codomain not correct")
    ID[affine_charts(Y)[i]] = pullback(p_res)(gen(I, i))
  end
  E = oscar.EffectiveCartierDivisor(Y, ID, trivializing_covering=domain(p_cov), check=false)
  set_attribute!(Y, :exceptional_divisor, E)
  set_attribute!(Bl_W, :exceptional_divisor, E)
  return Bl_W
end

function is_regular_sequence(g::Vector{T}) where {T<:RingElem}
  length(g) == 0 && return true
  R = parent(g[1])
  all(x->parent(x)===R, g) || error("elements do not belong to the correct ring")
  isunit(g[1]) && return false # See Bruns-Herzog: Cohen-Macaulay rings, section 1.1.
  is_zero_divisor(g[1]) && return false
  A, p = quo(R, ideal(R, g))
  return is_regular_sequence(p.(g[2:end]))
end


#function blow_up(W::Spec, I::MPolyQuoLocalizedIdeal;
#    var_names::Vector{Symbol}=Symbol.(["s$(i-1)" for i in 1:ngens(I)]),
#    verbose::Bool=false,
#    check::Bool=true,
#    is_regular_sequence::Bool=false
#  )
#  base_ring(I) == OO(W) || error("ideal does not belong to the correct ring")
#  r = ngens(I)-1
#  PW = projective_space(W, var_names)
#  PWC = affine_cone(PW)
#  prW = projection_to_base(PW)
#  WA1, pW, pA = product(W, affine_space(base_ring(base_ring(OO(W))), 1, variable_name="t"))
#  t = pullback(pA)(OO(codomain(pA))(base_ring(OO(codomain(pA)))[1]))
#  imgs = vcat((x->t*x).(pullback(pW).(gens(I))), pullback(pW).(gens(base_ring(OO(W)))))
#  inner_phi = hom(base_ring(OO(PWC)), OO(WA1), imgs)
#  phi = MPolyQuoLocalizedRingHom(OO(PWC), OO(WA1), inner_phi)
#  K = kernel(phi)
#  J = ideal(homog_poly_ring(PW), poly_to_homog(PW).(lifted_numerator.([f for f in gens(K) if !iszero(f)])))
#  return subscheme(PW, J)
#end

## blow up X in the center described by g using these explicit generators.
#function blow_up(
#    W::Spec, 
#    I::Vector{RingElemType}; 
#    var_names::Vector{Symbol}=Vector{Symbol}(), 
#    verbose::Bool=false,
#    check::Bool=true,
#    is_regular_sequence::Bool=false
#  ) where {RingElemType<:MPolyRingElem}
#
#  # some internal function
#  function _add_variables(R::RingType, v::Vector{Symbol}) where {RingType<:MPolyRing}
#    ext_R, _ = polynomial_ring(coefficient_ring(R), vcat(symbols(R), v))
#    n = ngens(R)
#    phi = AlgebraHomomorphism(R, ext_R, gens(ext_R)[1:n])
#    return ext_R, phi, gens(ext_R)[(ngens(R)+1):ngens(ext_R)]
#  end
#
#  A = OO(W)
#  R = base_ring(A)
#  #TODO: Check if for all i \in I parent(i) == R
#
#  m = length(I)
#  Pw = (length(var_names) > 0 ? projective_space(W, var_names) : projective_space(W,m-1))
#  S = homog_poly_ring(Pw)
#
#  CP = affine_cone(Pw)
#  Polyring = base_ring(OO(CP))
#  if !is_regular_sequence
#    At, embeddingAt, T =  _add_variables(R,[:t])
#    t = T[1]
#
#    #	@show vcat([t*embeddingAt(f) for f in I], gens(At)[1:end-1])
#    Phi = AlgebraHomomorphism(Polyring, At, vcat([t*embeddingAt(f) for f in I], gens(At)[1:end-1]))
#
#    Imod = modulus(A)
#    IW = ideal(At, [embeddingAt(f) for f in gens(Imod)])
#    @show "start groebner basis computation."
#    IWpre = preimage(Phi, IW)
#    @show "done. Proceed to process"
#    #SIWpre = ideal(S,[frac_to_homog(Pw,g) for g in gens(IWpre)])
#    SIWpre = ideal(S, poly_to_homog(Pw).(gens(IWpre)))
#    @show 0
#
#    projective_version = subscheme(Pw, SIWpre)
#    @show 1
#    covered_version = as_covered_scheme(projective_version)
#    @show 2
#    projection_map = get_attribute(projective_version, :covered_projection_to_base)
#    @show 3
#    E_dict = Dict{affine_patch_type(covered_version), Vector{RingElemType}}()
#    for i in 1:length(I)
#      @show i
#      U = default_covering(covered_version)[i]
#      E_dict[U] = [lifted_numerator(pullback(projection_map[U])(I[i]))]
#    end
#    @show 4
#    exc_div = IdealSheaf(covered_version, default_covering(covered_version), E_dict, check=false)
#    @show "processing done."
#    return projective_version, covered_version, projection_map, exc_div
#  else
#    A = zero_matrix(S, 2, ngens(S))
#    for i in 1:ngens(S)
#      A[1, i] = S[i]
#      A[2, i] = I[i]
#    end
#    IW = ideal(S, minors(A, 2))
#    projective_version = subscheme(Pw, IW)
#    covered_ambient = as_covered_scheme(Pw)
#    ambient_projection_map = get_attribute(Pw, :covered_projection_to_base)
#    C = default_covering(covered_ambient)
#    SpecType = affine_patch_type(covered_ambient)
#    PolyType = poly_type(SpecType)
#    Idict = Dict{SpecType, ideal_type(ring_type(SpecType))}()
#    @show I
#    for i in 1:npatches(C)
#      @show i
#      v = gens(OO(C[i]))
#      phi = pullback(ambient_projection_map[C[i]])
#      loc_eqns = vcat([v[j]*phi(I[i])-phi(I[j]) for j in 1:i-1], [v[j]*phi(I[i])-phi(I[j+1]) for j in i:length(I)-1])
#      @show loc_eqns
#      Idict[C[i]] = ideal(OO(C[i]), loc_eqns)
#      @show "ideal saved"
#    end
#    @show "computing ideal sheaf"
#    Itrans = IdealSheaf(covered_ambient, C, Idict, check=false)
#    @show "computation done; getting subscheme"
#    covered_version = subscheme(Itrans)
#    @show "computation done"
#    set_attribute!(projective_version, :as_covered_scheme, covered_version)
#    @show "computing a default covering"
#    set_attribute!(projective_version, :standard_covering, default_covering(covered_version))
#    @show "done"
#
#    proj_dict = Dict{SpecType, morphism_type(SpecType, SpecType)}()
#    @show "doing some other stuff"
#    for i in 1:length(I)
#      @show i
#      Z = patches(default_covering(covered_version))[i]
#      U = patches(default_covering(covered_ambient))[i]
#      proj_dict[Z] = restrict(ambient_projection_map[U], Z, codomain(ambient_projection_map[U]), check=false)
#    end
#
#    projection_map = CoveringMorphism(default_covering(covered_version), Covering(W), proj_dict)
#    set_attribute!(projective_version, :covered_projection_to_base, projection_map)
#    @show 3
#    E_dict = Dict{affine_patch_type(covered_version), ideal_type(ring_type(affine_patch_type(covered_version)))}()
#    for i in 1:length(I)
#      @show i
#      U = default_covering(covered_version)[i]
#      E_dict[U] = ideal(OO(U), pullback(projection_map[U])(I[i]))
#    end
#    @show 4
#    exc_div = IdealSheaf(covered_version, default_covering(covered_version), E_dict, check=false)
#    @show "processing done."
#    return projective_version, covered_version, projection_map, exc_div
#  end
#end

# This is a sample for how to use LazyProjectiveGlueings. 
# Originally, we probably had some body of a double for-loop iterating over 
# pairs of patches (P, Q) that we need to glue. We take that body which 
# computes the glueing of P and Q and move 
# it to an external function (here _compute_projective_glueing). Then we 
# go through all the local variables in the body of the for-loop which 
# are needed for the actual computation and create a tailor-made struct to 
# house them. We add an extraction section in the beginning of the compute 
# function to restore them and recreate the original setting within the 
# for-loops. 
#
# Finally, we can replace the inner part of the double for-loop with the 
# actual wrap-up of the local variables and feed everything to a constructor 
# for a `LazyProjectiveGlueing` as documented above. 
struct CoveredProjectiveGlueingData
  U::AbsSpec
  V::AbsSpec
  P::AbsProjectiveScheme
  Q::AbsProjectiveScheme
  G::AbsGlueing
  I::IdealSheaf
end

function _compute_projective_glueing(gd::CoveredProjectiveGlueingData)
  P = gd.P
  Q = gd.Q
  U = gd.U
  V = gd.V
  G = gd.G
  I = gd.I
  X = scheme(I)
  OX = StructureSheafOfRings(X)

  SP = homogeneous_coordinate_ring(P)
  SQ = homogeneous_coordinate_ring(Q)
  UV, VU = glueing_domains(G)
  f, g = glueing_morphisms(G)

  PUV, PUVtoP = fiber_product(OX(U, UV), P)
  QVU, QVUtoQ = fiber_product(OX(V, VU), Q)

  # to construct the identifications of PUV with QVU we need to 
  # express the generators of I(U) in terms of the generators of I(V)
  # on the overlap U ∩ V. 
  !(G isa Glueing) || error("method not implemented for this type of glueing")

  # The problem is that on a SpecOpen U ∩ V
  # despite I(U)|U ∩ V == I(V)|U ∩ V, we 
  # have no method to find coefficients aᵢⱼ such that fᵢ = ∑ⱼaᵢⱼ⋅gⱼ
  # for the generators fᵢ of I(U) and gⱼ of I(V): Even though 
  # we can do this locally on the patches of a SpecOpen, the result 
  # is not guaranteed to glue to global functions on the overlap.
  # Abstractly, we know that the intersection of affine charts 
  # in a separated scheme must be affine, but we do not have a 
  # model of this overlap as an affine scheme and hence no computational
  # backup. 

  # fᵢ the generators of I(U)
  # gⱼ the generators of I(V)
  # aᵢⱼ the coefficients for fᵢ = ∑ⱼ aᵢⱼ⋅gⱼ in VU
  # bⱼᵢ the coefficients for gⱼ = ∑ᵢ bⱼᵢ⋅fᵢ in UV
  # sᵢ the variables for the homogeneous ring over U
  # tⱼ the variables for the homogenesous ring over V
  A = [coordinates(OX(U, VU)(f), I(VU)) for f in gens(I(U))] # A[i][j] = aᵢⱼ
  B = [coordinates(OX(V, UV)(g), I(UV)) for g in gens(I(V))] # B[j][i] = bⱼᵢ
  SQVU = homogeneous_coordinate_ring(QVU)
  SPUV = homogeneous_coordinate_ring(PUV)
  # the induced map is ℙ(UV) → ℙ(VU), tⱼ ↦ ∑ᵢ bⱼᵢ ⋅ sᵢ 
  # and ℙ(VU) → ℙ(UV), sᵢ ↦ ∑ⱼ aᵢⱼ ⋅ tⱼ 
  fup = ProjectiveSchemeMor(PUV, QVU, hom(SQVU, SPUV, pullback(f), [sum([B[j][i]*SPUV[i] for i in 1:ngens(SPUV)]) for j in 1:length(B)], check=false), check=false)
  gup = ProjectiveSchemeMor(QVU, PUV, hom(SPUV, SQVU, pullback(g), [sum([A[i][j]*SQVU[j] for j in 1:ngens(SQVU)]) for i in 1:length(A)], check=false), check=false)
  return ProjectiveGlueing(G, PUVtoP, QVUtoQ, fup, gup, check=false)
end


@doc raw"""
    blow_up(X::AbsSpec, I::Ideal)

Return the blow-up morphism of blowing up ``X`` at ``I`` in ``OO(X)``.
"""
function blow_up(
    X::AbsSpec{<:Any, <:MPolyAnyRing},
    I::MPolyAnyIdeal)

  R = OO(X)
  @req R == base_ring(I) "I must be an ideal in the coordinate ring of X"
  ## prepare trivially covered scheme and ideal sheaf on it as input for blow_up(I::IdealSheaf)
  C = Covering([X])
  XX = CoveredScheme(C)
  Isheaf = ideal_sheaf(XX,X,gens(I))
  return blow_up(Isheaf)
end

@doc raw"""
    blow_up(I::IdealSheaf)

Return the blow-up morphism of blowing up of the underlying scheme of ``I``  at ``I``.
"""
function blow_up(
    I::IdealSheaf;
    verbose::Bool=false,
    check::Bool=true,
    var_name::VarName=:s,
    covering::Covering=default_covering(scheme(I))
  )
  X = space(I)
  local_blowups = IdDict{AbsSpec, AbsProjectiveScheme}()
  for U in patches(covering)
    local_blowups[U] = blow_up_chart(U, I(U), var_name=var_name)
  end
  projective_glueings = IdDict{Tuple{AbsSpec, AbsSpec}, AbsProjectiveGlueing}()

  # prepare for the projective glueings
  for (U, V) in keys(glueings(covering))
    P = local_blowups[U]
    Q = local_blowups[V]
    G = covering[U, V]
    gd = CoveredProjectiveGlueingData(U, V, P, Q, G, I)
    projective_glueings[U, V] = LazyProjectiveGlueing(P, Q, G, _compute_projective_glueing, gd)
  end
  Bl_I = CoveredProjectiveScheme(X, covering, local_blowups, projective_glueings)
  Y = covered_scheme(Bl_I)
  # Assemble the exceptional divisor as a Cartier divisor
  ID = IdDict{AbsSpec, RingElem}()
  p = covered_projection_to_base(Bl_I)
  p_cov = covering_morphism(p)
  for U in patches(domain(p_cov)) # These coincide with the patches on which the local E's are defined
    V = codomain(p_cov[U])
    if !has_attribute(local_blowups[V], :exceptional_divisor)
      error("exceptional divisor was not cached on local blowup.")
    end
    E_loc = get_attribute(local_blowups[V], :exceptional_divisor)::EffectiveCartierDivisor
    g = E_loc(U)
    isone(length(g)) || error("exceptional divisor must be given by one equation in these charts.")
    ID[U] = first(E_loc(U))
  end
  pr = BlowupMorphism(Bl_I, I)
  pr.exceptional_divisor = EffectiveCartierDivisor(Y, ID, trivializing_covering=domain(p_cov), check=false)
  return pr
end

# The following is a wrap-up of local variables necessary to compute glueings
# from morphisms of graded algebras. It will be used to fill in the 
# `LazyGlueing` together with the compute function following below.
struct ProjectiveGlueingData
  down_left::AbsSpec
  down_right::AbsSpec
  up_left::AbsSpec
  up_right::AbsSpec
  down_covering::Covering
  projective_scheme::CoveredProjectiveScheme
end

# This function actually computes the glueing with the data extracted 
# from the ProjectiveGlueingData. 
#
# Originally, this was in the body of a constructor. But we want to 
# postpone the actual computation, so we wrap up all necessary local 
# variables in a `ProjectiveGlueingData` and copy-paste the required 
# part of code into the body of this function. Together with a header 
# which is restoring our local variables. 
function _compute_glueing(gd::ProjectiveGlueingData)
  U = gd.down_left
  V = gd.down_right
  UW = gd.up_left
  VW = gd.up_right
  C = gd.down_covering
  # Now we have the following diagram 
  #
  #               fup, gup
  # P[U] ⊃ UW ⇢  UD ↔ VD ⇠ VW ⊂ P[V]
  #         ↓    ↓    ↓    ↓
  #         U ⊃  A  ↔  B ⊂ V
  #                f,g
  #
  # with UW = {sᵢ≠ 0} and VW = {tⱼ≠ 0}. 
  P = gd.projective_scheme
  (A, B) = glueing_domains(C[U, V])
  (f, g) = glueing_morphisms(C[U, V])
  (UD, VD) = glueing_domains(P[U, V])
  (fup, gup) = glueing_morphisms(P[U, V])
  (incU, incV) = inclusion_maps(P[U, V])
  S = homogeneous_coordinate_ring(P[U])
  T = homogeneous_coordinate_ring(P[V])
  i = P[U][UW][2]
  j = P[V][VW][2]
  s_i = gen(S, i)
  t_j = gen(T, j)
  AW = affine_charts(Oscar.covered_scheme(UD))[i]
  BW = affine_charts(Oscar.covered_scheme(VD))[j]

  hU = dehomogenization_map(UD, AW)(pullback(fup)(pullback(incV)(t_j)))
  hV = dehomogenization_map(VD, BW)(pullback(gup)(pullback(incU)(s_i)))

  # We need to construct the glueing 
  #
  #          f', g'
  #   UW ↩ AAW ↔ BBW ↪ VW
  #   
  # as follows. 
  # From the base glueing we already have that UW differs 
  # from AAW by the pullback of the equation `complement_equation(A)`.
  # But then, also `hU` cuts out another locus which needs to be 
  # removed before glueing. 
  #
  # The same applies to the other side. Finally, we need to track the 
  # coordinates through the homogenization-glueing-pullback-dehomogenization 
  # machinery to provide the glueing morphisms. 

  ptbUD = covered_projection_to_base(UD)
  ptbVD = covered_projection_to_base(VD)
  hhU = lifted_numerator(pullback(ptbUD[AW])(complement_equation(A)))
  hhU = hhU * lifted_numerator(hU)
  AAW = PrincipalOpenSubset(UW, OO(UW)(hhU))
  hhV = lifted_numerator(pullback(ptbVD[BW])(complement_equation(B)))
  hhV = hhV * lifted_numerator(hV)
  BBW = PrincipalOpenSubset(VW, OO(VW)(hhV))

  x = gens(ambient_coordinate_ring(AAW))
  y = gens(ambient_coordinate_ring(BBW))

  xh = homogenization_map(UD, AW).(OO(AW).(x))
  yh = homogenization_map(VD, BW).(OO(BW).(y))

  xhh = [(pullback(gup)(pp), pullback(gup)(qq)) for (pp, qq) in xh]
  yhh = [(pullback(fup)(pp), pullback(fup)(qq)) for (pp, qq) in yh]

  phi = dehomogenization_map(VD, BW)
  psi = dehomogenization_map(UD, AW) 
  yimgs = [OO(AAW)(psi(pp))*inv(OO(AAW)(psi(qq))) for (pp, qq) in yhh]
  ximgs = [OO(BBW)(phi(pp))*inv(OO(BBW)(phi(qq))) for (pp, qq) in xhh]
  ff = SpecMor(AAW, BBW, hom(OO(BBW), OO(AAW), yimgs, check=false), check=false)
  gg = SpecMor(BBW, AAW, hom(OO(AAW), OO(BBW), ximgs, check=false), check=false)

  return SimpleGlueing(UW, VW, ff, gg, check=false)
  pr.exceptional_divisor = EffectiveCartierDivisor(Y, ID, trivializing_covering=domain(p_cov))
  return pr
end

@attr function covered_scheme(P::CoveredProjectiveScheme)
  X = base_scheme(P)
  C = base_covering(P)
  new_patches = Vector{AbsSpec}()
  new_glueings = IdDict{Tuple{AbsSpec, AbsSpec}, AbsGlueing}()
  projection_dict = IdDict{AbsSpec, AbsSpecMor}()
  parts = IdDict{AbsSpec, AbsCoveredScheme}() 
  for U in patches(C)
    parts[U] = Oscar.covered_scheme(P[U])
  end
  # We first assemble a Covering where the preimage of every base 
  # patch appears as one connected component
  result_patches = vcat([affine_charts(parts[U]) for U in patches(C)]...)
  result_glueings = IdDict{Tuple{AbsSpec, AbsSpec}, AbsGlueing}()
  for U in patches(C)
    GG = glueings(parts[U])
    for (A, B) in keys(GG)
      result_glueings[(A, B)] = GG[(A, B)]
    end
  end
  result_covering = Covering(result_patches, result_glueings, check=false)

  # Now we need to add glueings
  for U in patches(C)
    for V in patches(C)
      U === V && continue
      for UW in affine_charts(parts[U])
        for VW in affine_charts(parts[V])
          GGG = LazyGlueing(UW, VW, _compute_glueing, 
                            ProjectiveGlueingData(U, V, UW, VW, C, P)
                           )
          result_glueings[(UW, VW)] = GGG
        end
      end
    end
  end
  result_covering = Covering(result_patches, result_glueings, check=false)
  result = CoveredScheme(result_covering)

  projection_dict = IdDict{AbsSpec, AbsSpecMor}()
  for U in patches(C)
    PP = P[U]
    p = covered_projection_to_base(PP)
    cov_mor = covering_morphism(p)
    cov_mor_dict = morphisms(cov_mor)
    for V in keys(cov_mor_dict)
      projection_dict[V] = cov_mor_dict[V]
    end
  end

  # TODO: Remove the internal checks in the constructors below
  covering_map = CoveringMorphism(result_covering, C, projection_dict) 
  set_attribute!(P, :covering_projection_to_base, covering_map)
  return result
end

@attr function covered_projection_to_base(P::CoveredProjectiveScheme) 
  if !has_attribute(P, :covering_projection_to_base)
    covered_scheme(P)
  end
  covering_pr = get_attribute(P, :covering_projection_to_base)::CoveringMorphism
  return CoveredSchemeMorphism(covered_scheme(P), base_scheme(P), covering_pr)
end


# function strict_transform(f::SpecMor, h::Vector{PolyType}, g::Vector{PolyType}) where{PolyType<:MPolyRingElem}
	#(IOw: Exc Div ^\infty)
	
#        X = domain(f)
#        Y = codomain(f)
#        R = base_ring(OO(X))
#	Excdiv = ideal(h)
#
#	Pf = pullback(f)
#        Iold = ideal(R, lifted_numerator.(Pf.(g))) + strict_modulus(X)
#	
#	while true
#          Inew = quotient(Iold, Excdiv)
#          Iold == Inew && break
#          Iold = Inew
#	end
#        return gens(Iold)
#end
#
#
#
#function total_transform(f::SpecMor, h::Vector{PolyType}, g::Vector{PolyType}) where{PolyType<:MPolyRingElem}
#	#IOw
#	
#        X = domain(f)
#        Y = codomain(f)
#        R = base_ring(OO(X))
#	Excdiv = ideal(h)
#
#	Pf = pullback(f)
#        Iold = ideal(R, lifted_numerator.(Pf.(g))) + strict_modulus(X)
#        return gens(Iold)
#end
#
#
#### NOT TESTED YET
#function weak_transform(f::SpecMor, h::Vector{PolyType}, g::Vector{PolyType}) where{PolyType<:MPolyRingElem}
#
#	X = domain(f)
#        Y = codomain(f)
#        R = base_ring(OO(X))
#	Excdiv = ideal(h)
#
#	Pf = pullback(f)
#        Iold = ideal(R, lifted_numerator.(Pf.(g))) + strict_modulus(X)
#
#	while true
#		Inew = quotient(Iold, Excdiv)
#		!(Iold == Excdiv * Inew) && break
#		Iold = Inew 
#	end
#	return gens(Iold)
#	#(IOw : Exc Div ^k), k maximal
#end
#
#### NOT TESTED YET
#function controlled_transform(f::SpecMor, h::Vector{PolyType}, g::Vector{PolyType}, i::Int) where{PolyType<:MPolyRingElem}
#	#(IOw : Exc Div ^i)
#
#	X = domain(f)
#        Y = codomain(f)
#        R = base_ring(OO(X))
#	Excdiv = ideal(h)
#
#	Pf = pullback(f)
#        Iold = ideal(R, lifted_numerator.(Pf.(g))) + strict_modulus(X)
#	
#
#	for j in 1:i
#		Inew = quotient(Iold,Excdiv)
#		Inew == 1 && break
#		Iold = Inew
#	end
#	
#	return gens(Iold)
#
#end
#
#function strict_transform(
#    P::AbsProjectiveScheme, 
#    E::IdealSheaf, 
#    I::Vector{PolyType}
#  ) where {PolyType<:MPolyRingElem}
#  X = as_covered_scheme(P)
#  C = covering(E) 
#  C in coverings(X) || error("covering not found")
#  Y = base_scheme(P)
#  length(I) > 0 || return IdealSheaf(X) # return the zero ideal sheaf
#  for i in 1:length(I)
#    parent(I[i]) == base_ring(OO(Y)) || error("polynomials do not belong to the correct ring")
#  end
#  f = covered_projection_to_base(P)
#  if domain(f) !== C 
#    f = compose(X[C, domain(f)], f)
#  end
#
#  SpecType = affine_patch_type(X)
#  trans_dict = Dict{SpecType, Vector{poly_type(SpecType)}}()
#  for U in patches(C)
#    println("computing the strict transform on the patch $U")
#    trans_dict[U] = strict_transform(f[U], E[U], I)
#  end
#  return IdealSheaf(X, C, trans_dict, check=false)
#  #return IdealSheaf(X, CX, trans_dict, check=false)
#end
#
#
#function strict_transform(f::CoveredSchemeMorphism, E::IdealSheaf, I::IdealSheaf)
#  X = domain(f)
#  Y = codomain(f)
#  CX = domain_covering(f)
#  CY = codomain_covering(f)
#  SpecType = affine_patch_type(X)
#  trans_dict = Dict{SpecType, Vector{poly_type(SpecType)}}()
#  for U in patches(CX)
#    trans_dict[U] = strict_transform(f[U], E[U], I[codomain(f[U])])
#  end
#  return IdealSheaf(X, CX, trans_dict, check=true)
#  #return IdealSheaf(X, CX, trans_dict, check=false)
#end
#
#function total_transform(f::CoveredSchemeMorphism, E::IdealSheaf, I::IdealSheaf)
#  X = domain(f)
#  Y = codomain(f)
#  CX = domain_covering(f)
#  CY = codomain_covering(f)
#  SpecType = affine_patch_type(X)
#  trans_dict = Dict{SpecType, Vector{poly_type(SpecType)}}()
#  for U in patches(CX)
#    trans_dict[U] = total_transform(f[U], E[U], I[codomain(f[U])])
#  end
#  return IdealSheaf(X, CX, trans_dict, check=true)
#end
#
#function prepare_smooth_center(I::IdealSheaf; check::Bool=true)
#  X = scheme(I)
#  C = covering(I)
#  res_dict = Dict{affine_patch_type(X), Vector{poly_type(affine_patch_type(X))}}()
#  center_dict = Dict{affine_patch_type(X), Vector{poly_type(affine_patch_type(X))}}()
#  for U in patches(C)
#    if check
#      is_equidimensional_and_smooth(U) || error("ambient space is not smooth")
#    end
#    merge!(res_dict, prepare_smooth_center(U, I[U], check=check))
#  end
#  Z = subscheme(I)
#  CZ = default_covering(Z)
#  if check
#    for U in patches(CZ)
#      is_equidimensional_and_smooth(U) || error("center is not smooth")
#    end
#  end
#  return I
#end
#
#function prepare_smooth_center(X::Spec, f::Vector{PolyType}; check::Bool=true) where {PolyType<:MPolyRingElem}
#  X_dict = as_smooth_local_complete_intersection(X, verbose=true)
#  return X_dict
#end
#
#function as_smooth_local_complete_intersection(I::IdealSheaf; check::Bool=true, verbose::Bool=false)
#  X = scheme(I)
#  C = covering(I)
#  SpecType = affine_patch_type(X)
#  PolyType = poly_type(SpecType)
#  res_dict = Dict{SpecType, Tuple{Vector{PolyType}, Vector{PolyType}, Vector{Int}}}()
#  for U in patches(C)
#    merge!(res_dict, as_smooth_local_complete_intersection(U, I[U], check=check, verbose=verbose))
#  end
#  return res_dict
#end
#
#function as_smooth_local_complete_intersection(
#    X::Spec, 
#    f::Vector{PolyType}; 
#    check::Bool=true, 
#    verbose::Bool=false
#  ) where {PolyType}
#  verbose && println("call with $X and $f")
#  R = base_ring(OO(X))
#  g = gens(modulus(OO(X)))
#  Dg = jacobi_matrix(g)
#  d = dim(X)
#  n = nvars(R)
#  verbose && println("generators:")
#  verbose && println(f)
#  verbose && println("jacobian matrix:")
#  if verbose
#    for i in 1:nrows(Dg)
#      println(Dg[i, 1:end])
#    end
#  end
#  X_dict = _as_smooth_lci_rec(X, X, PolyType[], (Int[], Int[]), Vector{Tuple{Int, Int}}(), g, Dg, Dg, d, n, check=check, verbose=verbose)
#  verbose && println("####################### done with the ambient space #######################")
#  SpecType = typeof(X)
#  f_dict = Dict{SpecType, Tuple{Vector{PolyType}, Vector{PolyType}, Vector{Int}}}()
#  for U in keys(X_dict)
#    verbose && println("proceed with finding a regular sequence on $U")
#    g = X_dict[U][2]
#    f_ext = vcat(g, f)
#    Df_ext = jacobi_matrix(f_ext)
#    good = X_dict[U][3]
#    Z = subscheme(U, f)
#    d = dim(Z)
#    merge!(f_dict, _as_smooth_lci_rec(X, Z, X_dict[U][1], (collect(1:length(g)), good), 
#                                      Vector{Tuple{Int, Int}}(), 
#                                      f_ext, Df_ext, Df_ext, d, n, 
#                                      check=check, verbose=verbose)
#          )
#  end
#  return f_dict
#end
#
#function as_smooth_lci_of_cod(X::CoveredScheme, c::Int; check::Bool=true, verbose::Bool=false)
#  C = default_covering(X)
#  SpecType = affine_patch_type(X)
#  PolyType = poly_type(SpecType)
#  res_dict = Dict{SpecType, Tuple{Vector{PolyType}, Vector{PolyType}, Vector{Int}}}()
#  for U in patches(C)
#    merge!(res_dict, as_smooth_lci_of_cod(U, d, check=check, verbose=verbose))
#  end
#  return res_dict
#end
#
#function as_smooth_lci_of_cod(X::Spec, c::Int; check::Bool=true, verbose::Bool=false)
#  R = base_ring(OO(X))
#  f = gens(modulus(OO(X)))
#  Df = jacobi_matrix(f)
#  n = nvars(R)
##  verbose && println("generators:")
##  verbose && println(f)
##  verbose && println("jacobian matrix:")
##  if verbose
##    for i in 1:nrows(Df)
##      println(Df[i, 1:end])
##    end
##  end
##    PolyType = poly_type(X)
##    f_part = PolyType[]
##    for i in 1:c+3
##      push!(f_part, sum([(rand(Int)%1000)*g for g in f]))
##    end
##    @show "first try"
##    while !isempty(degeneracy_locus(X, jacobi_matrix(f_part), c, verbose=true))
##      @show "new try"
##      push!(f_part, dot([rand(Int)%1000 for j in 1:length(f)], f))
##    end
##    return _as_smooth_lci_rec(X, X, poly_type(X)[], (Int[], Int[]), Vector{Tuple{Int, Int}}(), f_part, jacobi_matrix(f_part), n-c, n, check=check, verbose=verbose)
#  return _as_smooth_lci_rec(X, X, poly_type(X)[], (Int[], Int[]), Vector{Tuple{Int, Int}}(), f, Df, Df, n-c, n, check=check, verbose=verbose)
#end
#
#function as_smooth_local_complete_intersection(X::CoveredScheme; check::Bool=true, verbose::Bool=false)
#  C = default_covering(X)
#  SpecType = affine_patch_type(X)
#  PolyType = poly_type(SpecType)
#  res_dict = Dict{SpecType, Tuple{Vector{PolyType}, Vector{PolyType}, Vector{Int}}}()
#  for U in patches(C)
#    merge!(res_dict, as_smooth_local_complete_intersection(U, check=check, verbose=verbose))
#  end
#  return res_dict
#end
#
#function as_smooth_local_complete_intersection(X::Spec; check::Bool=true, verbose::Bool=false)
#  R = base_ring(OO(X))
#  f = gens(modulus(OO(X)))
#  Df = jacobi_matrix(f)
#  d = dim(X)
#  n = nvars(R)
#  verbose && println("generators:")
#  verbose && println(f)
#  verbose && println("jacobian matrix:")
#  if verbose
#    for i in 1:nrows(Df)
#      println(Df[i, 1:end])
#    end
#  end
#  return _as_smooth_lci_rec(X, X, poly_type(X)[], (Int[], Int[]), Vector{Tuple{Int, Int}}(), f, Df, Df, d, n, check=check, verbose=verbose)
#end
#
#### 
## Given a matrix A with polynomial entries, this routine returns a 
## list of equations hₖ and gₖ, and sets of columns Iₖ and rows Jₖ such that 
## the (Iₖ, Jₖ)-minor of A has the prescribed full rank r along the 
## hypersurface complement of hₖ of the zero locus of gₖ on C.
##
## Altogether, the latter sets form a disjoint union of C into locally 
## closed sets. In particular, we can derive an open covering of C 
## from the hypersurface complements of the listed (I, J)-minors.
##
## Input: 
##   C::Spec     the locus that needs to be covered
##   A::Matrix   the matrix that needs to have rank r
##   r::Int      the prescribed rank of the matrix.
##
## Output: 
##   the following lists with index i
##   loc_list::Vector{Vector{MPolyRingElem}} A list of partial
##               factorizations of polynomials that need to
##               be localized in order to arrive at the i-th 
##               patch.
##   row_list::Vector{Vector{Int}} the indices of the rows 
##               of the non-vanishing minor on the i-th patch.
##   column_list::Vector{Vector{Int}}    the indices of the 
##               columns of the non-vanishing minor on the 
##               i-th patch.
#function _non_degeneration_cover(
#    C::SpecType,
#    A::MatrixType,
#    r::Int; 
#    rec_count::Int=0,
#    check::Bool=true,
#    verbose::Bool=false,
#    restricted_rows::Vector{Vector{Int}}=[Int[]],
#    restricted_columns::Vector{Vector{Int}}=[Int[]]
#  ) where {SpecType<:Spec, MatrixType} #TODO: How specific can we be with the matrices?
#  indent_str = prod([ "#" for i in 1:rec_count ]) * " "
#  #verbose && println(indent_str * "call with $(X) and matrix $(A) for rank $r")
#  verbose && println(indent_str * "call to _non_degeneration_cover for rank $r")
#  verbose && println(indent_str * "on $C")
#  verbose && println(indent_str * "for $A")
#  verbose && println(indent_str * "with restricted rows $restricted_rows")
#  verbose && println(indent_str * "and restricted columns $restricted_columns")
#
#  # check for some other aborting conditions
#  m = nrows(A)
#  while length(restricted_rows) > 0 && length(restricted_rows[1]) == 0
#    popfirst!(restricted_rows)
#  end
#  length(restricted_rows) == 0 && (restricted_rows = [collect(1:m)])
#  n = ncols(A)
#  while length(restricted_columns) > 0 && length(restricted_columns[1]) == 0
#    popfirst!(restricted_columns)
#  end
#  length(restricted_columns) == 0 && (restricted_columns = [collect(1:n)])
#  R = base_ring(OO(C))
#  
#  verbose && println("checking for the scheme being empty...")
#  if isempty(C)
#    verbose && println("the scheme is empty; returning empty lists")
#    return (Vector{Vector{poly_type(C)}}(), Vector{Vector{poly_type(C)}}(), Vector{Vector{Int}}(), Vector{Vector{Int}}())
#  end
#
##  if r == 1 
##    verbose && println("reached the rank-1-case")
##    Y = subscheme(C, ideal(R, [A[i,j] for i in 1:nrows(A) for j in 1:ncols(A)]))
##    if check 
##      isempty(Y) || error("the prescribed rank can not be achieved on $Y")
##    end
##    ll = Vector{Vector{poly_type(C)}}()
##    dl = Vector{Vector{poly_type(C)}}()
##    rl = Vector{Vector{Int}}() 
##    cl = Vector{Vector{Int}}()
##    for i in 1:nrows(A)
##      for j in 1:ncols(A)
##        if !iszero(OO(C)(A[i,j]))
##          push!(ll, [A[i,j]])
##          push!(dl, poly_type(C)[])
##          push!(rl, [i])
##          push!(cl, [j])
##        end
##      end
##    end
##    #verbose && println("returning localizations at $ll modulo $dl")
##    return ll, dl, rl, cl
##  end
#
#  # find a suitable entry of lowest degree for the next recursion
#  verbose && println("reducing the entries...")
#  (k, l) = (1, 1)
#  d = maximum(total_degree.(A[restricted_rows[1], restricted_columns[1]]))+1
#  I = localized_modulus(OO(C))
#  allzero = true
#  W = localized_ring(OO(C))
#  verbose && print(indent_str*" reducing row number ")
#  for i in 1:m
#    verbose && print("$i")
#    v = W.([A[i,j] for j in 1:n])
#    verbose && print(".")
#    w = numerator.(reduce(v, groebner_basis(I))) # TODO: Implement and use reduction for matrices?
#    verbose && print(".")
#    for j in 1:n
#      A[i,j] = w[j]
#      if i in restricted_rows[1] && j in restricted_columns[1] 
#        allzero = allzero && iszero(A[i,j])
#        if total_degree(w[j]) <= d && !iszero(w[j])
#          d = total_degree(w[j])
#          (k, l) = (i, j)
#        end
#      end
#    end
#    verbose && print(".")
#  end
#  f = A[k,l]
#  verbose && println("")
#
#  # in case that after reduction all the matrix' entries are zero, quit
#  if r> 0 && allzero
#    error(indent_str * "All entries are zero")
#  end
#
#  verbose && println("oracle gives $(_degree_oracle_rec(total_degree.(A[restricted_rows[1], restricted_columns[1]]), rec_depth=r-1))")
#  A_part = A[restricted_rows[1], restricted_columns[1]]
#  deg_A_part = total_degree.(A[restricted_rows[1], restricted_columns[1]])
#  (i,j), e = _degree_oracle_rec(deg_A_part, rec_depth=r-1)
#  (k,l) = restricted_rows[1][i], restricted_columns[1][j]
#  f = A[k,l]
#  verbose && println(indent_str * "selected entry at ($k, $l): $f")
#
#  verbose && println("preparing the matrix for induction...")
#  B = copy(A)
#  for i in 1:k-1
#    multiply_row!(B, f, i)
#    add_row!(B, -A[i,l], k, i)
#  end
#  for i in k+1:m
#    multiply_row!(B, f, i)
#    add_row!(B, -A[i,l], k, i)
#  end
#
#  # set up the submatrix for induction on the open part
#  u = [i for i in 1:m if i != k]
#  v = [j for j in 1:n if j != l]
#  B_part = B[u, v]
#  new_restricted_rows = [[(i > k ? i-1 : i) for i in a if i != k] for a in restricted_rows]
#  new_restricted_columns = [[(j > l ? j-1 : j) for j in b if j != l] for b in restricted_columns]
#
#  verbose && println("do the induction steps.")
#
#  loc_list = Vector{Vector{poly_type(C)}}()
#  div_list = Vector{Vector{poly_type(C)}}()
#  row_list = Vector{Vector{Int}}() 
#  column_list = Vector{Vector{Int}}()
#
#  if r>1
#    # prepare for recursion on the open part
#    verbose && println("preparing the hypersurface complement")
#    U = hypersurface_complement(C, f, keep_cache=true)
#
#    # harvest from recursion on the open part
#    verbose && println("recursive call for the complement")
#    llU, dlU, rlU, clU = _non_degeneration_cover(U, B_part, r-1, 
#                                                 rec_count=rec_count+1, 
#                                                 check=check, 
#                                                 verbose=verbose,
#                                                 restricted_rows=new_restricted_rows,
#                                                 restricted_columns=new_restricted_columns
#                                                )
#
#    # process the output according to the preparations for the recursion
#    loc_list = [push!(Ul, f) for Ul in llU]
#    div_list = dlU
#
#    # the indices have to be adjusted according to the choice of the submatrix B from A
#    row_list = [push!(Ul, k) for Ul in [[(i < k ? i : i+1) for i in a] for a in rlU]]
#    column_list = [push!(Ul, l) for Ul in [[(i < l ? i : i+1) for i in a] for a in clU]]
#    #verbose && println("return value was non-trivial; composing as $loc_list and $div_list")
#    if check
#      verbose && println("checking the return values")
#      n = length(loc_list)
#      for i in 1:n
#        X = hypersurface_complement(subscheme(C, div_list[i]), prod(loc_list[i]))
#        D = A[row_list[i], column_list[i]]
#        g = det(D)
#        isunit(OO(X)(g)) || error("selected minor is not a unit")
#      end
#    end
#  else
#    push!(loc_list, [f])
#    push!(div_list, poly_type(C)[])
#    push!(row_list, [k])
#    push!(column_list, [l])
#  end
#
#  # prepare for recursion on the closed part
#  verbose && println("preparing the hypersurface subscheme")
#  Y = subscheme(C, f)
#
#  # harvest from recursion on the closed part
#  verbose && println("recursive call for the subscheme")
#  llY, dlY, rlY, clY = _non_degeneration_cover(Y, copy(A), r, 
#                                               rec_count=rec_count+1, 
#                                               check=check, 
#                                               verbose=verbose,
#                                               restricted_rows=restricted_rows,
#                                               restricted_columns=restricted_columns
#                                              )
#  
#  # process the output according to the preparations for the recursion
#  loc_list = vcat(loc_list, llY)
#  div_list = vcat(div_list, [push!(Ud, f) for Ud in dlY])
#  # no adjustment necessary, since this happens on the same recursion level
#  row_list = vcat(row_list, rlY)
#  column_list = vcat(column_list, clY)
#  #verbose && println("return value was non-trivial; adding $llY and $dlY")
#
#  return loc_list, div_list, row_list, column_list
#end
#    
#function _as_smooth_lci_rec(
#    X::SpecType, 
#    Z::SpecType, # the uncovered locus
#    h::Vector{PolyType}, # equations that were localized already
#    good::Tuple{Vector{Int}, Vector{Int}}, # partial derivatives that form the beginning of the solution sequence
#    bad::Vector{Tuple{Int, Int}}, # partial derivatives that vanish on the uncovered locus.
#    f::Vector{PolyType}, 
#    Df::MatrixType,
#    B::MatrixType,
#    d::Int, n::Int;
#    check::Bool=true,
#    verbose::Bool=false,
#    rec_depth::Int=0
#  ) where{
#          SpecType<:Spec,
#          PolyType<:MPolyRingElem,
#          MatrixType
#         }
#  recstring = prod(["#" for i in 0:rec_depth])
#  #verbose && println(recstring * "call with $X, $Z")
#  verbose && println(recstring * "selected minors: $(good[1]) x $(good[2])")
#  verbose && println(recstring * "bad positions: $bad")
#  # return format: 
#  # key: affine patch
#  # value: (equations that were localized from root, 
#  #         regular sequence, 
#  #         index of variables for non-vanishing minor)
#  res_dict = Dict{typeof(X), Tuple{Vector{PolyType}, Vector{PolyType}, Vector{Int}}}()
#
#  if isempty(Z)
#    verbose && println("this patch is already covered by others")
#    return res_dict
#  end
#
#  if length(good[1]) == n-d
#    verbose && println(recstring * "end of recursion reached")
#    #g = det(Df[good[1], good[2]])
#    #isempty(subscheme(Z, g)) || error("scheme is not smooth")
#    res_dict[X] = (h, f[good[2]], good[1])
#    if check
#      issubset(localized_modulus(OO(X)), ideal(OO(X), f[good[2]])) || error("the found generators are not correct")
#    end
#    verbose && println(recstring*"returning $X with generators $(res_dict[X])")
#    return res_dict
#  end
#
#  verbose && println("setting up matrix of minors")
#  R = base_ring(OO(X))
#  W = localized_ring(OO(X))
#  J = localized_modulus(OO(Z))
#  m = ncols(Df)
#  n = nrows(Df)
#  for i in [a for a in 1:n if !(a in good[1])]
#    for j in [b for b in 1:m if !(b in good[2])]
#      if !((i,j) in bad)
#        #B[i,j] = numerator(reduce(W(B[i,j]), groebner_basis(localized_modulus(OO(X)))))
#        #B[i,j] = numerator(reduce(W(det(Df[vcat(good[1], i), vcat(good[2], j)])), groebner_basis(localized_modulus(OO(X)))))
#        B[i,j] = det(Df[vcat(good[1], i), vcat(good[2], j)])
#      end
#    end
#  end
#  min = maximum([total_degree(a) for a in collect(B)])
#  (k, l) = (0, 0)
#  verbose && println("reducing matrix and selecting pivot...")
#  for i in [a for a in 1:n if !(a in good[1])]
#    for j in [b for b in 1:m if !(b in good[2])]
#      if !((i,j) in bad)
#        #B[i,j] = numerator(reduce(W(B[i,j]), groebner_basis(localized_modulus(OO(X)))))
#        #B[i,j] = numerator(reduce(W(det(Df[vcat(good[1], i), vcat(good[2], j)])), groebner_basis(localized_modulus(OO(X)))))
#        #Df[i,j] = normal_form(Df[i,j], groebner_basis(saturated_ideal(localized_modulus(OO(X)))))
#        if total_degree(B[i, j]) <= min && !iszero(OO(Z)(B[i, j]))
#          min = total_degree(B[i, j])
#          (k, l) = (i, j)
#        end
#      end
#    end
#  end
##  if verbose
##    for i in 1:nrows(Df)
##      println(Df[i, 1:end])
##    end
##  end
#  if (k, l) == (0, 0)
#    isempty(Z) && return res_dict
#    error("scheme is not smooth")
#  end
#  verbose && println(recstring*"selected the $((k, l))-th entry: $(Df[k,l])")
#  
#  good_ext = (vcat(good[1], k), vcat(good[2], l))
#  p = Df[k,l]
#  Bnew = copy(B)
## for i in [a for a in 1:n if !(a in good_ext[1])]
##   if !iszero(Bnew[i,l])
##     @show i
##     @show n
##     multiply_row!(Bnew, p, i)
##     add_row!(Bnew, -B[i,l], k, i)
##   end
## end
##   for j in [b for b in 1:m if !(b in good_ext[2])]
##     if !iszero(Bnew[k,j])
##       @show j
##       multiply_column!(Bnew, p, j)
##       @show m
##       add_column!(Bnew, -Df[k,j], j, l)
##     end
##   end
#  @show 1
#  g = det(Df[good_ext[1], good_ext[2]])
#  g = B[k,l]
#  #g = B[k,l]
#  @show total_degree(g)
#  @show has_attribute(localized_modulus(OO(X)), :saturated_ideal)
#  @show has_attribute(localized_modulus(OO(Z)), :saturated_ideal)
#  I = localized_modulus(OO(Z))
#  Isat = saturated_ideal(I)
#  @show has_attribute(I, :saturated_ideal)
#  @show 2
#  U = hypersurface_complement(X, g, keep_cache=false)
#  @show 2.5
#  Z_next = hypersurface_complement(Z, g, keep_cache=true)
##  U = X
##  Z_next = Z
##  for a in factor(g)
##    @show a
##    U = hypersurface_complement(U, a[1])
##    Z_next = hypersurface_complement(Z_next, a[1])
##  end
#  @show 3
#  U_dict = _as_smooth_lci_rec(U, Z_next, vcat(h, [g]), good_ext, bad, f, Df, B, d, n, check=check, verbose=verbose, rec_depth=rec_depth+1)
#  @show 4
#  bad_ext = vcat(bad, (k, l))
#  @show 5
#  V_dict = _as_smooth_lci_rec(X, subscheme(Z, g), h, good, bad_ext, f, Df, B, d, n, check=check, verbose=verbose, rec_depth=rec_depth+1)
#  @show 6
#  return merge!(U_dict, V_dict)
#end
#
#
#
#function weak_transform(f::CoveredSchemeMorphism, E::IdealSheaf, I::IdealSheaf)
#  X = domain(f)
#  Y = codomain(f)
#  CX = domain_covering(f)
#  CY = codomain_covering(f)
#  SpecType = affine_patch_type(X)
#  trans_dict = Dict{SpecType, Vector{poly_type(SpecType)}}()
#  for U in patches(CX)
#    trans_dict[U] = weak_transform(f[U], E[U], I[codomain(f[U])])
#  end
#  return IdealSheaf(X, CX, trans_dict, check=true)
#end
#
#function controlled_transform(f::CoveredSchemeMorphism, E::IdealSheaf, I::IdealSheaf, i::Int)
#  X = domain(f)
#  Y = codomain(f)
#  CX = domain_covering(f)
#  CY = codomain_covering(f)
#  SpecType = affine_patch_type(X)
#  trans_dict = Dict{SpecType, Vector{poly_type(SpecType)}}()
#  for U in patches(CX)
#    trans_dict[U] = controlled_transform(f[U], E[U], I[codomain(f[U])],i)
#  end
#  return IdealSheaf(X, CX, trans_dict, check=true)
#end
#
#function test_cover(C, f, hl, ql, rl, cl)
#  n = length(hl)
#  Df = jacobi_matrix(f)
#  for i in 1:n
#    h = prod(hl[i])
#    A = Df[rl[i], cl[i]]
#    g = det(A)
#    U = hypersurface_complement(subscheme(C, ql[i]), h)
#    @show isunit(OO(U)(g))
#  end
#
#
#  I = ideal(OO(C), [det(Df[rl[i], cl[i]]) for i in 1:length(rl)])
#  @show one(localized_ring(OO(C))) in I
#end
#
#function _degree_oracle_rec(M::Matrix{Int}; rec_depth::Int=0)
#  m = nrows(M)
#  n = ncols(M)
#  d = maximum(M)
#  minimal_entries = Vector{Tuple{Int, Int}}()
#  for i in 1:m
#    for j in 1:n
#      if M[i,j] >= 0 && M[i,j] < d
#        d = M[i,j]
#        minimal_entries = [(i,j)]
#      end
#      if M[i,j] >=0 && M[i,j] <= d
#        push!(minimal_entries, (i, j))
#      end
#    end
#  end
#
#  if rec_depth==0
#    return minimal_entries[1], d
#  end
#
#  oracles = Vector{Tuple{Tuple{Int, Int}, Int}}()
#  for p in 1:length(minimal_entries)
#    (k,l) = minimal_entries[p]
#    Mnew = Matrix{Int}(undef, m-1, n-1)
#    for i in 1:k-1
#      for j in 1:l-1
#        Mnew[i,j] = maximum([(M[i,j] >= 0 ? M[i,j] + M[k,l] : -1), (M[i,l]>=0 ? M[k,j] + M[i,l] : -1)])
#      end
#      for j in l+1:n
#        Mnew[i,j-1] = maximum([(M[i,j] >= 0 ? M[i,j] + M[k,l] : -1), (M[i,l]>=0 ? M[k,j] + M[i,l] : -1)])
#      end
#    end
#    for i in k+1:m
#      for j in 1:l-1
#        Mnew[i-1,j] = maximum([(M[i,j] >= 0 ? M[i,j] + M[k,l] : -1), (M[i,l]>=0 ? M[k,j] + M[i,l] : -1)])
#      end
#      for j in l+1:n
#        Mnew[i-1,j-1] = maximum([(M[i,j] >= 0 ? M[i,j] + M[k,l] : -1), 
#                                 (M[i,l]>=0 ? M[k,j] + M[i,l] : -1)])
#      end
#    end
#    push!(oracles, _degree_oracle_rec(Mnew))
#  end
#
#  e = minimum([a[2] for a in oracles])
#  for p in 1:length(minimal_entries)
#    if oracles[p][2] == e
#      return minimal_entries[p], maximum([d, e])
#    end
#  end
#  return minimal_entries[p], maximum([d, e])
#end
#
#

########################################################################
# Printing                                                             #
########################################################################
function Base.show(io::IO, X::CoveredProjectiveScheme)
  print(io, "relative projective scheme over $(base_scheme(X))")
end

