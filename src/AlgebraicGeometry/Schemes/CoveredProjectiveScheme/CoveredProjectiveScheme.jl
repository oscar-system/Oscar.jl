export CoveredProjectiveScheme
export ProjectiveGluing
export base_scheme
export blow_up
export controlled_transform
export covered_scheme
export projective_patches
export strict_transform
export weak_transform



### getters for the essential functionality
base_gluing(PG::AbsProjectiveGluing) = base_gluing(underlying_gluing(PG))
inclusion_maps(PG::AbsProjectiveGluing) = inclusion_maps(underlying_gluing(PG))
gluing_domains(PG::AbsProjectiveGluing) = gluing_domains(underlying_gluing(PG))
patches(PG::AbsProjectiveGluing) = patches(underlying_gluing(PG))
gluing_morphisms(PG::AbsProjectiveGluing) = gluing_morphisms(underlying_gluing(PG))

###############################################################################
#
#  Printing
#
###############################################################################


# Essential getters
patches(G::LazyProjectiveGluing) = G.patches
base_gluing(G::LazyProjectiveGluing) = G.base_gluing

# Everything else will trigger the computation
function underlying_gluing(G::LazyProjectiveGluing)
  if !isdefined(G, :underlying_gluing)
    G.underlying_gluing = G.compute_function(G.gluing_data)
  end
  return G.underlying_gluing
end

function Base.show(io::IO, PG::LazyProjectiveGluing)
  print(io, "Gluing of projective patches (not yet computed)")
end

function Base.show(io::IO, PG::ProjectiveGluing)
  print(io, "Gluing of projective patches")
end

function Base.show(io::IO, ::MIME"text/plain", PG::ProjectiveGluing)
  io = pretty(io)
  f = gluing_morphisms(PG)[1]
  PX, PY = patches(PG)
  PU, PV = gluing_domains(PG)
  println(io, "Gluing")
  println(io, Indent(), "of  ", Lowercase(), PX)
  println(io, "and ", Lowercase(), PY)
  println(io, Dedent(), "along the open subsets")
  println(io, Indent(), Lowercase(), PU)
  println(io, Lowercase(), PV)
  print(io, Dedent(), "defined by ", Lowercase())
  show(IOContext(io, :show_semi_compact => true), f)
end

### type getters
#=
TODO: Do we need these?
gluing_type(P::T) where {T<:ProjectiveScheme} = ProjectiveGluing{gluing_type(base_scheme_type(T)), T, morphism_type(T)}
gluing_type(::Type{T}) where {T<:ProjectiveScheme} = ProjectiveGluing{gluing_type(base_scheme_type(T)), T, morphism_type(T)}
=#
### essential getters

base_gluing(PG::ProjectiveGluing) = PG.G
inclusion_maps(PG::ProjectiveGluing) = (PG.inc_to_P, PG.inc_to_Q)
gluing_domains(PG::ProjectiveGluing) = (domain(PG.f), domain(PG.g))
patches(PG::ProjectiveGluing) = (codomain(PG.inc_to_P), codomain(PG.inc_to_Q))
gluing_morphisms(PG::ProjectiveGluing) = (PG.f, PG.g)


base_scheme(P::CoveredProjectiveScheme) = P.Y
base_covering(P::CoveredProjectiveScheme) = P.BC
base_patches(P::CoveredProjectiveScheme) = patches(P.BC)
projective_patches(P::CoveredProjectiveScheme) = values(P.patches)
getindex(P::CoveredProjectiveScheme, U::AbsAffineScheme) = (P.patches)[U]
getindex(P::CoveredProjectiveScheme, U::AbsAffineScheme, V::AbsAffineScheme) = (P.gluings)[(U, V)]

###############################################################################
#
#  Printing
#
###############################################################################

function Base.show(io::IO, CPS::CoveredProjectiveScheme)
  io = pretty(io)
  n = length(projective_patches(CPS))
  K = base_ring(base_scheme(CPS))
  if is_terse(io)
    print(io, "Relative projective scheme")
  else
    if length(projective_patches(CPS)) == 0
      print(io, "Empty relative projective scheme over ")
    else
      print(io, "Relative projective scheme over ")
    end
    print(io, Lowercase(), base_scheme(CPS))
    if n != 0
      print(io, " covered with ", ItemQuantity(n, "projective patch"))
    end
  end
end

function Base.show(io::IO, ::MIME"text/plain", CPS::CoveredProjectiveScheme)
  io = pretty(io)
  pp = projective_patches(CPS)
  n = length(pp)
  println(io, "Relative projective scheme")
  print(io, Indent(), "over ", Lowercase())
  show(IOContext(io, :show_semi_compact => true, :covering => base_covering(CPS)), base_scheme(CPS))
  println(io)
  print(io, Dedent(), "covered with ", ItemQuantity(n, "projective patch"))
  print(io, Indent())
  l = ndigits(n)
  for (i, ppp) in enumerate(pp)
    li = ndigits(i)
    println(io)
    print(io, " "^(l-li)*"$(i): ", Lowercase(), ppp)
  end
  print(io, Dedent())
end

@doc raw"""
    empty_covered_projective_scheme(R::T) where T <: Ring
                                                  -> CoveredProjectiveScheme{T}

Given a ring `R`, return the empty relative projective scheme over the
empty covered scheme over `R`.

# Examples
```jldoctest
julia> R, (x,y,z) = QQ[:x, :y, :z];

julia> Oscar.empty_covered_projective_scheme(R)
Relative projective scheme
  over empty covered scheme over R
covered with 0 projective patches
```
"""
function empty_covered_projective_scheme(R::T) where {T<:AbstractAlgebra.Ring}
  Y = empty_covered_scheme(R)
  C = default_covering(Y)
  #U = C[1]
  #ST = affine_patch_type(Y)
  pp = IdDict{AbsAffineScheme, AbsProjectiveScheme}()
  #P = projective_space(U, 0)
  #pp[U] = P
  tr = IdDict{Tuple{AbsAffineScheme, AbsAffineScheme}, AbsProjectiveGluing}()
  #W = AffineSchemeOpenSubscheme(U)
  #PW, inc = fiber_product(restriction_map(U, W), P)
  #tr[(U, U)] = ProjectiveGluing(Gluing(U, U, identity_map(W), identity_map(W)),
                                 #inc, inc, identity_map(PW), identity_map(PW))
  return CoveredProjectiveScheme(Y, C, pp, tr)
end

@doc raw"""
    blow_up_chart(W::AbsAffineScheme, I::Ideal)

Return the blowup of ``W`` at the ideal ``I``; this is a `ProjectiveScheme`
with `base_scheme` ``W``.

!!! note
    blow_up relies on this internal method for computing the blow ups of all chartsand appropriately assembles the returned projective schemes to a single coverec scheme.
"""

function blow_up_chart(W::AbsAffineScheme, I::Ideal; var_name::VarName = :s)
  error("method `blow_up_chart` not implemented for arguments of type $(typeof(W)) and $(typeof(I))")
end

########################################################################
# Blowups of affine schemes                                            #
########################################################################

function blow_up_chart(W::AbsAffineScheme{<:Field, <:MPolyRing}, I::MPolyIdeal;
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
    J = ideal(S, [t[i]*g[j] - t[j]*g[i] for i in 1:r for j in i+1:r+1])
    IPY = subscheme(IPW, J)
    # Compute the IdealSheaf for the exceptional divisor
    ID = IdDict{AbsAffineScheme, RingElem}()
    Y = covered_scheme(IPY)
    p = covered_projection_to_base(IPY)
    p_cov = covering_morphism(p)
    for i in 1:ngens(I)
      U = affine_charts(Y)[i]
      p_res = p_cov[U]
      W === codomain(p_res) || error("codomain not correct")
      ID[affine_charts(Y)[i]] = pullback(p_res)(gen(I, i))
    end
    E = Oscar.EffectiveCartierDivisor(Y, ID, trivializing_covering=domain(p_cov), check=false)
    set_attribute!(Y, :exceptional_divisor, E)
    set_attribute!(IPY, :exceptional_divisor, E)

    # Prepare the decomposition data
    decomp_dict = IdDict{AbsAffineScheme, Vector{RingElem}}()
    for k in 1:ngens(I)
      U = affine_charts(Y)[k]
      decomp_dict[U] = gens(OO(U))[1:k-1] # Relies on the projective variables coming first!
    end

    # Cache the isomorphism on the complement of the center
    p_res_dict = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
    for i in 1:ngens(I)
      UW = PrincipalOpenSubset(W, gen(I, i))
      V = affine_charts(Y)[i]
      VW = PrincipalOpenSubset(V, E(V))
      p_res_dict[VW] = restrict(p_cov[V], VW, UW, check=false)
      g = OO(UW).(gens(I))
      set_attribute!(p_res_dict[VW], :inverse,
                     morphism(UW, VW, vcat([g[j]*inv(g[i]) for j in 1:ngens(I) if j != i], gens(OO(UW))), check=false)
                    )
    end
    set_attribute!(p, :isos_on_complement_of_center, p_res_dict)
    set_decomposition_info!(default_covering(Y), decomp_dict)
    return IPY
  else
    # construct the blowup by elimination.
    CIPW, pullback_to_cone = affine_cone(IPW)
    R = base_ring(I)
    x = gens(R)
    kk = coefficient_ring(R)
    A, x_ext = polynomial_ring(kk, vcat(symbols(R), [:t]), cached=false)
    t = last(x_ext)
    inc = hom(R, A, x_ext[1:end-1], check=false)
    phi = hom(OO(CIPW), A, vcat([inc(g[i])*t for i in 1:r+1], x_ext[1:end-1], ), check=false) # the homogeneous variables come first
    J = kernel(phi)
    pb = inverse(pullback_to_cone)
    Jh = ideal(homogeneous_coordinate_ring(IPW), pb.(lifted_numerator.(gens(J))))
    IPY = subscheme(IPW, Jh)
    # Compute the IdealSheaf for the exceptional divisor
    ID = IdDict{AbsAffineScheme, RingElem}()
    Y = covered_scheme(IPY)
    p = covered_projection_to_base(IPY)
    p_cov = covering_morphism(p)
    for i in 1:ngens(I)
      U = affine_charts(Y)[i]
      p_res = p_cov[U]
      W === codomain(p_res) || error("codomain not correct")
      ID[affine_charts(Y)[i]] = pullback(p_res)(gen(I, i))
    end
    E = Oscar.EffectiveCartierDivisor(Y, ID, trivializing_covering=domain(p_cov), check=false)
    set_attribute!(Y, :exceptional_divisor, E)
    set_attribute!(IPY, :exceptional_divisor, E)
    # Cache the isomorphism on the complement of the center
    p_res_dict = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
    for i in 1:ngens(I)
      UW = PrincipalOpenSubset(W, gen(I, i))
      V = affine_charts(Y)[i]
      VW = PrincipalOpenSubset(V, E(V))
      p_res_dict[VW] = restrict(p_cov[V], VW, UW, check=false)
      g = OO(UW).(gens(I))
      set_attribute!(p_res_dict[VW], :inverse,
                     morphism(UW, VW, vcat([g[j]*inv(g[i]) for j in 1:ngens(I) if j != i], gens(OO(UW))), check=false)
                    )
    end

    # Prepare the decomposition data
    decomp_dict = IdDict{AbsAffineScheme, Vector{RingElem}}()
    for k in 1:ngens(I)
      U = affine_charts(Y)[k]
      decomp_dict[U] = gens(OO(U))[1:k-1] # Relies on the projective variables coming first!
    end

    set_decomposition_info!(default_covering(Y), decomp_dict)
    set_attribute!(p, :isos_on_complement_of_center, p_res_dict)
    return IPY
  end
end

function blow_up_chart(W::AbsAffineScheme{<:Field, <:RingType}, I::Ideal;
    var_name::VarName = :s
  ) where {RingType<:Union{MPolyQuoRing, MPolyLocRing, MPolyQuoLocRing}}
  base_ring(I) === OO(W) || error("ideal does not belong to the correct ring")

  # It follows the generic Proj construction
  R = OO(W)
  T, (t,) = polynomial_ring(R, [:t], cached=false)
  S, s = graded_polynomial_ring(R, [Symbol(var_name, i-1) for i in 1:ngens(I)])
  phi = hom(S, T, [t*g for g in gens(I)], check=false)
  K = kernel(phi)
  K = ideal(S, [g for g in gens(K) if !iszero(g)]) # clean up superfluous generators
  Bl_W = proj(S, K)
  set_base_scheme!(Bl_W, W)
  # Compute the IdealSheaf for the exceptional divisor
  ID = IdDict{AbsAffineScheme, RingElem}()
  Y = covered_scheme(Bl_W)
  p = covered_projection_to_base(Bl_W)
  p_cov = covering_morphism(p)
  for i in 1:ngens(I)
    @assert !iszero(gen(I, i))
    U = affine_charts(Y)[i]
    p_res = p_cov[U]
    W === codomain(p_res) || error("codomain not correct")
    ID[affine_charts(Y)[i]] = pullback(p_res)(gen(I, i))
  end
  E = Oscar.EffectiveCartierDivisor(Y, ID, trivializing_covering=domain(p_cov), check=false)
  set_attribute!(Y, :exceptional_divisor, E)
  set_attribute!(Bl_W, :exceptional_divisor, E)

  # Cache the isomorphism on the complement of the center
  p_res_dict = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
  for i in 1:ngens(I)
    UW = PrincipalOpenSubset(W, gen(I, i))
    V = affine_charts(Y)[i]
    VW = PrincipalOpenSubset(V, E(V))
    p_res_dict[VW] = restrict(p_cov[V], VW, UW, check=false)
    g = OO(UW).(gens(I))
    set_attribute!(p_res_dict[VW], :inverse,
                   morphism(UW, VW, vcat([g[j]*inv(g[i]) for j in 1:ngens(I) if j != i], gens(OO(UW))), check=false)
                  )
  end
  set_attribute!(p, :isos_on_complement_of_center, p_res_dict)

  # Prepare the decomposition data
  decomp_dict = IdDict{AbsAffineScheme, Vector{RingElem}}()
  for k in 1:ngens(I)
    U = affine_charts(Y)[k]
    decomp_dict[U] = gens(OO(U))[1:k-1] # Relies on the projective variables coming first!
  end
  set_decomposition_info!(default_covering(Y), decomp_dict)

  return Bl_W
end

function is_regular_sequence(g::Vector{T}) where {T<:RingElem}
  length(g) == 0 && return true
  R = parent(g[1])
  all(x->parent(x)===R, g) || error("elements do not belong to the correct ring")
  is_unit(g[1]) && return false # See Bruns-Herzog: Cohen-Macaulay rings, section 1.1.
  is_zero_divisor(g[1]) && return false
  A, p = quo(R, ideal(R, g[1]))
  red_seq = elem_type(A)[p(f) for f in g[2:end]]
  return is_regular_sequence(red_seq)
end


#function blow_up(W::AffineScheme, I::MPolyQuoLocalizedIdeal;
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
#  WA1, pW, pA = product(W, affine_space(base_ring(base_ring(OO(W))), 1, variable_name="t#"))
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
#    W::AffineScheme,
#    I::Vector{RingElemType};
#    var_names::Vector{Symbol}=Vector{Symbol}(),
#    verbose::Bool=false,
#    check::Bool=true,
#    is_regular_sequence::Bool=false
#  ) where {RingElemType<:MPolyRingElem}
#
#  # some internal function
#  function _add_variables(R::RingType, v::Vector{Symbol}) where {RingType<:MPolyRing}
#    ext_R, _ = polynomial_ring(coefficient_ring(R), vcat(symbols(R), v), cached=false)
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
#    #@show vcat([t*embeddingAt(f) for f in I], gens(At)[1:end-1])
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
#    AffineSchemeType = affine_patch_type(covered_ambient)
#    PolyType = poly_type(AffineSchemeType)
#    Idict = Dict{AffineSchemeType, ideal_type(ring_type(AffineSchemeType))}()
#    @show I
#    for i in 1:n_patches(C)
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
#    proj_dict = Dict{AffineSchemeType, morphism_type(AffineSchemeType, AffineSchemeType)}()
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

# This is a sample for how to use LazyProjectiveGluings.
# Originally, we probably had some body of a double for-loop iterating over
# pairs of patches (P, Q) that we need to glue. We take that body which
# computes the gluing of P and Q and move
# it to an external function (here _compute_projective_gluing). Then we
# go through all the local variables in the body of the for-loop which
# are needed for the actual computation and create a tailor-made struct to
# house them. We add an extraction section in the beginning of the compute
# function to restore them and recreate the original setting within the
# for-loops.
#
# Finally, we can replace the inner part of the double for-loop with the
# actual wrap-up of the local variables and feed everything to a constructor
# for a `LazyProjectiveGluing` as documented above.
struct CoveredProjectiveGluingData
  U::AbsAffineScheme
  V::AbsAffineScheme
  P::AbsProjectiveScheme
  Q::AbsProjectiveScheme
  G::AbsGluing
  I::AbsIdealSheaf
end

function _compute_projective_gluing(gd::CoveredProjectiveGluingData)
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
  UV, VU = gluing_domains(G)
  f, g = gluing_morphisms(G)

  # to construct the identifications of PUV with QVU we need to
  # express the generators of I(U) in terms of the generators of I(V)
  # on the overlap U ∩ V.
  !(G isa Gluing) || error("method not implemented for this type of gluing")

  QVU, QVUtoQ = fiber_product(OX(V, VU), Q)
  PUV, PUVtoP = fiber_product(OX(U, UV), P)
  # The problem is that on a AffineSchemeOpenSubscheme U ∩ V
  # despite I(U)|U ∩ V == I(V)|U ∩ V, we
  # have no method to find coefficients aᵢⱼ such that fᵢ = ∑ⱼaᵢⱼ⋅gⱼ
  # for the generators fᵢ of I(U) and gⱼ of I(V): Even though
  # we can do this locally on the patches of a AffineSchemeOpenSubscheme, the result
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
  #A = [coordinates(OX(U, VU)(f), I(VU)) for f in gens(I(U))] # A[i][j] = aᵢⱼ
  A = [coordinates(OX(U, VU)(f), ideal(OO(VU), OX(V, VU).(gens(I(V))))) for f in gens(I(U))] # A[i][j] = aᵢⱼ
  #B = [coordinates(OX(V, UV)(g), I(UV)) for g in gens(I(V))] # B[j][i] = bⱼᵢ
  B = [coordinates(OX(V, UV)(g), ideal(OO(UV), OX(U, UV).(gens(I(U))))) for g in gens(I(V))] # B[j][i] = bⱼᵢ
  SQVU = homogeneous_coordinate_ring(QVU)
  SPUV = homogeneous_coordinate_ring(PUV)
  # the induced map is ℙ(UV) → ℙ(VU), tⱼ ↦ ∑ᵢ bⱼᵢ ⋅ sᵢ
  # and ℙ(VU) → ℙ(UV), sᵢ ↦ ∑ⱼ aᵢⱼ ⋅ tⱼ
  fup = ProjectiveSchemeMor(PUV, QVU, hom(SQVU, SPUV, pullback(f), [sum([B[j][i]*SPUV[i] for i in 1:ngens(SPUV)]) for j in 1:length(B)], check=false), check=false)
  gup = ProjectiveSchemeMor(QVU, PUV, hom(SPUV, SQVU, pullback(g), [sum([A[i][j]*SQVU[j] for j in 1:ngens(SQVU)]) for i in 1:length(A)], check=false), check=false)
  return ProjectiveGluing(G, PUVtoP, QVUtoQ, fup, gup, check=false)
end


@doc raw"""
    blow_up(X::AbsAffineScheme, I::Ideal)

Return the blow-up morphism of blowing up ``X`` at ``I`` in ``OO(X)``.
"""
function blow_up(
    X::AbsAffineScheme{<:Any, <:MPolyAnyRing},
    I::MPolyAnyIdeal)
  R = OO(X)
  @req R == base_ring(I) "I must be an ideal in the coordinate ring of X"
  Isheaf = IdealSheaf(X, I)
  return blow_up(Isheaf)
end

@doc raw"""
    blow_up(I::AbsIdealSheaf)

Return the blow-up morphism of blowing up of the underlying scheme of ``I``  at ``I``.
"""
function blow_up(
    I::AbsIdealSheaf;
    verbose::Bool=false,
    check::Bool=true,
    var_name::VarName=:s,
    covering::Covering=default_covering(scheme(I))
  )
  X = space(I)
  local_blowups = IdDict{AbsAffineScheme, AbsProjectiveScheme}()
  comp_iso_dict = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
  for U in patches(covering)
    local_blowups[U] = blow_up_chart(U, I(U), var_name=var_name)
    # Gather the information on the isomorphism on the complement
    p = covered_projection_to_base(local_blowups[U])
    isos_on_complement_of_center = get_attribute(p, :isos_on_complement_of_center)::IdDict{<:AbsAffineScheme, <:AbsAffineSchemeMor}
    # manual merge because `merge` does not preserve IdDicts.
    for x in keys(isos_on_complement_of_center)
      comp_iso_dict[x] = isos_on_complement_of_center[x]
    end
    #comp_iso_dict = merge(comp_iso_dict, isos_on_complement_of_center)
  end
  projective_gluings = IdDict{Tuple{AbsAffineScheme, AbsAffineScheme}, AbsProjectiveGluing}()

  # prepare for the projective gluings
  for (U, V) in keys(gluings(covering))
    P = local_blowups[U]
    Q = local_blowups[V]
    G = covering[U, V]
    gd = CoveredProjectiveGluingData(U, V, P, Q, G, I)
    projective_gluings[U, V] = LazyProjectiveGluing(P, Q, G, _compute_projective_gluing, gd)
  end
  Bl_I = CoveredProjectiveScheme(X, covering, local_blowups, projective_gluings)
  Y = covered_scheme(Bl_I)
  # Assemble the exceptional divisor as a Cartier divisor
  ID = IdDict{AbsAffineScheme, RingElem}()
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
  # Store the information for the isomorphism on the complement
  set_attribute!(pr, :isos_on_complement_of_center, comp_iso_dict)
  pr.exceptional_divisor = EffectiveCartierDivisor(Y, ID, trivializing_covering=domain(p_cov), check=false)
  return pr
end

# The following is a wrap-up of local variables necessary to compute gluings
# from morphisms of graded algebras. It will be used to fill in the
# `LazyGluing` together with the compute function following below.
struct ProjectiveGluingData
  down_left::AbsAffineScheme
  down_right::AbsAffineScheme
  up_left::AbsAffineScheme
  up_right::AbsAffineScheme
  down_covering::Covering
  proj::CoveredProjectiveScheme
end

# This function actually computes the gluing with the data extracted
# from the ProjectiveGluingData.
#
# Originally, this was in the body of a constructor. But we want to
# postpone the actual computation, so we wrap up all necessary local
# variables in a `ProjectiveGluingData` and copy-paste the required
# part of code into the body of this function. Together with a header
# which is restoring our local variables.
function _compute_gluing(gd::ProjectiveGluingData)
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
  P = gd.proj
  (A, B) = gluing_domains(C[U, V])
  @assert ambient_coordinate_ring(A) === ambient_coordinate_ring(U)
  @assert ambient_coordinate_ring(B) === ambient_coordinate_ring(V)
  (f, g) = gluing_morphisms(C[U, V])
  (UD, VD) = gluing_domains(P[U, V])
  @assert coefficient_ring(homogeneous_coordinate_ring(UD)) === OO(A)
  @assert coefficient_ring(homogeneous_coordinate_ring(VD)) === OO(B)
  (fup, gup) = gluing_morphisms(P[U, V])
  (incU, incV) = inclusion_maps(P[U, V])
  @assert coefficient_ring(homogeneous_coordinate_ring(codomain(incU))) === OO(U)
  @assert coefficient_ring(homogeneous_coordinate_ring(codomain(incV))) === OO(V)
  S = homogeneous_coordinate_ring(P[U])
  T = homogeneous_coordinate_ring(P[V])
  i = P[U][UW][2]
  j = P[V][VW][2]
  s_i = gen(S, i)
  t_j = gen(T, j)
  AW = affine_charts(covered_scheme(UD))[i]
  # Since caching of polynomial rings is set to `false` by default, 
  # the following does not hold anymore:
  #@assert ambient_coordinate_ring(AW) === ambient_coordinate_ring(UW)
  # This leads to some unintuitive code below using `evaluate` instead of casting.
  BW = affine_charts(covered_scheme(VD))[j]
  hU = dehomogenization_map(UD, AW)(pullback(fup)(pullback(incV)(t_j)))
  hV = dehomogenization_map(VD, BW)(pullback(gup)(pullback(incU)(s_i)))

  # We need to construct the gluing
  #
  #          f', g'
  #   UW ↩ AAW ↔ BBW ↪ VW
  #
  # as follows.
  # From the base gluing we already have that UW differs
  # from AAW by the pullback of the equation `complement_equation(A)`.
  # But then, also `hU` cuts out another locus which needs to be
  # removed before gluing.
  #
  # The same applies to the other side. Finally, we need to track the
  # coordinates through the homogenization-gluing-pullback-dehomogenization
  # machinery to provide the gluing morphisms.

  ptbUD = covered_projection_to_base(UD)
  ptbVD = covered_projection_to_base(VD)
  phi1 = pullback(ptbUD[AW])
  hhU = lifted_numerator(phi1(domain(phi1)(complement_equation(A), check=false)))
  hhU = hhU * lifted_numerator(hU)
  AAW = PrincipalOpenSubset(UW, evaluate(hhU, gens(OO(UW))))
  phi2 = pullback(ptbVD[BW])
  hhV = lifted_numerator(phi2(domain(phi2)(complement_equation(B), check=false)))
  hhV = hhV * lifted_numerator(hV)
  BBW = PrincipalOpenSubset(VW, evaluate(hhV, gens(OO(VW))))

  x = gens(ambient_coordinate_ring(AAW))
  y = gens(ambient_coordinate_ring(BBW))

  xh = homogenization_map(UD, AW).([evaluate(x, gens(OO(AW))) for x in x])
  yh = homogenization_map(VD, BW).([evaluate(y, gens(OO(BW))) for y in y])

  xhh = [(pullback(gup)(pp), pullback(gup)(qq)) for (pp, qq) in xh]
  yhh = [(pullback(fup)(pp), pullback(fup)(qq)) for (pp, qq) in yh]

  phi = dehomogenization_map(VD, BW)
  psi = dehomogenization_map(UD, AW)

  pb_AW_to_AAW = hom(OO(AW), OO(AAW), gens(OO(AAW)), check=false)
  pb_BW_to_BBW = hom(OO(BW), OO(BBW), gens(OO(BBW)), check=false)
  yimgs = [pb_AW_to_AAW(psi(pp))*inv(pb_AW_to_AAW(psi(qq))) for (pp, qq) in yhh]
  ximgs = [pb_BW_to_BBW(phi(pp))*inv(pb_BW_to_BBW(phi(qq))) for (pp, qq) in xhh]
  ff = morphism(AAW, BBW, hom(OO(BBW), OO(AAW), yimgs, check=false), check=false)
  gg = morphism(BBW, AAW, hom(OO(AAW), OO(BBW), ximgs, check=false), check=false)

  return SimpleGluing(UW, VW, ff, gg, check=false)
  pr.exceptional_divisor = EffectiveCartierDivisor(Y, ID, trivializing_covering=domain(p_cov))
  return pr
end

@attr Any function covered_scheme(P::CoveredProjectiveScheme)
  X = base_scheme(P)
  C = base_covering(P)
  new_patches = Vector{AbsAffineScheme}()
  new_gluings = IdDict{Tuple{AbsAffineScheme, AbsAffineScheme}, AbsGluing}()
  projection_dict = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
  parts = IdDict{AbsAffineScheme, AbsCoveredScheme}()
  for U in patches(C)
    parts[U] = Oscar.covered_scheme(P[U])
  end
  # We first assemble a Covering where the preimage of every base
  # patch appears as one connected component
  result_patches = vcat([affine_charts(parts[U]) for U in patches(C)]...)
  result_gluings = IdDict{Tuple{AbsAffineScheme, AbsAffineScheme}, AbsGluing}()
  for U in patches(C)
    GG = gluings(parts[U])
    for (A, B) in keys(GG)
      result_gluings[(A, B)] = GG[(A, B)]
    end
  end
  result_covering = Covering(result_patches, result_gluings, check=false)

  # Now we need to add gluings
  for (U, V) in keys(gluings(C))
    U === V && continue
    for UW in affine_charts(parts[U])
      for VW in affine_charts(parts[V])
        GGG = LazyGluing(UW, VW, 
                         ProjectiveGluingData(U, V, UW, VW, C, P)
                        )
        result_gluings[(UW, VW)] = GGG
      end
    end
  end
  result_covering = Covering(result_patches, result_gluings, check=false)
  result = CoveredScheme(result_covering)

  projection_dict = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
  for U in patches(C)
    PP = P[U]
    p = covered_projection_to_base(PP)
    cov_mor = covering_morphism(p)
    cov_mor_dict = morphisms(cov_mor)
    for (V, phi) in cov_mor_dict
      projection_dict[V] = phi
    end
  end

  # Assemble the decomposition information if applicable
  decomp_dict = IdDict{AbsAffineScheme, Vector{RingElem}}()
  if has_decomposition_info(C)
    for V in patches(C)
      cov_part = default_covering(parts[V])
      if has_decomposition_info(cov_part)
        for U in patches(cov_part)
          pr = projection_dict[U]
          decomp_dict[U] = vcat(decomposition_info(cov_part)[U],
                                elem_type(OO(U))[pullback(pr)(a) for a in decomposition_info(C)[V]]
                               )
        end
      end
    end
  end
  if length(keys(decomp_dict)) == length(patches(result_covering))
    set_decomposition_info!(result_covering, decomp_dict)
  end

  # TODO: Remove the internal checks in the constructors below
  covering_map = CoveringMorphism(result_covering, C, projection_dict, check=false)
  set_attribute!(P, :covering_projection_to_base, covering_map)

  return result
end

@attr Any function covered_projection_to_base(P::CoveredProjectiveScheme)
  if !has_attribute(P, :covering_projection_to_base)
    covered_scheme(P)
  end
  covering_pr = get_attribute(P, :covering_projection_to_base)::CoveringMorphism
  return CoveredSchemeMorphism(covered_scheme(P), base_scheme(P), covering_pr; check=false)
end


# function strict_transform(f::AffineSchemeMor, h::Vector{PolyType}, g::Vector{PolyType}) where{PolyType<:MPolyRingElem}
#(IOw: Exc Div ^\infty)

#        X = domain(f)
#        Y = codomain(f)
#        R = base_ring(OO(X))
#        Excdiv = ideal(h)
#
#        Pf = pullback(f)
#        Iold = ideal(R, lifted_numerator.(Pf.(g))) + strict_modulus(X)
#
#        while true
#          Inew = quotient(Iold, Excdiv)
#          Iold == Inew && break
#          Iold = Inew
#        end
#        return gens(Iold)
#end
#
#
#
#function total_transform(f::AffineSchemeMor, h::Vector{PolyType}, g::Vector{PolyType}) where{PolyType<:MPolyRingElem}
#        #IOw
#
#        X = domain(f)
#        Y = codomain(f)
#        R = base_ring(OO(X))
#        Excdiv = ideal(h)
#
#        Pf = pullback(f)
#        Iold = ideal(R, lifted_numerator.(Pf.(g))) + strict_modulus(X)
#        return gens(Iold)
#end
#
#
#### NOT TESTED YET
#function weak_transform(f::AffineSchemeMor, h::Vector{PolyType}, g::Vector{PolyType}) where{PolyType<:MPolyRingElem}
#
#        X = domain(f)
#        Y = codomain(f)
#        R = base_ring(OO(X))
#        Excdiv = ideal(h)
#
#        Pf = pullback(f)
#        Iold = ideal(R, lifted_numerator.(Pf.(g))) + strict_modulus(X)
#
#        while true
#           Inew = quotient(Iold, Excdiv)
#           !(Iold == Excdiv * Inew) && break
#           Iold = Inew
#        end
#        return gens(Iold)
#        #(IOw : Exc Div ^k), k maximal
#end
#
#### NOT TESTED YET
#function controlled_transform(f::AffineSchemeMor, h::Vector{PolyType}, g::Vector{PolyType}, i::Int) where{PolyType<:MPolyRingElem}
#        #(IOw : Exc Div ^i)
#
#        X = domain(f)
#        Y = codomain(f)
#        R = base_ring(OO(X))
#        Excdiv = ideal(h)
#
#        Pf = pullback(f)
#        Iold = ideal(R, lifted_numerator.(Pf.(g))) + strict_modulus(X)
#
#
#        for j in 1:i
#          Inew = quotient(Iold,Excdiv)
#          Inew == 1 && break
#          Iold = Inew
#        end
#
#        return gens(Iold)
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
#  AffineSchemeType = affine_patch_type(X)
#  trans_dict = Dict{AffineSchemeType, Vector{poly_type(AffineSchemeType)}}()
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
#  AffineSchemeType = affine_patch_type(X)
#  trans_dict = Dict{AffineSchemeType, Vector{poly_type(AffineSchemeType)}}()
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
#  AffineSchemeType = affine_patch_type(X)
#  trans_dict = Dict{AffineSchemeType, Vector{poly_type(AffineSchemeType)}}()
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
#function prepare_smooth_center(X::AffineScheme, f::Vector{PolyType}; check::Bool=true) where {PolyType<:MPolyRingElem}
#  X_dict = as_smooth_local_complete_intersection(X, verbose=true)
#  return X_dict
#end
#
#function as_smooth_local_complete_intersection(I::IdealSheaf; check::Bool=true, verbose::Bool=false)
#  X = scheme(I)
#  C = covering(I)
#  AffineSchemeType = affine_patch_type(X)
#  PolyType = poly_type(AffineSchemeType)
#  res_dict = Dict{AffineSchemeType, Tuple{Vector{PolyType}, Vector{PolyType}, Vector{Int}}}()
#  for U in patches(C)
#    merge!(res_dict, as_smooth_local_complete_intersection(U, I[U], check=check, verbose=verbose))
#  end
#  return res_dict
#end
#
#function as_smooth_local_complete_intersection(
#    X::AffineScheme,
#    f::Vector{PolyType};
#    check::Bool=true,
#    verbose::Bool=false
#  ) where {PolyType}
#  verbose && println("call with $X and $f")
#  R = base_ring(OO(X))
#  g = gens(modulus(OO(X)))
#  Dg = jacobian_matrix(g)
#  d = dim(X)
#  n = nvars(R)
#  verbose && println("generators:")
#  verbose && println(f)
#  verbose && println("Jacobian matrix:")
#  if verbose
#    for i in 1:nrows(Dg)
#      println(Dg[i, 1:end])
#    end
#  end
#  X_dict = _as_smooth_lci_rec(X, X, PolyType[], (Int[], Int[]), Vector{Tuple{Int, Int}}(), g, Dg, Dg, d, n, check=check, verbose=verbose)
#  verbose && println("####################### done with the ambient space #######################")
#  AffineSchemeType = typeof(X)
#  f_dict = Dict{AffineSchemeType, Tuple{Vector{PolyType}, Vector{PolyType}, Vector{Int}}}()
#  for U in keys(X_dict)
#    verbose && println("proceed with finding a regular sequence on $U")
#    g = X_dict[U][2]
#    f_ext = vcat(g, f)
#    Df_ext = jacobian_matrix(f_ext)
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
#  AffineSchemeType = affine_patch_type(X)
#  PolyType = poly_type(AffineSchemeType)
#  res_dict = Dict{AffineSchemeType, Tuple{Vector{PolyType}, Vector{PolyType}, Vector{Int}}}()
#  for U in patches(C)
#    merge!(res_dict, as_smooth_lci_of_cod(U, d, check=check, verbose=verbose))
#  end
#  return res_dict
#end
#
#function as_smooth_lci_of_cod(X::AffineScheme, c::Int; check::Bool=true, verbose::Bool=false)
#  R = base_ring(OO(X))
#  f = gens(modulus(OO(X)))
#  Df = jacobian_matrix(f)
#  n = nvars(R)
##  verbose && println("generators:")
##  verbose && println(f)
##  verbose && println("Jacobian matrix:")
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
##    while !isempty(degeneracy_locus(X, jacobian_matrix(f_part), c, verbose=true))
##      @show "new try"
##      push!(f_part, dot([rand(Int)%1000 for j in 1:length(f)], f))
##    end
##    return _as_smooth_lci_rec(X, X, poly_type(X)[], (Int[], Int[]), Vector{Tuple{Int, Int}}(), f_part, jacobian_matrix(f_part), n-c, n, check=check, verbose=verbose)
#  return _as_smooth_lci_rec(X, X, poly_type(X)[], (Int[], Int[]), Vector{Tuple{Int, Int}}(), f, Df, Df, n-c, n, check=check, verbose=verbose)
#end
#
#function as_smooth_local_complete_intersection(X::CoveredScheme; check::Bool=true, verbose::Bool=false)
#  C = default_covering(X)
#  AffineSchemeType = affine_patch_type(X)
#  PolyType = poly_type(AffineSchemeType)
#  res_dict = Dict{AffineSchemeType, Tuple{Vector{PolyType}, Vector{PolyType}, Vector{Int}}}()
#  for U in patches(C)
#    merge!(res_dict, as_smooth_local_complete_intersection(U, check=check, verbose=verbose))
#  end
#  return res_dict
#end
#
#function as_smooth_local_complete_intersection(X::AffineScheme; check::Bool=true, verbose::Bool=false)
#  R = base_ring(OO(X))
#  f = gens(modulus(OO(X)))
#  Df = jacobian_matrix(f)
#  d = dim(X)
#  n = nvars(R)
#  verbose && println("generators:")
#  verbose && println(f)
#  verbose && println("Jacobian matrix:")
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
##   C::AffineScheme     the locus that needs to be covered
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
#    C::AffineSchemeType,
#    A::MatrixType,
#    r::Int;
#    rec_count::Int=0,
#    check::Bool=true,
#    verbose::Bool=false,
#    restricted_rows::Vector{Vector{Int}}=[Int[]],
#    restricted_columns::Vector{Vector{Int}}=[Int[]]
#  ) where {AffineSchemeType<:AffineScheme, MatrixType} #TODO: How specific can we be with the matrices?
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
#        is_unit(OO(X)(g)) || error("selected minor is not a unit")
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
#    X::AffineSchemeType,
#    Z::AffineSchemeType, # the uncovered locus
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
#          AffineSchemeType<:AffineScheme,
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
#  min = maximum([total_degree(a) for a in B])
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
#  AffineSchemeType = affine_patch_type(X)
#  trans_dict = Dict{AffineSchemeType, Vector{poly_type(AffineSchemeType)}}()
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
#  AffineSchemeType = affine_patch_type(X)
#  trans_dict = Dict{AffineSchemeType, Vector{poly_type(AffineSchemeType)}}()
#  for U in patches(CX)
#    trans_dict[U] = controlled_transform(f[U], E[U], I[codomain(f[U])],i)
#  end
#  return IdealSheaf(X, CX, trans_dict, check=true)
#end
#
#function test_cover(C, f, hl, ql, rl, cl)
#  n = length(hl)
#  Df = jacobian_matrix(f)
#  for i in 1:n
#    h = prod(hl[i])
#    A = Df[rl[i], cl[i]]
#    g = det(A)
#    U = hypersurface_complement(subscheme(C, ql[i]), h)
#    @show is_unit(OO(U)(g))
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

function fiber_product(
    i1::CoveredClosedEmbedding,
    i2::CoveredClosedEmbedding
  )
  X1 = domain(i1)
  X2 = domain(i2)
  Y = codomain(i1)
  Y === codomain(i2) || error("codomains do not coincide")
  i1_cov = covering_morphism(i1)
  i2_cov = covering_morphism(i2)
  codomain(i1_cov) === codomain(i2_cov) || error("case of different coverings in codomain not implemented")
  cod_cov = codomain(i1_cov)
  #=
  cod_ref, ref1, ref2 = common_refinement(codomain(i1_cov), codomain(i2_cov))
  dom_ref1, i1_res = fiber_product(i1, ref1)
  dom_ref2, i2_res = fiber_product(i1, ref2)
  # etc. etc.... This is roughly the generic code to come.
  =#
  I1 = image_ideal(i1)
  pb_I1 = pullback(i2, I1)
  I2 = image_ideal(i2)
  pb_I2 = pullback(i1, I2)

  j1 = Oscar.CoveredClosedEmbedding(domain(i2), pb_I1)
  Z = domain(j1)
  morphism_dict = IdDict{AbsAffineScheme, ClosedEmbedding}()
  for U in affine_charts(Z)
    V2 = codomain(j1[U])
    W = codomain(i2[V2])
    V1_candidates = maps_with_given_codomain(i1, W)
    @assert length(V1_candidates) == 1 "not the correct number of patches found"
    V1 = domain(first(V1_candidates))
    x = gens(OO(V1))
    lift_x = [preimage(pullback(i1[V1]), f) for f in x]
    pb_x = pullback(i2[V2]).(lift_x)
    pb_x = pullback(j1[U]).(pb_x)
    morphism_dict[U] = ClosedEmbedding(morphism(U, V1, pb_x, check=false), pb_I2(V1), check=false)
  end
  j2_cov = CoveringMorphism(default_covering(Z), domain(i1_cov), morphism_dict, check=false)
  j2 = Oscar.CoveredClosedEmbedding(Z, X1, j2_cov)
  return j1, j2
end
