export ProjectiveGlueing
export glueing_type

export CoveredProjectiveScheme
export base_scheme, base_covering, projective_patches, as_covered_scheme

export blow_up, empty_covered_projective_scheme

export strict_transform, total_transform, weak_transform, controlled_transform

export prepare_smooth_center, as_smooth_local_complete_intersection, as_smooth_lci_of_cod, _non_degeneration_cover

mutable struct ProjectiveGlueing{
                                 GlueingType<:Glueing,
                                 ProjectiveSchemeType<:ProjectiveScheme,
                                 MorphismType<:ProjectiveSchemeMor
                                }
  G::GlueingType # the underlying glueing of the base schemes
  P::ProjectiveSchemeType # the projective scheme over the first patch
  Q::ProjectiveSchemeType # the projective scheme over the second patch
  P_on_U_patches::Vector{MorphismType} # restrictions of these 
  Q_on_V_patches::Vector{MorphismType} # schemes to the patches
  f::Vector{MorphismType} # the glueing maps on the affine
  g::Vector{MorphismType} # patches of the glueing domain 

  ### 
  # Given two relative projective schemes and a glueing 
  #
  #       P           Q
  #     π ↓           ↓ π
  #   G : X ↩ U ≅ V ↪ Y 
  #           ∪   ∪
  #           Uᵢ  Vⱼ
  #
  # together with compatible glueing data 
  #
  #   fᵢ: π⁻¹(Uᵢ) → Q,   gⱼ: π⁻¹(Vⱼ) → P
  #
  # on the affine patches Uᵢ of U (resp. Vⱼ of V), 
  # glue together the projective schemes P and Q
  function ProjectiveGlueing(G::GlueingType, P::T, Q::T, f::Vector{M}, g::Vector{M}; check::Bool=true) where {GlueingType<:Glueing, T<:ProjectiveScheme, M<:ProjectiveSchemeMor}
    (X, Y) = patches(G)
    (U, V) = glueing_domains(G)
    (base_scheme(P) == X && base_scheme(Q) == Y) || error("base glueing is incompatible with the projective schemes")
    P_on_U_patches = Vector{M}()
    Q_on_V_patches = Vector{M}()
    
    length(f) == npatches(U) || error("base glueing is incompatible with the projective schemes")
    # for every patch Uᵢ of U we need to set up and store 
    # the inclusion map of π⁻¹(Uᵢ) ↪ P = π⁻¹(X)
    for i in 1:length(f)
      h = f[i]
      base_scheme(domain(h)) == U[i] || error("base glueing is incompatible with morphisms of projective schemes")
      push!(P_on_U_patches, inclusion_map(domain(h), P))
    end
    
    length(g) == npatches(V) || error("base glueing is incompatible with the projective schemes")
    # the same for V
    for j in 1:length(g)
      h = g[j]
      base_scheme(domain(h)) == V[j] || error("base glueing is incompatible with morphisms of projective schemes")
      push!(Q_on_V_patches, inclusion_map(domain(h), Q))
    end

    # in case a check is required, we assure that the maps on the 
    # patches Uᵢ glue together on the overlaps.
    if check
      for i in 1:length(f)-1
        for j = i+1:length(f)
          W = intersect(U[i], U[j])
          PW = fiber_product(W, P)
          h1 = inclusion_map(PW, domain(P_on_U_patches[i]))
          h2 = inclusion_map(PW, domain(P_on_U_patches[j]))
          compose(h1, P_on_U_patches[i]) == compose(h2, P_on_U_patches[j]) || error("maps do not coincide on overlaps")
        end
      end
      # same for V
      for i in 1:length(g)-1
        for j = i+1:length(g)
          W = intersect(V[i], V[j])
          QW = fiber_product(W, Q)
          h1 = inclusion_map(QW, domain(Q_on_V_patches[i]))
          h2 = inclusion_map(QW, domain(Q_on_V_patches[j]))
          compose(h1, Q_on_V_patches[i]) == compose(h2, Q_on_V_patches[j]) || error("maps do not coincide on overlaps")
        end
      end
    end
    return new{GlueingType, T, M}(G, P, Q, P_on_U_patches, Q_on_V_patches, f, g)
  end
end

### type getters

glueing_type(P::T) where {T<:ProjectiveScheme} = ProjectiveGlueing{glueing_type(base_scheme_type(T)), T, morphism_type(T)}
glueing_type(::Type{T}) where {T<:ProjectiveScheme} = ProjectiveGlueing{glueing_type(base_scheme_type(T)), T, morphism_type(T)}

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
# These structs are designed to accomodate blowups of 
# covered schemes along arbitrary centers, as well as 
# projective bundles. 

@attributes mutable struct CoveredProjectiveScheme{
    BaseSchemeType<:CoveredScheme, 
    CoveringType<:Covering,
    SpecType<:Spec,
    ProjectiveSchemeType<:ProjectiveScheme,
    ProjectiveGlueingType<:ProjectiveGlueing,
    BRT, BRET} <: Scheme{BRT, BRET}
  Y::BaseSchemeType # the base scheme
  BC::CoveringType # the reference covering of the base scheme
  patches::Dict{SpecType, ProjectiveSchemeType} # the projective spaces over the affine patches in the base covering
  glueings::Dict{Tuple{SpecType, SpecType}, ProjectiveGlueingType} # the transitions sitting over the affine patches in the glueing domains of the base scheme

  function CoveredProjectiveScheme(
      Y::BaseSchemeType,
      C::CoveringType,
      projective_patches::Dict{SpecType, ProjectiveSchemeType},
      projective_glueings::Dict{Tuple{SpecType, SpecType}, ProjectiveGlueingType};
      check::Bool=true
    ) where {
             BaseSchemeType<:CoveredScheme, 
             CoveringType<:Covering,
             SpecType<:Spec,
             ProjectiveSchemeType<:ProjectiveScheme,
             ProjectiveGlueingType<:ProjectiveGlueing,
            }
    C in coverings(Y) || error("covering not listed")
    for P in values(projective_patches)
      base_scheme(P) in patches(C) || error("base scheme not found in covering")
    end
    for (U, V) in keys(glueings(C))
      (U, V) in keys(projective_glueings) || error("not all projective glueings were provided")
    end
    return new{BaseSchemeType, CoveringType, SpecType, ProjectiveSchemeType, ProjectiveGlueingType, base_ring_type(SpecType), elem_type(base_ring_type(SpecType))}(Y, C, projective_patches, projective_glueings)
  end
end

base_scheme(P::CoveredProjectiveScheme) = P.Y
base_covering(P::CoveredProjectiveScheme) = P.BC
base_patches(P::CoveredProjectiveScheme) = patches(P.BC)
projective_patches(P::CoveredProjectiveScheme) = values(P.patches)
getindex(P::CoveredProjectiveScheme, U::Spec) = (P.patches)[U]
getindex(P::CoveredProjectiveScheme, U::Spec, V::Spec) = (P.glueings)[(U, V)]

function empty_covered_projective_scheme(R::T) where {T<:AbstractAlgebra.Ring}
  Y = empty_covered_scheme(R)
  C = default_covering(Y)
  U = C[1]
  ST = affine_patch_type(Y)
  pp = Dict{affine_patch_type(Y), projective_scheme_type(affine_patch_type(Y))}()
  P = projective_space(U, 0)
  pp[U] = P
  tr = Dict{Tuple{ST, ST}, Vector{morphism_type(pp[U])}}
  tr[(U, U)] = identity_map(P) 
  return CoveredProjectiveScheme(Y, C, pp, tr)
end

# blow up X in the center described by g using these explicit generators.
function blow_up(
    W::Spec, 
    I::Vector{RingElemType}; 
    var_names::Vector{Symbol}=Vector{Symbol}(), 
    verbose::Bool=false,
    check::Bool=true,
    is_regular_sequence::Bool=false
  ) where {RingElemType<:MPolyElem}

  # some internal function
  function _add_variables(R::RingType, v::Vector{Symbol}) where {RingType<:MPolyRing}
    ext_R, _ = PolynomialRing(coefficient_ring(R), vcat(symbols(R), v))
    n = length(gens(R))
    phi = AlgebraHomomorphism(R, ext_R, gens(ext_R)[1:n])
    return ext_R, phi, gens(ext_R)[(length(gens(R))+1):length(gens(ext_R))]
  end

  A = OO(W)
  R = base_ring(A)
  #TODO: Check if for all i \in I parent(i) == R

  m = length(I)
  Pw = (length(var_names) > 0 ? projective_space(W, var_names) : projective_space(W,m-1))
  S = homog_poly_ring(Pw)

  CP = affine_cone(Pw)
  Polyring = base_ring(OO(CP))
  if !is_regular_sequence
    At, embeddingAt, T =  _add_variables(R,[:t])
    t = T[1]

    #	@show vcat([t*embeddingAt(f) for f in I], gens(At)[1:end-1])
    Phi = AlgebraHomomorphism(Polyring, At, vcat([t*embeddingAt(f) for f in I], gens(At)[1:end-1]))

    Imod = modulus(A)
    IW = ideal(At, [embeddingAt(f) for f in gens(Imod)])
    @show "start groebner basis computation."
    IWpre = preimage(Phi, IW)
    @show "done. Proceed to process"
    #SIWpre = ideal(S,[frac_to_homog(Pw,g) for g in gens(IWpre)])
    SIWpre = ideal(S, poly_to_homog(Pw).(gens(IWpre)))
    @show 0

    projective_version = subscheme(Pw, SIWpre)
    @show 1
    covered_version = as_covered_scheme(projective_version)
    @show 2
    projection_map = get_attribute(projective_version, :covered_projection_to_base)
    @show 3
    E_dict = Dict{affine_patch_type(covered_version), Vector{RingElemType}}()
    for i in 1:length(I)
      @show i
      U = default_covering(covered_version)[i]
      E_dict[U] = [lifted_numerator(pullback(projection_map[U])(I[i]))]
    end
    @show 4
    exc_div = IdealSheaf(covered_version, default_covering(covered_version), E_dict, check=false)
    @show "processing done."
    return projective_version, covered_version, projection_map, exc_div
  else
    M = MatrixSpace(S, 2, ngens(S))
    A = zero(M)
    for i in 1:ngens(S)
      A[1, i] = S[i]
      A[2, i] = I[i]
    end
    IW = ideal(S, minors(A, 2))
    projective_version = subscheme(Pw, IW)
    covered_ambient = as_covered_scheme(Pw)
    ambient_projection_map = get_attribute(Pw, :covered_projection_to_base)
    C = default_covering(covered_ambient)
    SpecType = affine_patch_type(covered_ambient)
    PolyType = poly_type(SpecType)
    Idict = Dict{SpecType, Vector{PolyType}}()
    @show I
    for i in 1:npatches(C)
      @show i
      v = gens(OO(C[i]))
      phi = pullback(ambient_projection_map[C[i]])
      loc_eqns = vcat([v[j]*phi(I[i])-phi(I[j]) for j in 1:i-1], [v[j]*phi(I[i])-phi(I[j+1]) for j in i:length(I)-1])
      @show loc_eqns
      Idict[C[i]] = lifted_numerator.(loc_eqns)
    end
    Itrans = IdealSheaf(covered_ambient, C, Idict, check=false)
    covered_version = subscheme(Itrans)
    set_attribute!(projective_version, :as_covered_scheme, covered_version)
    set_attribute!(projective_version, :standard_covering, default_covering(covered_version))

    proj_dict = Dict{SpecType, morphism_type(SpecType, SpecType)}()
    for i in 1:length(I)
      Z = patches(default_covering(covered_version))[i]
      U = patches(default_covering(covered_ambient))[i]
      proj_dict[Z] = restrict(ambient_projection_map[U], Z, codomain(ambient_projection_map[U]), check=false)
    end

    projection_map = CoveringMorphism(default_covering(covered_version), Covering(W), proj_dict)
    set_attribute!(projective_version, :covered_projection_to_base, projection_map)
    @show 3
    E_dict = Dict{affine_patch_type(covered_version), Vector{RingElemType}}()
    for i in 1:length(I)
      @show i
      U = default_covering(covered_version)[i]
      E_dict[U] = [lifted_numerator(pullback(projection_map[U])(I[i]))]
    end
    @show 4
    exc_div = IdealSheaf(covered_version, default_covering(covered_version), E_dict, check=false)
    @show "processing done."
    return projective_version, covered_version, projection_map, exc_div
  end
end

function blow_up(
    I::IdealSheaf;
    verbose::Bool=false,
    check::Bool=true,
  )
  X = scheme(I)
  C = covering(I)
  local_blowups = [blow_up(U, I[U], is_regular_sequence=is_regular_sequence(I), verbose=verbose, check=check) for U in all_patches(C) if is_defined_on(I, U)]
  ProjectivePatchType = projective_scheme_type(affine_patch_type(X))
  projective_glueings = Dict{Tuple{affine_patch_type(X), affine_patch_type(X)}, glueing_type(ProjectivePatchType)}()

  # prepare for the projective glueings
  for (U, V) in keys(glueings(C))
    P = local_blowups[C[U]][1]
    base_scheme(P) == U || error()
    SP = homog_poly_ring(P)
    Q = local_blowups[C[V]][1]
    base_scheme(Q) == V || error()
    SQ = homog_poly_ring(Q)
    G = C[U, V]
    UV, VU = glueing_domains(G)
    f, g = glueing_morphisms(G)

    # set up the maps on the patches of the overlap in U
    P_on_U_patches_to_Q = morphism_type(ProjectivePatchType)[]
    for i in 1:npatches(UV) 
      W = UV[i]
      gensUW = OO(W).(I[U])
      gensVW = OO(W).(pullback(f[i]).(I[V])) # TODO: Why is this conversion necessary?!
      transitionVU = [write_as_linear_combination(g, gensUW) for g in gensVW]
      PW, _ = fiber_product(W, P)
      SPW = homog_poly_ring(PW)
      push!(P_on_U_patches_to_Q, 
            ProjectiveSchemeMor(PW, Q, 
                                hom(SQ, SPW,
                                    pullback(f[i]),
                                    [dot(gens(SPW), SPW.(OO(W).(transitionVU[j]))) for j in 1:ngens(SQ)]
                                   )
                               )
           )
    end
    
    # set up the maps on the patches of the overlap in V
    Q_on_V_patches_to_P = morphism_type(ProjectivePatchType)[]
    for i in 1:npatches(VU) 
      W = VU[i]
      gensVW = OO(W).(I[V])
      gensUW = OO(W).(pullback(g[i]).(I[U]))
      transitionUV = [write_as_linear_combination(g, gensVW) for g in gensUW]
      QW, _ = fiber_product(W, Q)
      SQW = homog_poly_ring(QW)
      push!(Q_on_V_patches_to_P, 
            ProjectiveSchemeMor(QW, P, 
                                hom(SP, SQW,
                                    pullback(g[i]),
                                    [dot(gens(SQW), SQW.(OO(W).(transitionUV[j]))) for j in 1:ngens(SP)]
                                   )
                               )
           )
    end
    projective_glueings[(U, V)] = ProjectiveGlueing(G, P, Q, P_on_U_patches_to_Q, Q_on_V_patches_to_P)
  end
  tmp = Dict{affine_patch_type(X), ProjectivePatchType}()
  for U in patches(C)
    tmp[U] = local_blowups[C[U]][1]
  end
  projective_version = CoveredProjectiveScheme(X, C, tmp, projective_glueings)

  ### At this point we're done with the projective version of the blowup.
  # It remains to construct the associated CoveredScheme and the ideal sheaf 
  # of the exceptional divisor.
  covered_version = as_covered_scheme(projective_version)
  projection_map = covered_projection_to_base(projective_version)
  SpecType = affine_patch_type(covered_version)
  E_dict = Dict{SpecType, Vector{poly_type(SpecType)}}()
  for i in 1:length(local_blowups)
    merge!(E_dict, ideal_dict(local_blowups[i][4]))
  end
  exc_div = IdealSheaf(covered_version, default_covering(covered_version), E_dict)

  # Manually construct the ideal sheaf in each one of the charts of the blowup. 
  return projective_version, covered_version, projection_map, exc_div
end

function as_covered_scheme(P::CoveredProjectiveScheme)
  if !has_attribute(P, :as_covered_scheme) 
    X = base_scheme(P)
    C = base_covering(P)
    SpecType = affine_patch_type(X)
    new_patches = Vector{SpecType}()
    new_glueings = Dict{Tuple{SpecType, SpecType}, glueing_type(SpecType)}()
    projection_dict = Dict{SpecType, morphism_type(SpecType, SpecType)}() 
    for U in patches(C)
      PU = as_covered_scheme(P[U])
      new_patches = vcat(new_patches, patches(coverings(PU)[1]))
      merge!(new_glueings, glueings(PU))
      PU_projection = get_attribute(P[U], :covered_projection_to_base)
      merge!(projection_dict, morphisms(PU_projection))
    end
    #TODO: extend the remaining glueings
    new_covering = Covering(new_patches, new_glueings, check=false)
    covered_scheme = CoveredScheme(new_covering)
    covered_map = CoveringMorphism(new_covering, C, projection_dict) 
    projection_map = CoveredSchemeMorphism(covered_scheme, X, covered_map)
    set_attribute!(P, :as_covered_scheme, covered_scheme)
    set_attribute!(P, :covered_projection_to_base, projection_map)
  end
  return get_attribute(P, :as_covered_scheme)
end

function covered_projection_to_base(P::CoveredProjectiveScheme) 
  if !has_attribute(P, :covered_projection_to_base)
    as_covered_scheme(P)
  end
  return get_attribute(P, :covered_projection_to_base)
end



function strict_transform(f::SpecMor, h::Vector{PolyType}, g::Vector{PolyType}) where{PolyType<:MPolyElem}
	#(IOw: Exc Div ^\infty)
	
        X = domain(f)
        Y = codomain(f)
        R = base_ring(OO(X))
	Excdiv = ideal(h)

	Pf = pullback(f)
        Iold = ideal(R, lifted_numerator.(Pf.(g))) + strict_modulus(X)
	
	while true
          Inew = quotient(Iold, Excdiv)
          Iold == Inew && break
          Iold = Inew
	end
        return gens(Iold)
end



function total_transform(f::SpecMor, h::Vector{PolyType}, g::Vector{PolyType}) where{PolyType<:MPolyElem}
	#IOw
	
        X = domain(f)
        Y = codomain(f)
        R = base_ring(OO(X))
	Excdiv = ideal(h)

	Pf = pullback(f)
        Iold = ideal(R, lifted_numerator.(Pf.(g))) + strict_modulus(X)
        return gens(Iold)
end


### NOT TESTED YET
function weak_transform(f::SpecMor, h::Vector{PolyType}, g::Vector{PolyType}) where{PolyType<:MPolyElem}

	X = domain(f)
        Y = codomain(f)
        R = base_ring(OO(X))
	Excdiv = ideal(h)

	Pf = pullback(f)
        Iold = ideal(R, lifted_numerator.(Pf.(g))) + strict_modulus(X)

	while true
		Inew = quotient(Iold, Excdiv)
		!(Iold == Excdiv * Inew) && break
		Iold = Inew 
	end
	return gens(Iold)
	#(IOw : Exc Div ^k), k maximal
end

### NOT TESTED YET
function controlled_transform(f::SpecMor, h::Vector{PolyType}, g::Vector{PolyType}, i::Int) where{PolyType<:MPolyElem}
	#(IOw : Exc Div ^i)

	X = domain(f)
        Y = codomain(f)
        R = base_ring(OO(X))
	Excdiv = ideal(h)

	Pf = pullback(f)
        Iold = ideal(R, lifted_numerator.(Pf.(g))) + strict_modulus(X)
	

	for j in 1:i
		Inew = quotient(Iold,Excdiv)
		Inew == 1 && break
		Iold = Inew
	end
	
	return gens(Iold)

end

function strict_transform(
    P::ProjectiveScheme, 
    E::IdealSheaf, 
    I::Vector{PolyType}
  ) where {PolyType<:MPolyElem}
  X = as_covered_scheme(P)
  C = covering(E) 
  C in coverings(X) || error("covering not found")
  Y = base_scheme(P)
  length(I) > 0 || return IdealSheaf(X) # return the zero ideal sheaf
  for i in 1:length(I)
    parent(I[i]) == base_ring(OO(Y)) || error("polynomials do not belong to the correct ring")
  end
  f = covered_projection_to_base(P)
  if domain(f) !== C 
    f = compose(X[C, domain(f)], f)
  end

  SpecType = affine_patch_type(X)
  trans_dict = Dict{SpecType, Vector{poly_type(SpecType)}}()
  for U in patches(C)
    println("computing the strict transform on the patch $U")
    trans_dict[U] = strict_transform(f[U], E[U], I)
  end
  return IdealSheaf(X, C, trans_dict, check=false)
  #return IdealSheaf(X, CX, trans_dict, check=false)
end


function strict_transform(f::CoveredSchemeMorphism, E::IdealSheaf, I::IdealSheaf)
  X = domain(f)
  Y = codomain(f)
  CX = domain_covering(f)
  CY = codomain_covering(f)
  SpecType = affine_patch_type(X)
  trans_dict = Dict{SpecType, Vector{poly_type(SpecType)}}()
  for U in patches(CX)
    trans_dict[U] = strict_transform(f[U], E[U], I[codomain(f[U])])
  end
  return IdealSheaf(X, CX, trans_dict, check=true)
  #return IdealSheaf(X, CX, trans_dict, check=false)
end

function total_transform(f::CoveredSchemeMorphism, E::IdealSheaf, I::IdealSheaf)
  X = domain(f)
  Y = codomain(f)
  CX = domain_covering(f)
  CY = codomain_covering(f)
  SpecType = affine_patch_type(X)
  trans_dict = Dict{SpecType, Vector{poly_type(SpecType)}}()
  for U in patches(CX)
    trans_dict[U] = total_transform(f[U], E[U], I[codomain(f[U])])
  end
  return IdealSheaf(X, CX, trans_dict, check=true)
end

function prepare_smooth_center(I::IdealSheaf; check::Bool=true)
  X = scheme(I)
  C = covering(I)
  res_dict = Dict{affine_patch_type(X), Vector{poly_type(affine_patch_type(X))}}()
  center_dict = Dict{affine_patch_type(X), Vector{poly_type(affine_patch_type(X))}}()
  for U in patches(C)
    if check
      is_equidimensional_and_smooth(U) || error("ambient space is not smooth")
    end
    merge!(res_dict, prepare_smooth_center(U, I[U], check=check))
  end
  Z = subscheme(I)
  CZ = default_covering(Z)
  if check
    for U in patches(CZ)
      is_equidimensional_and_smooth(U) || error("center is not smooth")
    end
  end
  return I
end

function prepare_smooth_center(X::Spec, f::Vector{PolyType}; check::Bool=true) where {PolyType<:MPolyElem}
  X_dict = as_smooth_local_complete_intersection(X, verbose=true)
  return X_dict
end

function as_smooth_local_complete_intersection(I::IdealSheaf; check::Bool=true, verbose::Bool=false)
  X = scheme(I)
  C = covering(I)
  SpecType = affine_patch_type(X)
  PolyType = poly_type(SpecType)
  res_dict = Dict{SpecType, Tuple{Vector{PolyType}, Vector{PolyType}, Vector{Int}}}()
  for U in patches(C)
    merge!(res_dict, as_smooth_local_complete_intersection(U, I[U], check=check, verbose=verbose))
  end
  return res_dict
end

function as_smooth_local_complete_intersection(
    X::Spec, 
    f::Vector{PolyType}; 
    check::Bool=true, 
    verbose::Bool=false
  ) where {PolyType}
  verbose && println("call with $X and $f")
  R = base_ring(OO(X))
  g = gens(modulus(OO(X)))
  Dg = jacobi_matrix(g)
  d = dim(X)
  n = nvars(R)
  verbose && println("generators:")
  verbose && println(f)
  verbose && println("jacobian matrix:")
  if verbose
    for i in 1:nrows(Dg)
      println(Dg[i, 1:end])
    end
  end
  X_dict = _as_smooth_lci_rec(X, X, PolyType[], (Int[], Int[]), Vector{Tuple{Int, Int}}(), g, Dg, Dg, d, n, check=check, verbose=verbose)
  verbose && println("####################### done with the ambient space #######################")
  SpecType = typeof(X)
  f_dict = Dict{SpecType, Tuple{Vector{PolyType}, Vector{PolyType}, Vector{Int}}}()
  for U in keys(X_dict)
    verbose && println("proceed with finding a regular sequence on $U")
    g = X_dict[U][2]
    f_ext = vcat(g, f)
    Df_ext = jacobi_matrix(f_ext)
    good = X_dict[U][3]
    Z = subscheme(U, f)
    d = dim(Z)
    merge!(f_dict, _as_smooth_lci_rec(X, Z, X_dict[U][1], (collect(1:length(g)), good), 
                                      Vector{Tuple{Int, Int}}(), 
                                      f_ext, Df_ext, Df_ext, d, n, 
                                      check=check, verbose=verbose)
          )
  end
  return f_dict
end

function as_smooth_lci_of_cod(X::CoveredScheme, c::Int; check::Bool=true, verbose::Bool=false)
  C = default_covering(X)
  SpecType = affine_patch_type(X)
  PolyType = poly_type(SpecType)
  res_dict = Dict{SpecType, Tuple{Vector{PolyType}, Vector{PolyType}, Vector{Int}}}()
  for U in patches(C)
    merge!(res_dict, as_smooth_lci_of_cod(U, d, check=check, verbose=verbose))
  end
  return res_dict
end

function as_smooth_lci_of_cod(X::Spec, c::Int; check::Bool=true, verbose::Bool=false)
  R = base_ring(OO(X))
  f = gens(modulus(OO(X)))
  Df = jacobi_matrix(f)
  n = nvars(R)
#  verbose && println("generators:")
#  verbose && println(f)
#  verbose && println("jacobian matrix:")
#  if verbose
#    for i in 1:nrows(Df)
#      println(Df[i, 1:end])
#    end
#  end
#    PolyType = poly_type(X)
#    f_part = PolyType[]
#    for i in 1:c+3
#      push!(f_part, sum([(rand(Int)%1000)*g for g in f]))
#    end
#    @show "first try"
#    while !isempty(degeneracy_locus(X, jacobi_matrix(f_part), c, verbose=true))
#      @show "new try"
#      push!(f_part, dot([rand(Int)%1000 for j in 1:length(f)], f))
#    end
#    return _as_smooth_lci_rec(X, X, poly_type(X)[], (Int[], Int[]), Vector{Tuple{Int, Int}}(), f_part, jacobi_matrix(f_part), n-c, n, check=check, verbose=verbose)
  return _as_smooth_lci_rec(X, X, poly_type(X)[], (Int[], Int[]), Vector{Tuple{Int, Int}}(), f, Df, Df, n-c, n, check=check, verbose=verbose)
end

function as_smooth_local_complete_intersection(X::CoveredScheme; check::Bool=true, verbose::Bool=false)
  C = default_covering(X)
  SpecType = affine_patch_type(X)
  PolyType = poly_type(SpecType)
  res_dict = Dict{SpecType, Tuple{Vector{PolyType}, Vector{PolyType}, Vector{Int}}}()
  for U in patches(C)
    merge!(res_dict, as_smooth_local_complete_intersection(U, check=check, verbose=verbose))
  end
  return res_dict
end

function as_smooth_local_complete_intersection(X::Spec; check::Bool=true, verbose::Bool=false)
  R = base_ring(OO(X))
  f = gens(modulus(OO(X)))
  Df = jacobi_matrix(f)
  d = dim(X)
  n = nvars(R)
  verbose && println("generators:")
  verbose && println(f)
  verbose && println("jacobian matrix:")
  if verbose
    for i in 1:nrows(Df)
      println(Df[i, 1:end])
    end
  end
  return _as_smooth_lci_rec(X, X, poly_type(X)[], (Int[], Int[]), Vector{Tuple{Int, Int}}(), f, Df, Df, d, n, check=check, verbose=verbose)
end

### 
# Given a matrix A with polynomial entries, this routine returns a 
# list of equations hₖ and gₖ, and sets of columns Iₖ and rows Jₖ such that 
# the (Iₖ, Jₖ)-minor of A has the prescribed full rank r along the 
# hypersurface complement of hₖ of the zero locus of gₖ on C.
#
# Altogether, the latter sets form a disjoint union of C into locally 
# closed sets. In particular, we can derive an open covering of C 
# from the hypersurface complements of the listed (I, J)-minors.
#
# Input: 
#   C::Spec     the locus that needs to be covered
#   A::Matrix   the matrix that needs to have rank r
#   r::Int      the prescribed rank of the matrix.
#
# Output: 
#   the following lists with index i
#   loc_list::Vector{Vector{MPolyElem}} A list of partial
#               factorizations of polynomials that need to
#               be localized in order to arrive at the i-th 
#               patch.
#   row_list::Vector{Vector{Int}} the indices of the rows 
#               of the non-vanishing minor on the i-th patch.
#   column_list::Vector{Vector{Int}}    the indices of the 
#               columns of the non-vanishing minor on the 
#               i-th patch.
function _non_degeneration_cover(
    C::SpecType,
    A::MatrixType,
    r::Int; 
    rec_count::Int=0,
    check::Bool=true,
    verbose::Bool=false,
    restricted_rows::Vector{Vector{Int}}=[Int[]],
    restricted_columns::Vector{Vector{Int}}=[Int[]]
  ) where {SpecType<:Spec, MatrixType} #TODO: How specific can we be with the matrices?
  indent_str = prod([ "#" for i in 1:rec_count ]) * " "
  #verbose && println(indent_str * "call with $(X) and matrix $(A) for rank $r")
  verbose && println(indent_str * "call to _non_degeneration_cover for rank $r")
  verbose && println(indent_str * "on $C")
  verbose && println(indent_str * "for $A")
  verbose && println(indent_str * "with restricted rows $restricted_rows")
  verbose && println(indent_str * "and restricted columns $restricted_columns")

  # check for some other aborting conditions
  m = nrows(A)
  while length(restricted_rows) > 0 && length(restricted_rows[1]) == 0
    popfirst!(restricted_rows)
  end
  length(restricted_rows) == 0 && (restricted_rows = [collect(1:m)])
  n = ncols(A)
  while length(restricted_columns) > 0 && length(restricted_columns[1]) == 0
    popfirst!(restricted_columns)
  end
  length(restricted_columns) == 0 && (restricted_columns = [collect(1:n)])
  R = base_ring(OO(C))
  
  verbose && println("checking for the scheme being empty...")
  if isempty(C)
    verbose && println("the scheme is empty; returning empty lists")
    return (Vector{Vector{poly_type(C)}}(), Vector{Vector{poly_type(C)}}(), Vector{Vector{Int}}(), Vector{Vector{Int}}())
  end

#  if r == 1 
#    verbose && println("reached the rank-1-case")
#    Y = subscheme(C, ideal(R, [A[i,j] for i in 1:nrows(A) for j in 1:ncols(A)]))
#    if check 
#      isempty(Y) || error("the prescribed rank can not be achieved on $Y")
#    end
#    ll = Vector{Vector{poly_type(C)}}()
#    dl = Vector{Vector{poly_type(C)}}()
#    rl = Vector{Vector{Int}}() 
#    cl = Vector{Vector{Int}}()
#    for i in 1:nrows(A)
#      for j in 1:ncols(A)
#        if !iszero(OO(C)(A[i,j]))
#          push!(ll, [A[i,j]])
#          push!(dl, poly_type(C)[])
#          push!(rl, [i])
#          push!(cl, [j])
#        end
#      end
#    end
#    #verbose && println("returning localizations at $ll modulo $dl")
#    return ll, dl, rl, cl
#  end

  # find a suitable entry of lowest degree for the next recursion
  verbose && println("reducing the entries...")
  (k, l) = (1, 1)
  d = maximum(total_degree.(A[restricted_rows[1], restricted_columns[1]]))+1
  I = localized_modulus(OO(C))
  allzero = true
  W = localized_ring(OO(C))
  verbose && print(indent_str*" reducing row number ")
  for i in 1:m
    verbose && print("$i")
    v = W.([A[i,j] for j in 1:n])
    verbose && print(".")
    w = numerator.(reduce(v, groebner_basis(I))) # TODO: Implement and use reduction for matrices?
    verbose && print(".")
    for j in 1:n
      A[i,j] = w[j]
      if i in restricted_rows[1] && j in restricted_columns[1] 
        allzero = allzero && iszero(A[i,j])
        if total_degree(w[j]) <= d && !iszero(w[j])
          d = total_degree(w[j])
          (k, l) = (i, j)
        end
      end
    end
    verbose && print(".")
  end
  f = A[k,l]
  verbose && println("")

  # in case that after reduction all the matrix' entries are zero, quit
  if r> 0 && allzero
    error(indent_str * "All entries are zero")
  end

  verbose && println("oracle gives $(_degree_oracle_rec(total_degree.(A[restricted_rows[1], restricted_columns[1]]), rec_depth=r-1))")
  A_part = A[restricted_rows[1], restricted_columns[1]]
  deg_A_part = total_degree.(A[restricted_rows[1], restricted_columns[1]])
  (i,j), e = _degree_oracle_rec(deg_A_part, rec_depth=r-1)
  (k,l) = restricted_rows[1][i], restricted_columns[1][j]
  f = A[k,l]
  verbose && println(indent_str * "selected entry at ($k, $l): $f")

  verbose && println("preparing the matrix for induction...")
  B = copy(A)
  for i in 1:k-1
    multiply_row!(B, f, i)
    add_row!(B, -A[i,l], k, i)
  end
  for i in k+1:m
    multiply_row!(B, f, i)
    add_row!(B, -A[i,l], k, i)
  end

  # set up the submatrix for induction on the open part
  u = [i for i in 1:m if i != k]
  v = [j for j in 1:n if j != l]
  B_part = B[u, v]
  new_restricted_rows = [[(i > k ? i-1 : i) for i in a if i != k] for a in restricted_rows]
  new_restricted_columns = [[(j > l ? j-1 : j) for j in b if j != l] for b in restricted_columns]

  verbose && println("do the induction steps.")

  loc_list = Vector{Vector{poly_type(C)}}()
  div_list = Vector{Vector{poly_type(C)}}()
  row_list = Vector{Vector{Int}}() 
  column_list = Vector{Vector{Int}}()

  if r>1
    # prepare for recursion on the open part
    verbose && println("preparing the hypersurface complement")
    U = hypersurface_complement(C, f, keep_cache=true)

    # harvest from recursion on the open part
    verbose && println("recursive call for the complement")
    llU, dlU, rlU, clU = _non_degeneration_cover(U, B_part, r-1, 
                                                 rec_count=rec_count+1, 
                                                 check=check, 
                                                 verbose=verbose,
                                                 restricted_rows=new_restricted_rows,
                                                 restricted_columns=new_restricted_columns
                                                )

    # process the output according to the preparations for the recursion
    loc_list = [push!(Ul, f) for Ul in llU]
    div_list = dlU

    # the indices have to be adjusted according to the choice of the submatrix B from A
    row_list = [push!(Ul, k) for Ul in [[(i < k ? i : i+1) for i in a] for a in rlU]]
    column_list = [push!(Ul, l) for Ul in [[(i < l ? i : i+1) for i in a] for a in clU]]
    #verbose && println("return value was non-trivial; composing as $loc_list and $div_list")
    if check
      verbose && println("checking the return values")
      n = length(loc_list)
      for i in 1:n
        X = hypersurface_complement(subscheme(C, div_list[i]), prod(loc_list[i]))
        D = A[row_list[i], column_list[i]]
        g = det(D)
        isunit(OO(X)(g)) || error("selected minor is not a unit")
      end
    end
  else
    push!(loc_list, [f])
    push!(div_list, poly_type(C)[])
    push!(row_list, [k])
    push!(column_list, [l])
  end

  # prepare for recursion on the closed part
  verbose && println("preparing the hypersurface subscheme")
  Y = subscheme(C, f)

  # harvest from recursion on the closed part
  verbose && println("recursive call for the subscheme")
  llY, dlY, rlY, clY = _non_degeneration_cover(Y, copy(A), r, 
                                               rec_count=rec_count+1, 
                                               check=check, 
                                               verbose=verbose,
                                               restricted_rows=restricted_rows,
                                               restricted_columns=restricted_columns
                                              )
  
  # process the output according to the preparations for the recursion
  loc_list = vcat(loc_list, llY)
  div_list = vcat(div_list, [push!(Ud, f) for Ud in dlY])
  # no adjustment necessary, since this happens on the same recursion level
  row_list = vcat(row_list, rlY)
  column_list = vcat(column_list, clY)
  #verbose && println("return value was non-trivial; adding $llY and $dlY")

  return loc_list, div_list, row_list, column_list
end
    
function _as_smooth_lci_rec(
    X::SpecType, 
    Z::SpecType, # the uncovered locus
    h::Vector{PolyType}, # equations that were localized already
    good::Tuple{Vector{Int}, Vector{Int}}, # partial derivatives that form the beginning of the solution sequence
    bad::Vector{Tuple{Int, Int}}, # partial derivatives that vanish on the uncovered locus.
    f::Vector{PolyType}, 
    Df::MatrixType,
    B::MatrixType,
    d::Int, n::Int;
    check::Bool=true,
    verbose::Bool=false,
    rec_depth::Int=0
  ) where{
          SpecType<:Spec,
          PolyType<:MPolyElem,
          MatrixType
         }
  recstring = prod(["#" for i in 0:rec_depth])
  #verbose && println(recstring * "call with $X, $Z")
  verbose && println(recstring * "selected minors: $(good[1]) x $(good[2])")
  verbose && println(recstring * "bad positions: $bad")
  # return format: 
  # key: affine patch
  # value: (equations that were localized from root, 
  #         regular sequence, 
  #         index of variables for non-vanishing minor)
  res_dict = Dict{typeof(X), Tuple{Vector{PolyType}, Vector{PolyType}, Vector{Int}}}()

  if isempty(Z)
    verbose && println("this patch is already covered by others")
    return res_dict
  end

  if length(good[1]) == n-d
    verbose && println(recstring * "end of recursion reached")
    #g = det(Df[good[1], good[2]])
    #isempty(subscheme(Z, g)) || error("scheme is not smooth")
    res_dict[X] = (h, f[good[2]], good[1])
    if check
      issubset(localized_modulus(OO(X)), ideal(OO(X), f[good[2]])) || error("the found generators are not correct")
    end
    verbose && println(recstring*"returning $X with generators $(res_dict[X])")
    return res_dict
  end

  verbose && println("setting up matrix of minors")
  R = base_ring(OO(X))
  W = localized_ring(OO(X))
  J = localized_modulus(OO(Z))
  m = ncols(Df)
  n = nrows(Df)
  for i in [a for a in 1:n if !(a in good[1])]
    for j in [b for b in 1:m if !(b in good[2])]
      if !((i,j) in bad)
        #B[i,j] = numerator(reduce(W(B[i,j]), groebner_basis(localized_modulus(OO(X)))))
        #B[i,j] = numerator(reduce(W(det(Df[vcat(good[1], i), vcat(good[2], j)])), groebner_basis(localized_modulus(OO(X)))))
        B[i,j] = det(Df[vcat(good[1], i), vcat(good[2], j)])
      end
    end
  end
  min = maximum([total_degree(a) for a in collect(B)])
  (k, l) = (0, 0)
  verbose && println("reducing matrix and selecting pivot...")
  for i in [a for a in 1:n if !(a in good[1])]
    for j in [b for b in 1:m if !(b in good[2])]
      if !((i,j) in bad)
        #B[i,j] = numerator(reduce(W(B[i,j]), groebner_basis(localized_modulus(OO(X)))))
        #B[i,j] = numerator(reduce(W(det(Df[vcat(good[1], i), vcat(good[2], j)])), groebner_basis(localized_modulus(OO(X)))))
        #Df[i,j] = normal_form(Df[i,j], groebner_basis(saturated_ideal(localized_modulus(OO(X)))))
        if total_degree(B[i, j]) <= min && !iszero(OO(Z)(B[i, j]))
          min = total_degree(B[i, j])
          (k, l) = (i, j)
        end
      end
    end
  end
#  if verbose
#    for i in 1:nrows(Df)
#      println(Df[i, 1:end])
#    end
#  end
  if (k, l) == (0, 0)
    isempty(Z) && return res_dict
    error("scheme is not smooth")
  end
  verbose && println(recstring*"selected the $((k, l))-th entry: $(Df[k,l])")
  
  good_ext = (vcat(good[1], k), vcat(good[2], l))
  p = Df[k,l]
  Bnew = copy(B)
# for i in [a for a in 1:n if !(a in good_ext[1])]
#   if !iszero(Bnew[i,l])
#     @show i
#     @show n
#     multiply_row!(Bnew, p, i)
#     add_row!(Bnew, -B[i,l], k, i)
#   end
# end
#   for j in [b for b in 1:m if !(b in good_ext[2])]
#     if !iszero(Bnew[k,j])
#       @show j
#       multiply_column!(Bnew, p, j)
#       @show m
#       add_column!(Bnew, -Df[k,j], j, l)
#     end
#   end
  @show 1
  g = det(Df[good_ext[1], good_ext[2]])
  g = B[k,l]
  #g = B[k,l]
  @show total_degree(g)
  @show has_attribute(localized_modulus(OO(X)), :saturated_ideal)
  @show has_attribute(localized_modulus(OO(Z)), :saturated_ideal)
  I = localized_modulus(OO(Z))
  Isat = saturated_ideal(I)
  @show has_attribute(I, :saturated_ideal)
  @show 2
  U = hypersurface_complement(X, g, keep_cache=false)
  @show 2.5
  Z_next = hypersurface_complement(Z, g, keep_cache=true)
#  U = X
#  Z_next = Z
#  for a in factor(g)
#    @show a
#    U = hypersurface_complement(U, a[1])
#    Z_next = hypersurface_complement(Z_next, a[1])
#  end
  @show 3
  U_dict = _as_smooth_lci_rec(U, Z_next, vcat(h, [g]), good_ext, bad, f, Df, B, d, n, check=check, verbose=verbose, rec_depth=rec_depth+1)
  @show 4
  bad_ext = vcat(bad, (k, l))
  @show 5
  V_dict = _as_smooth_lci_rec(X, subscheme(Z, g), h, good, bad_ext, f, Df, B, d, n, check=check, verbose=verbose, rec_depth=rec_depth+1)
  @show 6
  return merge!(U_dict, V_dict)
end



function weak_transform(f::CoveredSchemeMorphism, E::IdealSheaf, I::IdealSheaf)
  X = domain(f)
  Y = codomain(f)
  CX = domain_covering(f)
  CY = codomain_covering(f)
  SpecType = affine_patch_type(X)
  trans_dict = Dict{SpecType, Vector{poly_type(SpecType)}}()
  for U in patches(CX)
    trans_dict[U] = weak_transform(f[U], E[U], I[codomain(f[U])])
  end
  return IdealSheaf(X, CX, trans_dict, check=true)
end

function controlled_transform(f::CoveredSchemeMorphism, E::IdealSheaf, I::IdealSheaf, i::Int)
  X = domain(f)
  Y = codomain(f)
  CX = domain_covering(f)
  CY = codomain_covering(f)
  SpecType = affine_patch_type(X)
  trans_dict = Dict{SpecType, Vector{poly_type(SpecType)}}()
  for U in patches(CX)
    trans_dict[U] = controlled_transform(f[U], E[U], I[codomain(f[U])],i)
  end
  return IdealSheaf(X, CX, trans_dict, check=true)
end

function test_cover(C, f, hl, ql, rl, cl)
  n = length(hl)
  Df = jacobi_matrix(f)
  for i in 1:n
    h = prod(hl[i])
    A = Df[rl[i], cl[i]]
    g = det(A)
    U = hypersurface_complement(subscheme(C, ql[i]), h)
    @show isunit(OO(U)(g))
  end


  I = ideal(OO(C), [det(Df[rl[i], cl[i]]) for i in 1:length(rl)])
  @show one(localized_ring(OO(C))) in I
end

function _degree_oracle_rec(M::Matrix{Int}; rec_depth::Int=0)
  m = nrows(M)
  n = ncols(M)
  d = maximum(M)
  minimal_entries = Vector{Tuple{Int, Int}}()
  for i in 1:m
    for j in 1:n
      if M[i,j] >= 0 && M[i,j] < d
        d = M[i,j]
        minimal_entries = [(i,j)]
      end
      if M[i,j] >=0 && M[i,j] <= d
        push!(minimal_entries, (i, j))
      end
    end
  end

  if rec_depth==0
    return minimal_entries[1], d
  end

  oracles = Vector{Tuple{Tuple{Int, Int}, Int}}()
  for p in 1:length(minimal_entries)
    (k,l) = minimal_entries[p]
    Mnew = Matrix{Int}(undef, m-1, n-1)
    for i in 1:k-1
      for j in 1:l-1
        Mnew[i,j] = maximum([(M[i,j] >= 0 ? M[i,j] + M[k,l] : -1), (M[i,l]>=0 ? M[k,j] + M[i,l] : -1)])
      end
      for j in l+1:n
        Mnew[i,j-1] = maximum([(M[i,j] >= 0 ? M[i,j] + M[k,l] : -1), (M[i,l]>=0 ? M[k,j] + M[i,l] : -1)])
      end
    end
    for i in k+1:m
      for j in 1:l-1
        Mnew[i-1,j] = maximum([(M[i,j] >= 0 ? M[i,j] + M[k,l] : -1), (M[i,l]>=0 ? M[k,j] + M[i,l] : -1)])
      end
      for j in l+1:n
        Mnew[i-1,j-1] = maximum([(M[i,j] >= 0 ? M[i,j] + M[k,l] : -1), 
                                 (M[i,l]>=0 ? M[k,j] + M[i,l] : -1)])
      end
    end
    push!(oracles, _degree_oracle_rec(Mnew))
  end

  e = minimum([a[2] for a in oracles])
  for p in 1:length(minimal_entries)
    if oracles[p][2] == e
      return minimal_entries[p], maximum([d, e])
    end
  end
  return minimal_entries[p], maximum([d, e])
end


