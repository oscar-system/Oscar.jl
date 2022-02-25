export ProjectiveGlueing
export glueing_type

export CoveredProjectiveScheme
export base_scheme, base_covering, projective_patches, as_covered_scheme

export blow_up, empty_covered_projective_scheme

export strict_transform, total_transform, weak_transform, controlled_transform

export prepare_smooth_center, as_smooth_local_complete_intersection

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
function blow_up(W::Spec, I::Vector{RingElemType}) where {RingElemType<:MPolyElem}
  @show W
  @show I

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
  Pw = projective_space(W,m-1)
  S = homogeneous_coordinate_ring(Pw)

  CP = affine_cone(Pw)
  Polyring = base_ring(OO(CP))
  At, embeddingAt, T =  _add_variables(R,[:t])
  t = T[1]

  #	@show vcat([t*embeddingAt(f) for f in I], gens(At)[1:end-1])
  Phi = AlgebraHomomorphism(Polyring, At, vcat([t*embeddingAt(f) for f in I], gens(At)[1:end-1]))
  kernel(Phi)

  Imod = modulus(A)
  IW = ideal(At, [embeddingAt(f) for f in gens(Imod)])
  IWpre = preimage(Phi, IW)
  #SIWpre = ideal(S,[frac_to_homog(Pw,g) for g in gens(IWpre)])
  SIWpre = ideal(S, poly_to_homog(Pw).(gens(IWpre)))

  @show "done. Proceed to process"
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
end

function blow_up(I::IdealSheaf)
  X = scheme(I)
  C = covering(I)
  local_blowups = [blow_up(U, I[U]) for U in patches(C)]
  ProjectivePatchType = projective_scheme_type(affine_patch_type(X))
  projective_glueings = Dict{Tuple{affine_patch_type(X), affine_patch_type(X)}, glueing_type(ProjectivePatchType)}()

  # prepare for the projective glueings
  for (U, V) in keys(glueings(C))
    P = local_blowups[C[U]][1]
    base_scheme(P) == U || error()
    SP = homogeneous_coordinate_ring(P)
    Q = local_blowups[C[V]][1]
    base_scheme(Q) == V || error()
    SQ = homogeneous_coordinate_ring(Q)
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
      SPW = homogeneous_coordinate_ring(PW)
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
      SQW = homogeneous_coordinate_ring(QW)
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
        Iold = ideal(R, lifted_numerator.(Pf.(g)))
	
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
        Iold = ideal(R, lifted_numerator.(Pf.(g)))
        return gens(Iold)
end


### NOT TESTED YET
function weak_transform(f::SpecMor, h::Vector{PolyType}, g::Vector{PolyType}) where{PolyType<:MPolyElem}

	X = domain(f)
        Y = codomain(f)
        R = base_ring(OO(X))
	Excdiv = ideal(h)

	Pf = pullback(f)
        Iold = ideal(R, lifted_numerator.(Pf.(g)))

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
        Iold = ideal(R, lifted_numerator.(Pf.(g)))
	

	for j in 1:i
		Inew = quotient(Iold,Excdiv)
		Inew == 1 && break
		Iold = Inew
	end
	
	return gens(Iold)

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

function as_smooth_local_complete_intersection(X::CoveredScheme; check::Bool=true, verbose::Bool=false)
  C = default_covering(X)
  SpecType = affine_patch_type(X)
  PolyType = poly_type(SpecType)
  res_dict = Dict{SpecType, Vector{PolyType}}()
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
  return _as_smooth_lci_rec(X, X, (Int[], Int[]), Vector{Tuple{Int, Int}}(), f, Df, d, n, check=check, verbose=verbose)
end

function _as_smooth_lci_rec(
    X::SpecType, 
    Z::SpecType, # the uncovered locus
    good::Tuple{Vector{Int}, Vector{Int}}, # partial derivatives that form the beginning of the solution sequence
    bad::Vector{Tuple{Int, Int}}, # partial derivatives that vanish on the uncovered locus.
    f::Vector{PolyType}, 
    Df::MatrixType,
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
  verbose && println(recstring * "call with $X, $Z")
  verbose && println(recstring * "selected minors: $(good[1]) x $(good[2])")
  verbose && println(recstring * "bad positions: $bad")
  res_dict = Dict{typeof(X), Vector{PolyType}}()

  if isempty(Z)
    verbose && println("this patch is already covered by others")
    return res_dict
  end

  if length(good[1]) == n-d
    verbose && println(recstring * "end of recursion reached")
    g = det(Df[good[1], good[2]])
    isempty(subscheme(Z, g)) || error("scheme is not smooth")
    U = hypersurface_complement(X, g)
    res_dict[U] = f[good[2]]
    if check
      issubset(localized_modulus(OO(U)), ideal(OO(U), f[good[2]])) || error("the found generators are not correct")
    end
    verbose && println(recstring*"returning $U with generators $(res_dict[U])")
    return res_dict
  end

  R = base_ring(OO(X))
  W = localized_ring(OO(X))
  J = localized_modulus(OO(Z))
  m = ncols(Df)
  n = nrows(Df)
  min = maximum([total_degree(a) for a in collect(Df)])
  (k, l) = (0, 0)
  for i in [a for a in 1:n if !(a in good[1])]
    for j in [b for b in 1:m if !(b in good[2])]
      if !((i,j) in bad)
        if !iszero(OO(X)(Df[i, j])) && total_degree(Df[i, j]) <= min
          min = total_degree(Df[i, j])
          (k, l) = (i, j)
        end
      end
    end
  end
  if verbose
    for i in 1:nrows(Df)
      println(Df[i, 1:end])
    end
  end
  if (k, l) == (0, 0)
    error("scheme is not smooth")
  end
  verbose && println(recstring*"selected the $((k, l))-th entry: $(Df[k,l])")

  good_ext = (vcat(good[1], k), vcat(good[2], l))
  g = det(Df[good_ext[1], good_ext[2]])
  U = hypersurface_complement(X, g)
  U_dict = _as_smooth_lci_rec(U, intersect(U, Z), good_ext, bad, f, Df, d, n, check=check, verbose=verbose, rec_depth=rec_depth+1)
  bad_ext = vcat(bad, (k, l))
  V_dict = _as_smooth_lci_rec(X, subscheme(Z, g), good, bad_ext, f, Df, d, n, check=check, verbose=verbose, rec_depth=rec_depth+1)
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
