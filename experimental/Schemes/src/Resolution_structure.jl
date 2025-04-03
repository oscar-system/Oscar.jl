export find_refinement_with_local_system_of_params
export inclusion_morphisms
export embedded_desingularization
export desingularization
export desingularization_only_blowups
export exceptional_locus
export NormalizationMorphism
export locus_of_maximal_order


##################################################################################################
# getters
##################################################################################################
underlying_morphism(phi::NormalizationMorphism) = phi.underlying_morphism
inclusion_morphisms(phi::NormalizationMorphism) = phi.inclusions
morphisms(phi::AbsDesingMor) = copy(phi.maps)
morphism(phi::AbsDesingMor,i::Int) = phi.maps[i]
last_map(phi::AbsDesingMor) = phi.maps[end]
normalization_steps(phi::MixedBlowUpSequence) = phi.normalization_steps

exceptional_divisor_list(phi::BlowUpSequence) = phi.ex_div
exceptional_divisor_as_ideal_sheafs(phi::MixedBlowUpSequence) = exceptional_divisor_list(phi,true)

function exceptional_divisor_as_ideal_sheafs(phi::BlowUpSequence)
  !is_empty(phi.ex_div) || return phi.ex_div
  E_list = phi.ex_div
  ret_list = Vector{AbsIdealSheaf}()
  for i in 1:length(E_list)
    E = (!(E_list[i] isa AbsIdealSheaf) ? ideal_sheaf(E_list[i]) : E_list[i])
    push!(ret_list, E)
  end
  return ret_list
end

## entries of ex_div corresponding to normalization steps are only exceptional divisors at the very end
## so only return them at the very end or on specific demand 
function exceptional_divisor_list(phi::MixedBlowUpSequence, seq_unclean::Bool=false)
  phi.resolves_sing || seq_unclean || error("exceptional divisor list not available for intermediate steps.")
  return phi.ex_div
end

 
embeddings(phi::BlowUpSequence) = phi.embeddings

function underlying_morphism(phi::AbsDesingMor)
  if !isdefined(phi, :underlying_morphism)
    phi.underlying_morphism = CompositeCoveredSchemeMorphism(reverse(morphisms(phi)))
  end
  return phi.underlying_morphism
end

@doc raw"""
    exceptional_divisor(f::AbsBlowUpSequence)   --> CartierDivisor

Return a `CartierDivisor` on the `domain` of `f` which is the
exceptional divisor of the sequence of blow-ups `f`.

# Examples
```jldoctest
julia> R,(x,y) = polynomial_ring(QQ,2);

julia> I=ideal(R,[x^2-y^5]);

julia> W = AffineScheme(R);

julia> IS = IdealSheaf(W,I);

julia> X = subscheme(IS);

julia> U = first(affine_charts(X));

julia> phi = desingularization_only_blowups(X);

julia> exceptional_divisor(phi)
Cartier divisor
  on scheme over QQ covered with 3 patches
with coefficients in integer ring
defined by the formal sum of
  1 * effective cartier divisor on scheme over QQ covered with 3 patches
```
"""
function exceptional_divisor(f::BlowUpSequence)
  f.is_embedded && return _exceptional_divisor_in_ambient(f)
  return _exceptional_divisor_non_embedded(f)
end

@doc raw"""
    exceptional_locus(f::AbsBlowUpSequence)

Return a `WeilDivisor` on the `domain` of `f` which is the
exceptional divisor of the sequence of blow-ups `f`.

# Examples
```jldoctest
julia> R,(x,y) = polynomial_ring(QQ,2);

julia> I=ideal(R,[x^2-y^5]);

julia> W = AffineScheme(R);

julia> IS = IdealSheaf(W,I);

julia> X = subscheme(IS);

julia> U = first(affine_charts(X));

julia> phi = desingularization_only_blowups(X);

julia> exceptional_locus(phi)
Effective weil divisor
  on scheme over QQ covered with 3 patches
with coefficients in integer ring
given as the formal sum of
  1 * sheaf of ideals

```
"""
function exceptional_locus(phi::BlowUpSequence)
  if !isdefined(phi, :exceptional_locus)
    phi.exceptional_locus = weil_divisor(exceptional_divisor(phi))
  end
  return phi.exceptional_locus
end

@doc raw"""
    exceptional_divisor_with_multiplicities(f::BlowUpSequence) --> CartierDivisor

Return a `CartierDivisor` `C` on the `domain` of the embedded desingularization morphism `f`
which is the exceptional divisor of the sequence of blow-ups `f` in the ambient scheme.
"""
function exceptional_divisor_with_multiplicities(f::BlowUpSequence)
  f.is_embedded || error("only available for embedded desingularization")
  f.transform_type != :strict || error("only available for weak and controlled transforms")

  ex_div_list = exceptional_divisor_list(f)
  C = sum(f.ex_mult[i] * cartier_divisor(ex_div_list[i]) for i in 1:length(ex_div_list))
  return C
end

function _exceptional_divisor_in_ambient(f::BlowUpSequence)
  f.is_embedded || error("only for embedded desingularization")
  !isdefined(f,:exceptional_divisor) || return f.exceptional_divisor
  f.exceptional_divisor = sum(f.ex_div; init=CartierDivisor(domain(f),ZZ))
  return f.exceptional_divisor
end

function _exceptional_divisor_non_embedded(f::MixedBlowUpSequence)
  !isdefined(f,:exceptional_divisor) || return f.exceptional_divisor

  ex_div_list = exceptional_divisor_list(f)
  C = WeilDivisor(scheme(ex_div_list[1]),ZZ)
  for i in 2:length(ex_div_list)
    dim(ex_div_list[i]) == -inf && continue            # kick out empty ones
    C = C + weil_divisor(ex_div_list[i])
  end

  f.exceptional_divisor = C
  return C
end

function _exceptional_divisor_non_embedded(f::BlowUpSequence)
  !isdefined(f,:exceptional_divisor_on_X) || return f.exceptional_divisor_on_X

  ex_div_list = exceptional_divisor_list(f)
  C = CartierDivisor(ambient_scheme(ex_div_list[1]),ZZ)
  for i in 1:length(ex_div_list)
# do we want to introduce is_empty for divisors?
    dim(ideal_sheaf(ex_div_list[i])) == -inf && continue          # kick out empty ones
                                                               # caution: dim(CartierDivisor) is not computed,
                                                               #          but inferred
                                                               #          ==> need to pass to ideal_sheaf first
    C = C + cartier_divisor(ex_div_list[i])
  end

  f.exceptional_divisor_on_X = C
  return C
end

function exceptional_divisor(f::MixedBlowUpSequence)
  !isdefined(f,:exceptional_divisor_on_X) || return f.exceptional_divisor_on_X
  f.resolves_sing || error("exceptional locus need not be a divisor for intermediate steps -- use exceptional_locus")

  return _exceptional_divisor_non_embedded(f)
end

function exceptional_locus(f::MixedBlowUpSequence)
  !isdefined(f,:exceptional_locus) || return f.exceptional_locus

  ex_div_list = exceptional_divisor_list(f)         # these are IdealSheaves for MixedBlowUpSequences
                                                    # they might even have the wrong dimension
  C = AlgebraicCycle(scheme(ex_div_list[1]),ZZ)     # ==> we cannot expect to obtain a divisor, only a cycle
  for E in ex_div_list
    dim(E) != -inf || continue                        # kick out empty ones
    C = C + algebraic_cycle(E)
  end

  return C
end

##################################################################################################
# setting values in DesingMors -- Watch out: only place with direct access to fields!!!
##################################################################################################
function add_map!(f::BlowUpSequence, phi::BlowupMorphism)
  f = update_dont_meet_pts!(f,center(phi))
  push!(f.maps, phi)
  ex_div = [strict_transform(phi,E) for E in f.ex_div[1:end]]
  push!(ex_div, exceptional_divisor(phi))
  f.ex_div = ex_div
  if isdefined(f, :underlying_morphism)
    f.underlying_morphism = CompositeCoveredSchemeMorphism(reverse(morphisms(f)))
  end  
  return f
end

function add_map!(f::MixedBlowUpSequence, phi::BlowupMorphism)
  f = update_dont_meet_pts!(f,center(phi))
  push!(f.maps, phi)
  ex_div = (AbsIdealSheaf)[strict_transform(phi,E) for E in f.ex_div]
  push!(ex_div, ideal_sheaf(exceptional_divisor(phi)))
  f.ex_div = ex_div
  if isdefined(f, :underlying_morphism)
    f.underlying_morphism = CompositeCoveredSchemeMorphism(reverse(morphisms(f)))
  end
  return f
end

function add_map!(f::MixedBlowUpSequence, phi::NormalizationMorphism)
  push!(f.maps, phi)
  sl = ideal_sheaf_of_singular_locus(codomain(phi))
  ex_div = (AbsIdealSheaf)[pullback(phi,E) for E in exceptional_divisorlist(f,true)]
  push!(ex_div,pullback(phi, sl))
  f.ex_div = ex_div
  push!(normalization_steps(f),length(f.maps))
  if isdefined(f, :underlying_morphism)
    f.underlying_morphism = CompositeCoveredSchemeMorphism(reverse(morphisms(f)))
  end
  return f
end

function add_map_embedded!(f::BlowUpSequence, phi::BlowupMorphism)
  push!(f.maps, phi)

  strict_list = [strict_transform_with_multiplicity(phi,E) for E in f.ex_div]
  ex_div = [a[1] for a in strict_list]
  push!(ex_div, exceptional_divisor(phi))
  f.ex_div = ex_div

  excess_mult = sum(a[2] for a in strict_list)
  if f.transform_type == :strict
    X_strict, inc_strict,_ = strict_transform(phi, f.embeddings[end])
    push!(f.embeddings, inc_strict)
  elseif f.transform_type == :weak
    I_trans,b = weak_transform_with_multiplicity(phi, f.controlled_transform)
    push!(f.ex_mult,b + excess_mult)
    f.controlled_transform = I_trans
    f.control = b
  else
    I_trans = controlled_transform(phi, f.controlled_transform, f.control)
    f.controlled_transform = I_trans
    push!(f.ex_mult, f.ex_mult[end])
  end
  if isdefined(f, :underlying_morphism)
    f.underlying_morphism = CompositeCoveredSchemeMorphism(reverse(morphisms(f)))
  end
  return f
end

function extend!(f::MixedBlowUpSequence, g::Union{BlowUpSequence,MixedBlowUpSequence})
  for g_i in morphisms(g)
    add_map!(f,g_i)
  end
  return f
end

function _blow_up_at_all_points_embedded(f::BlowUpSequence, I_all::AbsIdealSheaf)
  @assert domain(f) === scheme(I_all)

  decomp=maximal_associated_points(I_all)
  while !is_empty(decomp)
    I = small_generating_set(pop!(decomp))
    f = _do_blow_up_embedded!(f,I)
    if !is_empty(decomp)
      decomp = StrictTransformIdealSheaf[strict_transform(last_map(f),J) for J in decomp]
    end
  end
  return f
end

function _blow_up_at_all_points_embedded(f::BlowUpSequence, V::Vector{<:AbsIdealSheaf})
   !is_empty(V) || return(f)
   I = small_generating_set(pop!(V))
   f = _do_blow_up_embedded!(f,I)
   if length(V) > 0
     tempV = [strict_transform(last_map(f),J) for J in V]
     f = _blow_up_at_all_points_embedded(f,tempV)
   end
   return f
end

function _blow_up_at_all_points(f::Union{BlowUpSequence,MixedBlowUpSequence}, I_all::AbsIdealSheaf)
  @assert domain(f) === scheme(I_all)

  decomp=maximal_associated_points(I_all)
  while !is_empty(decomp)
    I = small_generating_set(pop!(decomp))
    f = _do_blow_up!(f,I)
    if length(decomp)>0
      decomp = [strict_transform(last_map(f),J) for J in decomp]
    end
  end
  return f
end

function _blow_up_at_all_points(f::Union{BlowUpSequence,MixedBlowUpSequence}, V::Vector{<:AbsIdealSheaf})
   !is_empty(V) || return(f)
   I = small_generating_set(pop!(V))
   f = _do_blow_up!(f,I)
   if length(V) > 0
     tempV = StrictTransformIdealSheaf[strict_transform(last_map(f),J) for J in V]
     f = _blow_up_at_all_points(f,tempV)
   end
   return f
end


function extend!(f::BlowUpSequence, g:: BlowUpSequence)
  for g_i in morphisms(g)
    add_map!(f,g_i)
  end
  return f
end

function update_dont_meet_pts!(f::Union{BlowUpSequence, MixedBlowUpSequence}, I::AbsIdealSheaf)
  dim(I) == 0 || return f
  is_prime(I) || return f
  patches = Oscar.patches(default_covering(scheme(I)))
  dont_meet = isdefined(f,:dont_meet) ? f.dont_meet : Vector{Tuple{Int,Int}}()
  i=1

  ## I describes a point, find a patch where it is visible
  i = findfirst(j ->!is_one(I(patches[j])), 1:length(patches))
  U = patches[i]

  ## now check containment for exceptional divisor
  C_list = exceptional_divisor_as_ideal_sheafs(f)      ## this looks clumsy, but we need the position in the list
  for i in 1:length(C_list)
    if !all(radical_membership(x,I(U)) for x in gens(C_list[i](U)))
      push!(dont_meet,(i,length(C_list)+1))            ## add the pair (position, new_one) to dont_meet
    end
  end

  f.dont_meet = dont_meet
  return f
end

function update_caution_multi_charts!(f::Union{BlowUpSequence, MixedBlowUpSequence},U :: AbsAffineScheme)
  dim(domain(phi)) == 2 || error("feature only available in the surface case")
  cent = center(phi)
  dim(cent) == 0 || return f

  ## find the charts originating from phi and prepare for using decomposition_info
  dom_covered = covered_scheme(domain(phi))
  patches_codom = patches(simplified_covering(codomain(phi)))
  chart_ind = findfirst(U -> !is_one(cent(U)), patches_codom)
  U = oatches_codom[chart_ind]
  U_blown_up = covered_scheme(P[U])
  has_decomposition_info(U_blown_up) || return f

  ## check whether a single chart suffices for the intersections not in dont_meet
  ex_div = exceptional_divisor_as_ideal_sheafs(f)
  cov_above_U = simplifed_covering(U_blown_up)
  patches_above_U = patches(cov_above_U)
  caution_multi_charts_new = Vector{Tuple{Int,Int}}
  E_last = ex_div[end]
  for i in 1:length(ex_div)-1
    findfirst(==(i,length(ex_div)), f.dont_meet) !== nothing || continue
    if length(findall(V -> !is_one(ex_div[i](V)+E_last(V)+decomposition_info(cov_above_U)[V]), patches_above_U)) != 1
      push!(caution_multi_charts_new, (i,length(ex_div)))
    end
  end
  append!(f.caution_multi_charts, caution_multi_charts_new)
  return f
end

function initialize_blow_up_sequence(phi::BlowupMorphism)
  f = BlowUpSequence([phi])
  f.ex_div = [exceptional_divisor(phi)]
  f.is_trivial = is_one(center(phi))
  f.is_strong = false
  f.resolves_sing = false                                # we have no information, whether we are done
                                                         # without further computation
  f.is_embedded = false
  return f
end

function initialize_mixed_blow_up_sequence(phi::NormalizationMorphism, I::AbsIdealSheaf)
  f = MixedBlowUpSequence([phi])
  f.ex_div = [pullback(phi,I)]
  f.is_trivial = is_one(I)
  f.is_strong = false
  f.resolves_sing = false                                # we have no information, whether we are done
                                                         # without further computation
  return f
end

function initialize_mixed_blow_up_sequence(phi::BlowupMorphism)
  return mixed_blow_up_sequence(initialize_blow_up_sequence(phi))
end

function initialize_embedded_blowup_sequence(phi::BlowupMorphism, inc::CoveredClosedEmbedding)
  f = BlowUpSequence([phi])
  f.ex_div = [exceptional_divisor(phi)]
  f.is_embedded = true
  f.transform_type = :strict
  if !is_one(center(phi))
    f.is_trivial = false
    X_strict,inc_strict,_ = strict_transform(phi,inc)
    f.embeddings = [f, inc_strict]
    f.is_strong = false                                  # we have no information, whether ncr will be ensured
    f.resolves_sing = false                              # we have no information, whether we are done
                                                         # without further computation
  else
    f.is_trivial = true
    f.embeddings = [inc, inc]
    f.is_strong = false                                  # should be set elsewhere
    f.resolves_sing = false                              # should be set elsewhere
  end
  return f
end

function initialize_embedded_blowup_sequence(phi::BlowupMorphism, I::AbsIdealSheaf, b::Int)
  f = BlowUpSequence([phi])
  f.ex_div = [exceptional_divisor(phi)]
  f.is_embedded = true
  if !is_one(center(phi))
    f.is_trivial = false
    if b == 0
      I_trans, b = weak_transform_with_multiplicity(phi,I)
      f.transform_type = :weak
    elseif b > 0
      I_trans = controlled_transform(phi, I, b)
      f.transform_type = :controlled
    end
    f.controlled_transform = I_trans
    f.ex_mult = [b]
    f.control = b
    f.is_strong = false
    f.resolves_sing = false                              # we have no information, whether we are done
                                                         # without further computation
  else
    f.is_trivial = true
    f.controlled_transform = I
    f.transform_type = :weak
    f.ex_mult = [0]
    f.is_strong = false
    f.resolves_sing = false                              # should be set elsewhere
  end
  return f
end

function forget_embedding(f::BlowUpSequence)
## create new BlowUpSequence, where maps contains maps between the domains of f
## - set is_embedded to false
## - inherit resolves_sing, is_trivial
## - new ex_div is induced divisor arising from f.ex_div on domain(f.maps[end])
## forget last blow-ups arising after the X has become smooth
  error("not implemented yet")
end

function mixed_blow_up_sequence(f::BlowUpSequence)
  phi = MixedBlowUpSequence(morphisms(f))
  phi.resolves_sing = f.resolves_sing
  phi.is_trivial = f.is_trivial
  phi.is_strong = f.is_strong
  phi.ex_div = [ideal_sheaf(E) for E in f.ex_div]
  phi.normalization_steps = Vector{Int}[]
  if isdefined(f, :underlying_morphism)
    phi.underlying_morphism = f.underlying_morphism
  end
  if isdefined(f, :exceptional_divisor)
    phi.exceptional_divisor = f.exceptional_divisor
  end
  return phi
end


##################################################################################################
# desingularization workers
##################################################################################################
function principalization(I::AbsIdealSheaf; algorithm::Symbol=:BEV)
  is_one(ideal_sheaf_of_singular_locus(scheme(I))) || error("principalization only available on smooth ambient schemes")

  ## trivial case: domain(f) was already smooth
  if is_one(I)
    id_W = identity_blow_up(scheme(I))
    phi = initialize_embedded_blowup_sequence(id_W,I,zero(ZZ))
    phi.resolves_sing = true
    phi.is_strong = true
    return phi
  end

  ## I non-empty, we need to do something
  dimX = dim(space(I))
  if dimX == 1
    return _principalize_curve(I)
#  elseif ((dimX == 2) && (algorithm == :CJS))
#    return _desing_CJS(f)
#  elseif (algorithm == :BEV)
#    return _desing_BEV(f)
  end
# here the keyword algorithm ensures that the desired method is called
  error("not implemented yet")
end

function embedded_desingularization(f::CoveredClosedEmbedding; algorithm::Symbol=:BEV)
  I_sl = ideal_sheaf_of_singular_locus(domain(f))

  ## trivial case: domain(f) was already smooth
  if is_one(I_sl)
    id_W = identity_blow_up(codomain(f))
    phi = initialize_embedded_blowup_sequence(id_W,f)
    phi.resolves_sing = true
    phi.is_strong = true
    return phi
  end

  ## I_sl non-empty, we need to do something
  dimX = dim(domain(f))
  if dimX == 1
    return _desing_emb_curve(f,I_sl)
#  elseif ((dimX == 2) && (algorithm == :CJS))
#    return _desing_CJS(f)
#  elseif (algorithm == :BEV)
#    return _desing_BEV(f)
  end
# here the keyword algorithm ensures that the desired method is called
  error("not implemented yet")
end

function embedded_desingularization(inc::ClosedEmbedding; algorithm::Symbol=:BEV)
  return embedded_desingularization(CoveredClosedEmbedding(inc); algorithm)
end

function CoveredClosedEmbedding(inc::ClosedEmbedding; domain=CoveredScheme(domain(inc)), codomain=CoveredScheme(codomain(inc)))
  mor_dict = IdDict{AbsAffineScheme, ClosedEmbedding}(domain[1][1] => inc)
  cov_mor = CoveringMorphism(default_covering(domain), default_covering(codomain), mor_dict; check=false)
  return CoveredClosedEmbedding(domain, codomain, cov_mor; check=false)
end

function desingularization(X::AbsCoveredScheme; algorithm::Symbol=:Lipman)
  I_sl = ideal_sheaf_of_singular_locus(X)
  
  ## trivial case: X is already smooth
  if is_one(I_sl)
    id_X = identity_blow_up(X)
    maps = [id_X] 
    return_value = BlowUpSequence(maps)
    return_value.resolves_sing = true
    return_value.is_trivial = true
    return_value.is_strong = true
    return mixed_blow_up_sequence(return_value)
  end

  Xnorm, phi = normalization(X)
  incs = phi.inclusions
  f = initialize_mixed_blow_up_sequence(phi,I_sl)
  I_sl = ideal_sheaf_of_singular_locus(Xnorm)

  ## I_sl non-empty, we need to do something 
# here the keyword algorithm ensures that the desired method is called
  dimX = dim(X)
  if dimX == 1
    return_value = f
    return_value.resolves_sing = true
    return_value.is_strong = true                # trivially, because disjoint points do not meet
  elseif ((dimX == 2) && (algorithm==:Lipman))
    return_value = _desing_lipman(Xnorm, I_sl, f)
  elseif ((dimX == 2) && (algorithm==:Jung))
    error("not implemented yet")
#    second_seq = _desing_jung(Xnorm,f)
#    return_value = extend!(phi, second_seq)
  else
    error("not implemented yet")
#    second_seq = forget_embedding(_desing_BEV(Xnorm))
#    return_value = extend!(phi, second_seq)
  end

  return return_value
end

function desingularization_only_blowups(X::AbsCoveredScheme; algorithm::Symbol=:Jung)
  I_sl = ideal_sheaf_of_singular_locus(X)

  ## trivial case: X is already smooth
  if is_one(I_sl)
    id_X = identity_blow_up(X)
    maps = [id_X]
    return_value = BlowUpSequence(maps)
    return_value.resolves_sing = true
    return_value.is_trivial = true
    return_value.is_strong = true
    return mixed_blow_up_sequence(return_value)
  end


  ## I_sl non-empty, we need to do something
  # here the keyword :algorithm ensures that the desired method is called
  dimX = dim(X)
  if dimX == 1
    return_value = _desing_curve(X, I_sl)
#  elseif ((dimX == 2) && (algorithm==:Jung))
#    error("not implemented yet")
#    return_value = _desing_jung(X)
  else
    error("not implemented yet")
    return_value = forget_embedding(_desing_BEV(X))
  end

  return return_value
end


function desingularization(X::AbsAffineScheme; algorithm::Symbol=:BEV)
  return desingularization(CoveredScheme(X); algorithm)
end

function weak_to_strong_desingularization_surface(phi::BlowUpSequence)
  return weak_to_strong_desingularization_surface(mixed_blow_up_sequence(phi))
end

function weak_to_strong_desingularization_surface(phi::MixedBlowUpSequence)
  dim(domain(phi)) == 2 || error("not implemented yet")

  !phi.is_strong || return phi
  ex_divs = phi.ex_div
  for i in normalization_steps(phi)
    !is_one(ex_divs[i]) || continue
    scheme_E,inc_E = sub(radical(ex_divs[i]))                             # here radical is cheap and necessary
    I_sl = ideal_sheaf_of_singular_locus(scheme_E)
    !is_one(I_sl) || continue
    psi = _desing_emb_curve(inc_E,I_sl)
    phi = extend!(phi,psi)
  end

  for i in 1:length(ex_divs)
    !(i in normalization_steps(phi)) || continue
    ex_divs[i] = radical(ex_divs[i])                                      # multiplicities here are artefacts
    I_sl, countA1s = curve_sing_A1_or_beyond(radical(ex_divs[i]))
    set_attribute!(phi.ex_div[i],:A1count,countA1s)
    !is_one(I_sl) || continue
    _,inc_E = sub(ex_divs[i])
    phi = _blow_up_at_all_points(phi,pushforward(inc_E)(I_sl))
  end

  return phi
end

function weak_to_strong_desingularization(phi::BlowUpSequence)
  error("not implemented yet")
end

function _principalize_curve(I::AbsIdealSheaf)

  dim(I) == 1 || error("IdealSheaf does not describe a curve.")

#  I_amb,I_bad,cod,b,MC = _nu_star_not_one(I)
#  has_dimension_leq_zero(I_bad) || error("IdealSheaf describes non-reduced subscheme")
#
#  decomp = maximal_associated_points(I_bad)                
#  cent = small_generating_set(pop!(decomp))
#  current_blow_up = blow_up(J)
#  initialize_embedded_blowup_sequence(current_blow_up,I,
#HIER HIER HIER
#  b==1 && return ... 

# ADD OTHER FRAGMENT HERE IN SEPARATE PR....
end

function _desing_curve(X::AbsCoveredScheme, I_sl::AbsIdealSheaf)
  ## note: I_sl not unit_ideal_sheaf, because this has been caught before in desingularization(X) 
  decomp = maximal_associated_points(I_sl)
  I = small_generating_set(pop!(decomp))
  current_blow_up = blow_up(I)
  phi = initialize_blow_up_sequence(current_blow_up)::BlowUpSequence
  phi = _blow_up_at_all_points(phi,decomp)
  
  I_sl_temp = ideal_sheaf_of_singular_locus(domain(phi))
  while !is_one(I_sl_temp)
    phi = _blow_up_at_all_points(phi,I_sl_temp)
    I_sl_temp = ideal_sheaf_of_singular_locus(domain(last_map(phi)))
  end

  phi.resolves_sing = true
  phi.is_strong = true                # trivially, as points do not meet
  return phi
end

function _desing_lipman(X::AbsCoveredScheme, I_sl::AbsIdealSheaf, f::MixedBlowUpSequence)
  dim(X) == 2 || error("Lipman's algorithm is not applicable")

  if dim(I_sl) == 1
  # not normal, do a normalization first
    f = _do_normalization!(f)
    Xnorm = domain(f)
    I_sl_temp = ideal_sheaf_of_singular_locus(Xnorm)
  else
  # already normal
    I_sl_temp = I_sl
  end    

  # now iterate this
  while !is_one(I_sl_temp)
    f = _blow_up_at_all_points(f,I_sl_temp)
    I_sl_temp = ideal_sheaf_of_singular_locus(domain(last_map(f)))
    if is_one(dim(I_sl_temp))
      f = _do_normalization!(f)
      I_sl_temp = ideal_sheaf_of_singular_locus(domain(last_map(f)))
    end
    decomp = maximal_associated_points(I_sl_temp)
  end

  f.resolves_sing = true
  f.is_strong = false        # not control of exceptional curves from normalization
  return f
end

function _desing_emb_curve(f::CoveredClosedEmbedding, I_sl::AbsIdealSheaf)
  ## note: I_sl not unit_ideal_sheaf, because this has been caught before in embedded_desingularization(f)
  decomp = maximal_associated_points(pushforward(f)(I_sl))
  I = small_generating_set(pop!(decomp))
  current_blow_up = blow_up(I)
  phi = initialize_embedded_blowup_sequence(current_blow_up,f)::BlowUpSequence
  decomp = [strict_transform(current_blow_up,J) for J in decomp]
  I_sl_temp = I_sl
  while !is_one(I_sl_temp)
    while !is_empty(decomp)
      I = small_generating_set(pop!(decomp))
      phi = _do_blow_up_embedded!(phi,I)
      if !is_empty(decomp)
        decomp = [strict_transform(last_map(phi),J) for J in decomp]
      end
    end
    last_emb = embeddings(phi)[end]
    I_sl_temp = pushforward(last_emb, ideal_sheaf_of_singular_locus(domain(last_emb)))
    decomp = maximal_associated_points(I_sl_temp)
  end

  phi = _ensure_ncr!(phi)
  phi.resolves_sing = true
  phi.is_strong = true
  return phi
end

function _ensure_ncr!(f::AbsDesingMor)
  current_divs = exceptional_divisor_list(f)

# this first step can only come into play, if the centers were not determined algorithmically
# it is ensured by all standard desingularization algorithms
  I_bad = non_snc_locus(current_divs)
  while !is_one(I_bad)
    f = _blow_up_at_all_points_embedded(f,I_bad)
    I_bad = non_snc_locus(exceptional_divisor_list(f))
  end

# now ensure the transversal intersections with the strict transform
  current_divs = copy(exceptional_divisor_list(f))
  I_X = image_ideal(f.embeddings[end])
  while !is_empty(current_divs)
    one_div = popfirst!(current_divs)       ## we need a FIFO here, not the usual LIFO
    I_temp = I_X + ideal_sheaf(one_div)
    last_emb = embeddings(f)[end]
    inc_temp = CoveredClosedEmbedding(scheme(I_temp), I_temp)
    next_locus = ideal_sheaf_of_singular_locus(domain(inc_temp))
    new_divs = length(f.ex_div) + 1
    f = _blow_up_at_all_points_embedded(f,pushforward(inc_temp, next_locus))
    append!(current_divs,[f.ex_div[i] for i in new_divs:length(f.ex_div)])
    I_X = image_ideal(f.embeddings[end])
  end

# finally make sure not too many exceptional divisors meet the strict transform in the same point
  n_max = dim(I_X)
  current_divs = copy(exceptional_divisor_list(f))
  _,inter_div = divisor_intersections_with_X(current_divs,I_X)
  f = _blow_up_at_all_points_embedded(f, inter_div)
  return f
end

function _do_blow_up!(f::AbsDesingMor, cent::AbsIdealSheaf)
  old_sequence = morphisms(f)
  X = domain(old_sequence[end])
  X === scheme(cent) || error("center needs to be defined on same scheme")
  current_blow_up = blow_up(cent,var_name=string("v", length(old_sequence), "_"))
  add_map!(f, current_blow_up)
  return(f)
end

function _do_normalization!(f::MixedBlowUpSequence)
  Xnorm, phi = normalization(domain(last_map(f)))
  add_map!(f,phi)
  return f
end
  
function _do_blow_up_embedded!(phi::AbsDesingMor,cent::AbsIdealSheaf)
  old_sequence = morphisms(phi)
  X = domain(old_sequence[end])
  X === scheme(cent) || error("center needs to be defined on same scheme")
  current_blow_up = blow_up(cent,var_name=string("v", length(old_sequence), "_"))
  add_map_embedded!(phi, current_blow_up)
  return(phi)
end


###################################################################################################
# Should go to IdealSheaf.jl, when PR is ready to merge
###################################################################################################

function unit_ideal_sheaf(X::AbsCoveredScheme)
  dd = IdDict{AbsAffineScheme, Ideal}(U=>ideal(OO(U), [one(OO(U))]) for U in affine_charts(X))
  return IdealSheaf(X, dd, check=false)
end

function zero_ideal_sheaf(X::AbsCoveredScheme)
  dd = IdDict{AbsAffineScheme, Ideal}(U=>ideal(OO(U), elem_type(OO(U))[]) for U in affine_charts(X))
  return IdealSheaf(X, dd, check=false)
end

function identity_blow_up(X::AbsCoveredScheme)
  f = BlowupMorphism(X, unit_ideal_sheaf(X))
  return f
end

########################################################################
# Refinements to find local systems of parameters
########################################################################

function find_refinement_with_local_system_of_params(W::AbsAffineScheme{<:Field, <:MPolyRing}; check::Bool=true)
  U = PrincipalOpenSubset(W, one(OO(W)))
  res_cov = Covering([U])
  R = ambient_coordinate_ring(W)
  minor_dict = IdDict{AbsAffineScheme, Tuple{Vector{Int}, Vector{Int}, elem_type(R)}}()
  minor_dict[U] = (Int[], Int[], one(R))
  return res_cov, minor_dict
end

function find_refinement_with_local_system_of_params(W::AbsAffineScheme; check::Bool=true)
  @check is_smooth(W) "scheme must be smooth"
  @check is_equidimensional(W) "scheme must be equidimensional"
  mod_gens = lifted_numerator.(gens(modulus(OO(W))))::Vector{<:MPolyRingElem}
  # We run into difficulties for the zero ideal as a modulus.
  # To get the matrix of minors of jac(I) we use `induced_map_on_exterior_power` below.
  # It is mathematical convention that ⋀⁰R⁰ = R¹. But that's unfortunately not 
  # coherent with the generic indexing we use in the implementation below. 
  # Should this assertion lead to problems, one can still replace throwing an error 
  # by inserting the shortcut from the method above and returning that. 
  @assert !isempty(mod_gens) "method not implemented for empty modulus; try to create the scheme without modulus instead"
  R = ambient_coordinate_ring(W)
  M = jacobian_matrix(R, mod_gens)

  n = nrows(M) # the number of variables in the ambient_ring
  r = ncols(M) # the number of generators

  Rn = FreeMod(R, n)
  Rr = FreeMod(R, r)
  phi = hom(Rn, Rr, M)
  codim = n - dim(W)
  phi_cod = induced_map_on_exterior_power(phi, codim)
  M_ext = matrix(phi_cod)
  n_cod = nrows(M_ext)
  r_cod = ncols(M_ext)

  all_entries = Vector{Int}[[i, j] for i in 1:n_cod for j in 1:r_cod]
  M_ext_vec = elem_type(R)[M[i, j] for i in 1:n_cod for j in 1:r_cod]
  min_id = ideal(OO(W), M_ext_vec)
  lambda_vec = coordinates(one(OO(W)), min_id)
  lambda = elem_type(OO(W))[lambda_vec[(i-1)*r_cod + j] for i in 1:n_cod, j in 1:r_cod]
  
  nonzero_indices_linear = [k for k in 1:length(lambda_vec) if !is_zero(lambda_vec[k])]
  non_zero_indices = [[i, j] for i in 1:n_cod, j in 1:r_cod if !is_zero(lambda_vec[(i-1)*r_cod + j])]

  ref_patches = AbsAffineScheme[]
  minor_dict = IdDict{AbsAffineScheme, Tuple{Vector{Int}, Vector{Int}, elem_type(R)}}()
  for (i, j) in non_zero_indices
    h_ij = M_ext[i, j]
    U_ij = PrincipalOpenSubset(W, OO(W)(h_ij))
    I = ordered_multi_index(i, codim, n)
    J = ordered_multi_index(j, codim, r)
    push!(ref_patches, U_ij)
    minor_dict[U_ij] = (indices(I), indices(J), M_ext[i, j])
  end
  res_cov = Covering(ref_patches)
  inherit_gluings!(res_cov, Covering(W))
  return res_cov, minor_dict
end

function find_refinement_with_local_system_of_params_rec(
    W::AbsAffineScheme, 
    mod_gens::Vector{PolyType} = lifted_numerator.(gens(modulus(OO(W)))),
    row_ind::Vector{Int} = Int[],
    col_ind::Vector{Int} = Int[],
    trans_mat::MatrixElem{RingElemType} = change_base_ring(OO(W), jacobi_matrix(mod_gens));
    check::Bool=true
  ) where {PolyType <: MPolyRingElem, RingElemType <: RingElem}

  # End of recursion
  n = dim(ambient_coordinate_ring(W))
  if length(row_ind) == n - dim(W)
    return [(W, row_ind, col_ind, prod(trans_mat[row_ind[k], col_ind[k]] for k in 1:dim(W); init=one(OO(W))))]
  end

  # generate the unit ideal of OO(W) with the entries of trans_mat
  n = nrows(trans_mat)
  r = ncols(trans_mat)
  all_entries_ind = [[i, j] for i in 1:n if !(i in row_ind) for j in 1:r if !(j in col_ind)]
  all_entries = elem_type(OO(W))[trans_mat[i, j] for (i, j) in all_entries_ind]
  entry_id = ideal(OO(W), all_entries)
  lambda = coordinates(one(OO(W)), entry_id)

  non_zero_entries = [k for k in 1:length(lambda) if !is_zero(lambda[k])]

  loc_results = Tuple{<:AbsAffineScheme, Vector{Int}, Vector{Int}, <:RingElem}[]
  for k in non_zero_entries
    i, j = all_entries_ind[k]
    h_ij = trans_mat[i, j]
    U_ij = PrincipalOpenSubset(W, OO(W)(h_ij))
    res_mat = change_base_ring(OO(U_ij), trans_mat) # TODO: Avoid checks here
    new_row_ind = vcat(row_ind, [i])
    new_col_ind = vcat(col_ind, [j])
    
    # Do Gaussian elimination on the matrix to kill off the other entries in this row
    u = res_mat[i, j]
    inv_u = inv(u)
    for l in 1:r
      l in new_col_ind && continue
      res_mat = add_column!(res_mat, -inv_u * res_mat[i, l], j, l)
    end
    loc_results = vcat(loc_results, 
          find_refinement_with_local_system_of_params_rec(
              U_ij, mod_gens, new_row_ind, new_col_ind, res_mat; check=check
             )
         )
  end
  return loc_results
end

##################################################################################################
#  locus of order at least b and of maximal order
##################################################################################################

function locus_of_maximal_order(I::AbsIdealSheaf)
  delta_list = _delta_list(I)
  return (delta_list[end],length(delta_list))
end

function locus_of_order_geq_b(I::AbsIdealSheaf, b::Int)
  return _delta_list(I,b)[end]
end

function _delta_ideal_for_order(inc::CoveredClosedEmbedding)

  W = codomain(inc)
  I = image_ideal(inc)

  Delta_dict = IdDict{AbsAffineScheme,Ideal}()

  for U in affine_charts(W)
    _, inc_U = sub(U,I(U)) 
    Cov,Chart_dict = find_refinement_with_local_system_of_params(U)
    Delta_dict[U] = (_delta_ideal_for_order(CoveredClosedEmbedding(inc_U), Cov, Chart_dict))(U)
  end

  return IdealSheaf(W,Delta_dict;check=false)
end

function _delta_list(inc::CoveredClosedEmbedding)
  I = image_ideal(inc)
  return _delta_list(I)
end

function _delta_list(I::AbsIdealSheaf, b::Int=0)
  W = scheme(I)
  is_smooth(W) || error("ambient scheme needs to be smooth")
  Delta_list = AbsIdealSheaf[]
  j = 0
  while (( !is_one(I) && b == 0 ) || ( j < b ))
    push!(Delta_list, I)
    inc = CoveredClosedEmbedding(scheme(I),I)
    I = _delta_ideal_for_order(inc)
    j = j + 1
  end

  return Delta_list
end

function _delta_ideal_for_order(inc::CoveredClosedEmbedding, Cov::Covering, 
       ambient_param_data::IdDict{<:AbsAffineScheme,
                                  <:Tuple{Vector{Int64},Vector{Int64},<:RingElem}};
       check::Bool=true)

  W = codomain(inc)                                
  @check is_smooth(W) "codomain of embedding needs to be smooth"
#  @check is_equidimensional(W) "codomain of embedding needs to be equidimensional"
  I_X = small_generating_set(image_ideal(inc))         # ideal sheaf describing X on W

  ## run through all charts, compute the jacobian matrices w.r.t. the respective system of parameters
  ## and flatten the matrix to append it to gens(I)
  Delta_dict = IdDict{AbsAffineScheme,Ideal}()
  for U in Cov
    I = I_X(U)
    if is_one(I)
      Delta_dict[U] = I
      continue
    end

    result_mat = _jacobian_matrix_wrt_system_of_param(U,I,ambient_param_data[U];check=false)
    Ivec = copy(gens(I))
    append!(Ivec,[a for a in OO(U).(collect(result_mat))])
    Delta_dict[U] = ideal(OO(U),Ivec)
  end

  return small_generating_set(IdealSheaf(W,Delta_dict))
end

##############################################################################################
#  jacobian matrix with respect to a consistent choice of local system of parameters on a    #
#  sufficiently small chart                                                                  #
##############################################################################################
function _jacobian_matrix_wrt_system_of_param(W::AbsAffineScheme{<:Any, <:MPolyRing}, I::Ideal,
                    ambient_param_data::Tuple{Vector{Int64},Vector{Int64},<:RingElem};
       check::Bool=true)
  return jacobian_matrix(OO(W),gens(I))
end

function  _jacobian_matrix_wrt_system_of_param(W::AbsAffineScheme{<:Any, <:MPolyLocRing}, I::Ideal,
                    ambient_param_data::Tuple{Vector{Int64},Vector{Int64},<:RingElem};
       check::Bool=true)
  I_gens = lifted_numerator.(gens(I))
  return jacobian_matrix(base_ring(OO(W)),I_gens)
end

function _jacobian_matrix_wrt_system_of_param(W::AbsAffineScheme,I::Ideal,
                    ambient_param_data::Tuple{Vector{Int64},Vector{Int64},<:RingElem};
       check::Bool=true)

  @check is_smooth(W) "ambient space W needs to be smooth"
  @check base_ring(I) == OO(W) "ideal not defined in the correct ring"

  ## non-trivial ambient scheme W
  R = base_ring(OO(W))
  I_X = small_generating_set(ideal(R,lifted_numerator.(gens(I))))
                                                     # ideal sheaf describing X on W by ideal in OOP
  mod_gens = lifted_numerator.(gens(modulus(OO(W)))) # generators of ideal of W
  amb_row,amb_col,h = ambient_param_data             # data of specific invertible minor

  ## set up pseudoinverse of jacobian matrix of generators of W
  JM = jacobian_matrix(R, mod_gens)
  if length(amb_col) < length(mod_gens)
      JM_essential = JM[:, amb_col]
  else
      JM_essential = JM
  end
  submat_for_minor = JM[amb_row, amb_col]
  Ainv, h2 = pseudo_inv(submat_for_minor)

  ## sanity check make sure that we are inverting the right minor
  h == h2 || error("inconsistent input data")
  JM_essential = JM_essential * Ainv                 # appropriate submatrix transformed to h*unit_matrix
                                                     ## TODO: should be cached!!

  ## now do Gaussian reduction on I_gens w.r.t. JM_essential
  I_gens = lifted_numerator.(gens(I))
  JI = jacobian_matrix(I_gens)
  result_mat = h*JI                                  # h is a unit (we only have h*unit_matrix, not unit_matrix)

  for i in 1:length(amb_col)
    for j in 1:ncols(result_mat)
      result_mat[:,j] = result_mat[:,j] - [JM_essential[k,i] * JI[amb_row[i],j] for k in 1:nrows(JM_essential)]
    end
  end

  return result_mat
end

########################################################################
# test for snc                                                         #
########################################################################

# check for singularities worse than A1 (they lie at ret_ideal_sheaf)
# and count the A1 encountered on the way (their count is total_number)
function curve_sing_A1_or_beyond(I::AbsIdealSheaf)
  !is_one(I) || return(I,0)
  @assert is_one(dim(I))
  I_scheme,I_inc = sub(I)
  I_sl = pushforward(I_inc)(ideal_sheaf_of_singular_locus(I_scheme))
  decomp = maximal_associated_points(I_sl)    # zero-dimensional, as I describes a curve
  init_ideal_sheaf = unit_ideal_sheaf(scheme(I_sl))
  ret_ideal_sheaf = init_ideal_sheaf
  total_number = 0
  for J in decomp
    check_val,local_number = is_A1_at_point_curve(I,J)
    if check_val
      total_number = total_number + local_number
    else
      ret_ideal_sheaf = (ret_ideal_sheaf != init_ideal_sheaf ? ret_ideal_sheaf * J : J)
    end
  end
  return ret_ideal_sheaf, total_number
end

function is_A1_at_point_curve(IX::AbsIdealSheaf,Ipt::AbsIdealSheaf)
  @assert scheme(IX) === scheme(Ipt)
  @assert dim(scheme(IX)) == 2
  @assert dim(Ipt) == 0

  patches_scheme = patches(simplified_covering(scheme(Ipt)))
  found_index = findfirst(V -> !is_one(Ipt(V)), patches_scheme)
  U = patches_scheme[found_index]
  IXsat = saturated_ideal(IX(U))
  R = base_ring(IXsat)
  IXU = ideal(R,small_generating_set(IXsat))
  IptU = saturated_ideal(Ipt(U))

  # absolutely irreducible point
  if vector_space_dimension(quo(R,IptU)[1]) == 1
    return (check_A1_at_point_curve(IXU,IptU) ? (true,1) : (false,0))
  else
    decomp = absolute_primary_decomposition(IptU)
    length(decomp) == 1 || error("decomposition problem of point")
    mult = decomp[1][4]
    I_kbar = decomp[1][3]
    r_changed = base_ring(I_kbar)
    kk = coefficient_ring(r_changed)
    IXU_changed = ideal(r_changed,  [change_coefficient_ring(kk,a, parent = r_changed) for a=gens(IXU)])
    IptU_changed = ideal(r_changed,  [change_coefficient_ring(kk,a, parent = r_changed) for a=gens(IptU)])
    return (check_A1_at_point_curve(IXU_changed,IptU_changed) ? (true, mult) : (false,0))
  end
end

function check_A1_at_point_curve(IX::Ideal, Ipt::Ideal)
## only call this from higher functions
## it assumes: IX singular at Ipt, germ of IX at Ipt contact equivalent to hypersurface singularity
  R = base_ring(IX)
  is_one(dim(IX)) || error("not applicable: not a curve")
  R == base_ring(Ipt) || error("basering mismatch")
  kk = base_ring(R)
  characteristic(kk) == 0 || error("only available in characteristic zero")

  JM = jacobian_matrix(gens(IX))

  ## localize at point
  a = rational_point_coordinates(Ipt)
  U = complement_of_point_ideal(R,a)
  RL, loc_map = localization(R,U)
  IX_loc = loc_map(IX)
  JM_loc =  map(loc_map, JM[:,:])

  if !all(iszero(a))
    F_loc = free_module(RL,ngens(IX))
    Jm = sub(F_loc,JM_loc)[1]
    Jm = Jm + (loc_map(IX)*F_loc)[1]
    Jm_shifted = shifted_module(Jm)[1]
    F_shifted = ambient_free_module(Jm_shifted)
  else
    F = free_module(R,ngens(IX))
    Jm = sub(F,JM)[1]
    Jm_shifted = Jm + (IX * F)[1]
    F_shifted = F
  end

  o = negdegrevlex(R)*lex(F_shifted)
  F1 = leading_module(Jm_shifted,o)
  F1quo = quo(F_shifted, F1)[1]

  return vector_space_dimension(F1quo) == 1
end

function divisor_intersections_with_X(current_div::Vector{<:EffectiveCartierDivisor}, I_X::AbsIdealSheaf)
  return divisor_intersections_with_X(ideal_sheaf.(current_div), I_X)
end
  
function divisor_intersections_with_X(current_div::Vector{<:AbsIdealSheaf}, I_X::AbsIdealSheaf)
  scheme(I_X) == scheme(current_div[1]) || error("underlying schemes do not match")
  n_max = dim(I_X)

  inter_div_dict = Dict{Vector{Int},Tuple{AbsIdealSheaf,Int}}()
  old_keys = Vector{Int}[]
  empty_keys = Vector{Int}[]
  essential_inter = AbsIdealSheaf[]

# initialization: each divisor + I_X
  for k in 1:length(current_div)
    Idiv = current_div[k] + I_X
    if !is_one(Idiv)
      inter_div_dict[[k]] = (Idiv,0)
      push!(old_keys, [k])
    end
  end
  new_keys = copy(empty_keys)

# add intersections
  while !is_empty(old_keys)
    old_keyvec = popfirst!(old_keys)       # this is the intersection to which we add a new divisor
    for i in (old_keyvec[end]+1):length(current_div)
      haskey(inter_div_dict, [i]) || continue     # divisor i meets X at all
      copykey = copy(old_keyvec)
      push!(copykey,i)                     # here we add it
      Idiv = inter_div_dict[[i]][1] + inter_div_dict[old_keyvec][1]
      !is_one(Idiv) || continue            # it intersects

      subsetlist = subsets(old_keyvec,length(old_keyvec)-1)
      subsetlist = [push!(a,i) for a in subsetlist]
      if (sum([inter_div_dict[a][2] for a in subsetlist])> 0)
                                           # offending intersection, known before
        inter_div_dict[copykey] = (Idiv,2)
      elseif dim(Idiv) == n_max - length(copykey)
                                           # non-offending intersection
        inter_div_dict[copykey] = (Idiv,0)
      else
                                           # offending intersection, new
        inter_div_dict[copykey] = (Idiv,1)
        push!(essential_inter, Idiv)
      end
      push!(new_keys,copykey)
    end

    # go to next higher number of intersecting divisors
    if is_empty(old_keys)
      old_keys = copy(new_keys)
      new_keys = copy(empty_keys)
    end
  end

  return inter_div_dict, essential_inter
end

function non_snc_locus(divs::Vector{<:EffectiveCartierDivisor})
  is_empty(divs) && error("list of divisors must not be empty")
  X = ambient_scheme(first(divs))
  @assert all(d->ambient_scheme(d) === X, divs)
  @assert is_smooth(X)
  r = length(divs)
  triv_cov = trivializing_covering.(divs)

  com_ref, incs = common_refinement(triv_cov, default_covering(X))

  ideal_dict = IdDict{AbsAffineScheme, Ideal}() # ideal sheaf of the non_snc_locus
  for U in patches(com_ref)
    loc_eqns = elem_type(OO(U))[]
    for k in 1:length(incs)
      I = ideal_sheaf(divs[k])
      inc = incs[k][U]
      V = codomain(inc)
      h = gen(I(V),1)
      hh = pullback(inc)(h)
      is_unit(hh) && continue # Not every divisor needs to be visible here
      push!(loc_eqns, hh)
    end
    if isempty(loc_eqns) || is_regular_sequence(loc_eqns)
      ideal_dict[U] = ideal(OO(U), one(OO(U))) # Nothing visible here
      continue
    end

    # Determine the non-snc-locus
    K = koszul_complex(KoszulComplex, loc_eqns)
    k = findfirst(k->!is_zero(homology(K, k)[1]), 1:length(loc_eqns))
    ideal_dict[U] = annihilator(homology(K, k)[1])
  end
  return IdealSheaf(X, ideal_dict; check=false)
end

function common_refinement(list::Vector{<:Covering}, def_cov::Covering)
  isempty(list) && error("list of coverings must not be empty")

  if length(list) == 1
    result = first(list)
    return result, [identity_map(result)]
  end
  patch_list = AbsAffineScheme[]
  anc_list = AbsAffineScheme[]
  to_U_dict = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
  to_V_dict = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()

  if length(list) == 2
    for U in patches(list[1])
      match_found = false
      for V in patches(list[2])
        success, W = _have_common_ancestor(U, V)
        !success && continue
        match_found = true
        push!(anc_list, W)
        #inc_U = _flatten_open_subscheme(U, W)
        #inc_V = _flatten_open_subscheme(V, W)
        inc_U, h_U = _find_chart(U, W)
        inc_U = PrincipalOpenEmbedding(inc_U, h_U; check=false)
        inc_V, h_V = _find_chart(V, W)
        inc_V = PrincipalOpenEmbedding(inc_V, h_V; check=false)

        UV, to_U, to_V = fiber_product(inc_U, inc_V) 
        push!(patch_list, UV)
        to_U_dict[UV] = to_U
        to_V_dict[UV] = to_V
      end
      !match_found && error("no common ancestor found for $U and $V")
    end
    #anc_cov = Covering(anc_list)
    #inherit_gluings!(anc_cov, def_cov)
    result = Covering(patch_list)
    inherit_gluings!(result, def_cov)

    tot_inc1 = CoveringMorphism(result, list[1], to_U_dict; check=false)
    tot_inc2 = CoveringMorphism(result, list[2], to_V_dict; check=false)
    return result, [tot_inc1, tot_inc2]
  end

  # More than two entries
  n = length(list)
  k = div(n, 2)
  res1, inc1 = common_refinement(list[1:k], def_cov)
  res2, inc2 = common_refinement(list[k+1:end], def_cov) 

  result, inc_tot = common_refinement([res1, res2], def_cov)
  return result, vcat([compose(inc_tot[1], inc1[k]) for k in 1:length(inc1)], 
                      [compose(inc_tot[2], inc2[k]) for k in 1:length(inc2)]
                     )
end

function strict_transform(bl::AbsSimpleBlowupMorphism, inc::ClosedEmbedding)
  B = codomain(bl)
  @assert length(affine_charts(B)) == 1 && first(affine_charts(B)) === codomain(inc)
  inc_cov = CoveredClosedEmbedding(inc, codomain=B)
  return strict_transform(bl, inc_cov)
end

is_graded(R::Ring) = false

########################################################################
# Implement the interface specified in 
# experimental/Schemes/BlowupMorphism.jl
########################################################################

# The following two methods will be ambiguous in general
# so we need to repeat the procedure for the specific types 
# of the second argument.
function strict_transform(phi::Union{BlowUpSequence,MixedBlowUpSequence}, a::Any)
  b = a
  for psi in morphisms(phi)
    b = strict_transform(psi, b)
  end
  return b
end

function total_transform(phi::Union{BlowUpSequence,MixedBlowUpSequence}, a::Any)
  for psi in morphisms(phi)
    a = total_transform(psi, a)
  end
  return a
end

function total_transform(phi::NormalizationMorphism, a::Any)
  return pullback(phi,a)
end

function strict_transform(phi::NormalizationMorphism, a::Any)
  return pullback(phi,a)
end

function strict_transform(phi::BlowUpSequence, a::AbsIdealSheaf)
  for psi in morphisms(phi)
    a = strict_transform(psi, a)
  end
  return a
end

function total_transform(phi::BlowUpSequence, a::AbsIdealSheaf)
  for psi in morphisms(phi)
    a = total_transform(psi, a)
  end
  return a
end

function last_center(phi::BlowUpSequence)
  return center(last_map(phi))
end

function last_center(phi::MixedBlowUpSequence)
  lm = last_map(phi)
  lm isa BlowupMorphism || return unit_ideal_sheaf(codomain(lm))
  return center(lm)
end
