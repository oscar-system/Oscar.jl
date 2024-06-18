export find_refinement_with_local_system_of_params
export inclusion_morphisms
export embedded_desingularization
export desingularization
export desingularization_only_blowups
export exceptional_locus
export NormalizationMorphism
export locus_of_maximal_order

##############################################################################
## Concrete Type for normalization
## very similar to CoveredSchemeMorphism, but allowing disjoint handling
## of disjoint components
##############################################################################
@doc raw"""
    NormalizationMorphism{
                  DomainType<:AbsCoveredScheme,
                  CodomainType<:AbsCoveredScheme
      } <:AbsCoveredSchemeMorphism{
                                 DomainType,
                                 CodomainType,
                                 Nothing,
                                 NormalizationMorphism,
                                }
A datastructure to encode normalizations of covered schemes.

It is described as the morphism from the new scheme to the original one, containing
information on the decomposition of the new scheme into disjoint components.
(This is the type of the return value of  `normalization(X::AbsCoveredScheme)`.)
"""
@attributes mutable struct NormalizationMorphism{
    DomainType<:AbsCoveredScheme,
    CodomainType<:AbsCoveredScheme
   } <:AbsCoveredSchemeMorphism{
                                 DomainType,
                                 CodomainType,
                                 Nothing,
                                 NormalizationMorphism,
                                }

  underlying_morphism::CoveredSchemeMorphism
  inclusions::Vector{<:AbsCoveredSchemeMorphism}

  function NormalizationMorphism(
      f::CoveredSchemeMorphism,
      inclusions::Vector{<:AbsCoveredSchemeMorphism};
      check::Bool=true
    ) 
    @check is_normal(X) "not a normalization morphism"
    @assert all(inc->codomain(inc) === domain(f), inclusions) "domains and codomains do not match"
    ret_value = new{typeof(domain(f)),typeof(codomain(f))}(f,inclusions)
    return ret_value    
  end

  function NormalizationMorphism(
      f::CoveredSchemeMorphism;
      check::Bool=true
    )
    @check is_normal(X) "not a normalization morphism"
    ret_value = new{typeof(domain(f)),typeof(codomain(f))}(f,[identity_map(X)])
    return ret_value    
  end

end

#####################################################################################################
# Desingularization morphism: birational map between covered schemes with smooth domain
#####################################################################################################
@doc raw"""
    MixedBlowUpSequence{
              DomainType<:AbsCoveredScheme,
              CodomainType<:AbsCoveredScheme
            }<:AbsDesingMor{  DomainType,
                              CodomainType, 
                              MixedBlowUpSequence{DomainType, CodomainType}
                                              }
A datastructure to encode sequences of blow-ups and normalizations of covered schemes
as needed for desingularization of non-embedded schemes by the approaches of Zariski and of
Lipman. 
"""

@attributes mutable struct MixedBlowUpSequence{
                                                DomainType<:AbsCoveredScheme,
                                                CodomainType<:AbsCoveredScheme
                                              }<:AbsDesingMor{  DomainType,
                                                                CodomainType, 
                                                                MixedBlowUpSequence{DomainType, CodomainType}
                                              }
  maps::Vector{Union{<:BlowupMorphism,<:NormalizationMorphism}}       # count right to left:
                                                 # original scheme is codomain of map 1
  # boolean flags
  resolves_sing::Bool                            # domain not smooth yet?
  is_trivial::Bool                               # codomain already smooth?

  # fields for caching, to be filled during desingularization
  # always carried along to domain(maps[end])) using strict_transform
  ex_div::Vector{AbsIdealSheaf}      # list of exc. divisors arising from individual steps

  # keep track of the normalization steps
  normalization_steps::Vector{Int}

  # fields for caching to be filled a posteriori (on demand, only if partial_res==false)
  underlying_morphism::CompositeCoveredSchemeMorphism{DomainType, CodomainType}
  exceptional_divisor::AbsWeilDivisor
  exceptional_locus::AbsAlgebraicCycle

  function MixedBlowUpSequence(maps::Vector{<:AbsCoveredSchemeMorphism})
    n = length(maps)
    for i in 1:n
      @assert all(x->((x isa BlowupMorphism) || (x isa NormalizationMorphism)), maps) "only blow-ups and normalizations allowed"
    end
    for i in 1:n-1
      @assert domain(maps[i]) === codomain(maps[i+1]) "not a sequence of morphisms"
    end
    resi = new{typeof(domain(maps[end])),typeof(codomain(first(maps)))}(maps)
    resi.normalization_steps = [i for i in 1:n if maps[i] isa NormalizationMorphism]
    return resi
  end

end


##################################################################################################
# getters
##################################################################################################
underlying_morphism(phi::NormalizationMorphism) = phi.underlying_morphism
inclusion_morphisms(phi::NormalizationMorphism) = phi.inclusions
normalization_steps(phi::NormalizationMorphism) = phi.normalization_steps
morphisms(phi::AbsDesingMor) = copy(phi.maps)
morphism(phi::AbsDesingMor,i::Int) = phi.maps[i]
last_map(phi::AbsDesingMor) = phi.maps[end]

exceptional_divisor_list(phi::BlowUpSequence) = phi.ex_div

## entries of ex_div corresponding to normalization steps are only exceptional divisors at the very end
## so only return them at the very end or on specific demand 
function exceptional_divisor_list(phi::MixedBlowUpSequence, seq_unclean::Bool=false)
  phi.resolves_sing || seq_unclean || error("excpetional divisor list not available for intermediate steps.")
  return phi.ex_div
end

normalization_steps(phi::MixedBlowUpSequence) = phi.normalization_steps
 
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

# Example
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

# Example
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

Return a `CartierDivisor` `C` on the `domain` of the emmbedded desingularization morphism `f` 
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
    dim(ex_div_list[i])== -1 && continue            # kick out empty ones
    C = C + weil_divisor(ex_div_list[i])
  end

  f.exceptional_divisor = C
  return C
end

function _exceptional_divisor_non_embedded(f::BlowUpSequence)
  !isdefined(f,:exceptional_divisor_on_X) || return f.exceptional_divisor_on_X

  ex_div_list = exceptional_divisor_list(f)
  C = CartierDivisor(scheme(ex_div_list[1]),ZZ)
  for i in 1:length(ex_div_list)
# do we want to introduce is_empty for divisors?
    dim(ideal_sheaf(ex_div_list[i]))== -1 && continue          # kick out empty ones
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
    dim(E) != -1 || continue                        # kick out empty ones
    C = C + algebraic_cycle(E)
  end

  return C
end

##################################################################################################
# setting values in DesingMors -- Watch out: only place with direct access to fields!!!
##################################################################################################
function add_map!(f::BlowUpSequence, phi::BlowupMorphism)
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
  push!(f.normalization_steps,length(f.maps))
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

function initialize_blow_up_sequence(phi::BlowupMorphism)
  f = BlowUpSequence([phi])
  f.ex_div = [exceptional_divisor(phi)]
  f.is_trivial = is_one(center(phi))
  f.resolves_sing = false                                # we have no information, wether we are done
                                                         # without further computation
  f.is_embedded = false
  return f
end

function initialize_mixed_blow_up_sequence(phi::NormalizationMorphism, I::AbsIdealSheaf)
  f = MixedBlowUpSequence([phi])
  f.ex_div = [pullback(phi,I)]
  f.is_trivial = is_one(I)
  f.resolves_sing = false                                # we have no information, wether we are done
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
    f.resolves_sing = false                              # we have no information, whether we are done
                                                         # without further computation
  else
    f.is_trivial = true
    f.embeddings = [inc, inc]
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
    f.resolves_sing = false                              # we have no information, whether we are done
                                                         # without further computation
  else
    f.is_trivial = true
    f.controlled_transform = I
    f.transform_type = :weak
    f.ex_mult = [0]
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
function embedded_desingularization(f::CoveredClosedEmbedding; algorithm::Symbol=:BEV)
  I_sl = ideal_sheaf_of_singular_locus(domain(f))

  ## trivial case: domain(f) was already smooth
  if is_one(I_sl)
    id_W = identity_blow_up(codomain(f))
    phi = initialize_embedded_blowup_sequence(id_W,f)
    phi.resolves_sing = true
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
  elseif ((dimX == 2) && (algorithm==:Lipman))
    return_value = _desing_lipman(Xnorm, I_sl, f)
  elseif ((dimX == 2) && (algorithm==:Jung))
    error("not implemented yet")
#    second_seq = _desing_jung(Xnorm,f)
#    return_value = extend(phi, second_seq)
  else
    error("not implemented yet")
#    second_seq = forget_embedding(_desing_BEV(Xnorm))
#    return_value = extend(phi, second_seq)
  end

  return return_value
end

function desingularization_only_blowups(X::AbsCoveredScheme; algorithm::Symbol=:Lipman)
  I_sl = ideal_sheaf_of_singular_locus(X)

  ## trivial case: X is already smooth
  if is_one(I_sl)
    id_X = identity_blow_up(X)
    maps = [id_X]
    return_value = BlowUpSequence(maps)
    return_value.resolves_sing = true
    return_value.is_trivial = true
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

function _desing_curve(X::AbsCoveredScheme, I_sl::AbsIdealSheaf)
  ## note: I_sl not unit_ideal_sheaf, because this has been caught before in desingularization(X) 
  decomp = maximal_associated_points(I_sl)
  I = small_generating_set(pop!(decomp))
  current_blow_up = blow_up(I)
  phi = initialize_blow_up_sequence(current_blow_up)::BlowUpSequence
  decomp = [strict_transform(current_blow_up,J) for J in decomp]
  
  I_sl_temp = I_sl
  while !is_one(I_sl_temp)
    while length(decomp) > 0
      I = small_generating_set(pop!(decomp))
      phi = _do_blow_up!(phi,I)
      if length(decomp)>0 
        decomp = [strict_transform(last_map(phi),J) for J in decomp]
      end
    end
    I_sl_temp = ideal_sheaf_of_singular_locus(domain(last_map(phi)))
    decomp = maximal_associated_points(I_sl_temp)
  end

  phi.resolves_sing = true
  return phi
end

function _desing_lipman(X::AbsCoveredScheme, I_sl::AbsIdealSheaf, f::MixedBlowUpSequence)
  dim(X) == 2 || error("Lipman's algorithm is not applicable")

  if dim(I_sl) == 1                         # called for a non-normal X
    Xnorm, phi = normalization(X)
    incs = phi.inclusions
    f = initialize_mixed_blow_up_sequence(phi,I_sl)
    I_sl_temp = ideal_sheaf_of_singular_locus(Xnorm)
  else
    I_sl_temp = I_sl
  end    
  
  decomp = maximal_associated_points(I_sl_temp)

  while !is_one(I_sl_temp)
    while length(decomp) > 0
      I = small_generating_set(pop!(decomp))
      f = _do_blow_up!(f,I)
      if length(decomp)>0 
        decomp = [strict_transform(last_map(f),J) for J in decomp]
      end
    end
    I_sl_temp = ideal_sheaf_of_singular_locus(domain(last_map(f)))
    if dim(I_sl_temp) == 1
      f = _do_normalization!(f)
      I_sl_temp = ideal_sheaf_of_singular_locus(domain(last_map(f)))
    end
    decomp = maximal_associated_points(I_sl_temp)
  end

  f.resolves_sing = true
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
    while length(decomp) > 0
      I = small_generating_set(pop!(decomp))
      phi = _do_blow_up_embedded!(phi,I)
      if length(decomp)>0
        decomp = [strict_transform(last_map(phi),J) for J in decomp]
      end
    end
    last_emb = embeddings(phi)[end]
    I_sl_temp = pushforward(last_emb, ideal_sheaf_of_singular_locus(domain(last_emb)))
    decomp = maximal_associated_points(I_sl_temp)
  end

  phi = _ensure_ncr!(phi)
  phi.resolves_sing = true
  return phi
end

function _ensure_ncr!(f::AbsDesingMor)
  current_divs = exceptional_divisor_list(f)

# this first step can only come into play, if the centers were not determined algorithmically
# it is ensured by all standard desingularization algorithms
  I_bad = non_snc_locus(current_divs)
  while !is_one(I_bad)
    decomp = maximal_associated_points(I_bad)
    while length(decomp)>0
      I = small_generating_set(pop!(decomp))
      f =_do_blow_up_embedded!(f,I)
      if length(decomp)>0
        decomp = [strict_transform(last_map(f),J) for J in decomp]
      end
    end
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
    decomp = maximal_associated_points(pushforward(inc_temp, next_locus ))
    while !is_empty(decomp)
      I = small_generating_set(pop!(decomp))
      f =_do_blow_up_embedded!(f,I)
      I_X = image_ideal(f.embeddings[end])
      current_divs = [strict_transform(last_map(f),J) for J in current_divs]
      push!(current_divs, exceptional_divisor_list(f)[end])
    end
  end

# finally make sure not too many exceptional divisors meet the strict transform in the same point
  n_max = dim(I_X)
  current_divs = copy(exceptional_divisor_list(f))
  _,inter_div = divisor_intersections_with_X(current_divs,I_X)
  while !is_empty(inter_div)
    cent = small_generating_set(pop!(inter_div))
    f =_do_blow_up_embedded!(f,cent)
    if length(inter_div)>0
      inter_div = [strict_transform(last_map(f),J) for J in inter_div]
    end
  end
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
  return _delta_list(I)[end]
end

function locus_of_order_geq_b(I::AbsIdealSheaf, b::Int)
  return _delta_list(I,b)[end]
end

function _delta_ideal_for_order(inc::CoveredClosedEmbedding)

  W = codomain(inc)
  I = image_ideal(inc)

  Delta_dict = IdDict{AbsAffineScheme,Ideal}()

  for U in affine_charts(W)
    XU, inc_U = sub(U,I(U))
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
                                 <: Tuple{Vector{Int64},Vector{Int64},<:RingElem}};
       check::Bool=true)

  W = codomain(inc)                                
  @check is_smooth(W) "codomain of embedding needs to be smooth"
#  @check is_equidimensional(W) "codomain of embedding needs to be equidimensional"
  I_X = small_generating_set(image_ideal(inc))         # ideal sheaf describing X on W

  Delta_dict = IdDict{AbsAffineScheme,Ideal}()
  for U in Cov
    I = I_X(U)
    if is_one(I)
      Delta_dict[U] = I
      continue
    end

    amb_row,amb_col,h = ambient_param_data[U]
    mod_gens = lifted_numerator.(gens(modulus(OO(U))))
    R = ambient_coordinate_ring(U)
    JM = jacobian_matrix(R, mod_gens)
    if length(amb_col) < length(mod_gens)
      JM_essential = JM[:, amb_col]
    else
      JM_essential = JM
    end
    submat_for_minor = JM[amb_row, amb_col]
    Ainv, h2 = pseudo_inv(submat_for_minor)
    h == h2 || error("inconsistent input data")
    JM_essential = JM_essential * Ainv
    I_gens = lifted_numerator.(gens(I))
    JI = jacobian_matrix(I_gens)
    result_mat = h*JI
    for i in 1:length(amb_col)
      for j in 1:ncols(result_mat)
        result_mat[:,j] = result_mat[:,j] - [JM_essential[k,i] * JI[amb_row[i],j] for k in 1:nrows(JM_essential)]
      end
    end
    Ivec = copy(gens(I))
    append!(Ivec,[a for a in OO(U).(collect(result_mat))])
    Delta_dict[U] = ideal(OO(U),Ivec)
  end

  return small_generating_set(IdealSheaf(W,Delta_dict))
end
 
########################################################################
# test for snc                                                         #
########################################################################

function divisor_intersections_with_X(current_div, I_X)
  scheme(I_X) == scheme(current_div[1]) || error("underlying schemes do not match")
  n_max = dim(I_X)

  inter_div_dict = Dict{Vector{Int},Tuple{AbsIdealSheaf,Int}}()
  old_keys = Vector{Int}[]
  empty_keys = Vector{Int}[]
  essential_inter = AbsIdealSheaf[]

# initialization: each divisor + I_X
  for k in 1:length(current_div)
    Idiv = ideal_sheaf(current_div[k]) + I_X
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
  X = scheme(first(divs))
  @assert all(d->scheme(d) === X, divs)
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

function strict_transform(bl::AbsSimpleBlowdownMorphism, inc::ClosedEmbedding)
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
  for psi in morphisms(phi)
    a = strict_transform(psi, a)
  end
  return a
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
