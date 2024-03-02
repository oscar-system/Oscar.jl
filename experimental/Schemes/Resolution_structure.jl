## Warnung: show auf desingMor geht noch nicht!!!
export _desing_curve

#####################################################################################################
# Desingularization morphism: birational map between covered schemes with smooth domain
#####################################################################################################

@doc raw"""
    BlowUpSequence{
    DomainType<:AbsCoveredScheme,
    CodomainType<:AbsCoveredScheme
   } <: AbsDesingMor{
                                 DomainType,
                                 CodomainType,
                                }


"""
@attributes mutable struct BlowUpSequence{
                                          DomainType<:AbsCoveredScheme,
                                          CodomainType<:AbsCoveredScheme
                                         }<:AbsBlowdownMorphism{
                                                                DomainType, CodomainType, 
                                                                BlowUpSequence{DomainType, CodomainType}
                                                               }
  maps::Vector{<:BlowupMorphism}                 # count right to left:
                                                 # original scheme is codomain of map 1
  
  embeddings::Vector{<:AbsCoveredSchemeMorphism} # if set,
                                                 # assert codomain(maps[i])===codomain(embeddings[i]) 
  # boolean flags
  is_embedded::Bool                              # do not set embeddings, ex_mult, controlled_transform etc
                                                 #     if is_embedded == false
  resolves_sing::Bool                            # domain(maps[end]) smooth?
  is_trivial::Bool                               # codomain already smooth?
  transform_type::Symbol                         # can be :strict, :weak or :control
                                                 #     only relevant for is_embedded == true

  # fields for caching, may be filled during computation
  ex_div::Vector{<:EffectiveCartierDivisor}      # list of exc. divisors arising from individual steps
                                                 # lives in domain(maps[end])
  ex_mult::Vector{Int64]                         # multiplicities of exceptional divisors removed from
                                                 # controlled or weak transform, not set for is_embedded == false
                                                 # and transform_type == strict
  controlled_transform::IdealSheaf               # holds weak or controlled transform according to transform_type

  # fields for caching to be filled a posteriori (on demand, only if partial_res==false)
  composed_map::AbsCoveredSchemeMorphism        
  exceptional_divisor::WeilDivisor               # exceptional divisor of composed_map
  exceptional_divisor_on_X::WeilDivisor          # exceptional divisor of composed_map
                                                 # restricted to domain(embeddings[end])

  function BlowUpSequence(maps::Vector{<:BlowupMorphism})
    n = length(maps)
    for i in 1:n-1
      @assert domain(maps[i]) === codomain(maps[i+1]) "not a sequence of morphisms"
    end
    return new{typeof(domain(maps[end])),typeof(codomain(first(maps)))}(maps)
  end
end

# Fehlt: NormalizationMorphism fuer Schemata -- muessten wir haben, sobald wir Lipman machen wollen
#
#@attributes mutable struct LipmanStyleSequence{
#    DomainType<:AbsCoveredScheme,
#    CodomainType<:AbsCoveredScheme
#   } <: AbsDesingMor{
#                                 DomainType,
#                                 CodomainType,
#                    }
#  maps::Vector{:<union{BlowupMorphism,NormalizationMorphism}}        # count right to left:
#                                                 # original scheme is codomain of map 1
#  # boolean flags
#  resolves_sing::Bool                            # domain not smooth yet?
#  is_trivial::Bool                               # codomain already smooth?
#
#  # fields for caching, may be filled during computation
#  ex_div::Vector{<:EffectiveCartierDivisor}      # list of exc. divisors arising from individual steps
                                                  # in domain(maps[end])
#
#  # fields for caching to be filled a posteriori (on demand, only if partial_res==false)
#  composed_map::AbsCoveredSchemeMorphism        
#  exceptional_divisor::WeilDivisor          
#
#  function LipmanStyleSequence(maps::Vector{<:AbsCoveredSchemeMorphism})
#    n = length(maps)
#    for i in 1:n-1
#      @assert domain(maps[i]) === codomain(maps[i+1]) "not a sequence of morphisms"
#    end
#    return new{typeof(domain(maps[end])),typeof(codomain(first(maps)))}(maps)
#  end
#end


##################################################################################################
# getters
##################################################################################################
maps(phi::AbsDesingMor) = phi.maps
last_map(phi::AbsDesingMor) = phi.maps[end]
exceptional_divisor_list(phi::BlowUpSequence) = phi.ex_div  ## derzeit Liste von Eff. Cartier Div.

## do not use!!! (for forwarding and certain emergenies)
function underlying_morphism(phi::AbsDesingMor)
  if !isdefined(phi, :composed_map)
    len=length(maps(phi))
    result=underlying_morphism(maps(phi)[1])
    for i in 2:len
      result = compose(underlying_morphism(maps(phi)[i]), result)
    end
    phi.composed_map = result
  end
  return phi.composed_map
end

##################################################################################################
# setting values in DesingMors -- Watch out: only place with direct access to fields!!!
##################################################################################################
function add_map!(f::AbsDesingMor,phi::BlowupMorphism)
  push!(f.maps, phi)
  ex_div = [strict_transform(phi,E) for E in f.ex_div[1:end]]
  push!(ex_div,Oscar.exceptional_divisor(phi))
  f.ex_div = ex_div
  return f
end

function initialize_blow_up_sequence(phi::BlowupMorphism)
  f = BlowUpSequence([phi])
  f.ex_div = [Oscar.exceptional_divisor(phi)]
  if !is_one(center(phi))
    f.is_trivial = false
  else
    f.is_trivial = true
  end
  f.resolves_sing = false                                # we have no information, wether we are done
                                                         # without further computation
  f.is_embedded = false
  return f
end

function add_map_embedded!(f::AbsDesingMor, phi::BlowupMorphism)
  push!(f.maps, phi)
  ex_div = [strict_transform(phi,E) for E in f.ex_div[1:end]]
  push!(ex_div,Oscar.exceptional_divisor(phi))
  f.ex_div = ex_div
  if f.transform_type == :strict
    X_strict,inc_strict,_ = strict_transform(phi,f.embeddings[end])
    push!(f.embeddings, inc_strict)
  ifelse f.transform_type == :weak
    I_trans,b = weak_transform_with_multiplicity(phi,f.controlled_transform)
    push!(f.ex_mult,b)
    f.controlled_transform = I_trans
  else
    I_trans = controlled_transform(phi, f.controlled_transform, f.ex_mult[end])
    f.controlled_transform = I_trans
    push!(f.ex_mult,f.ex_mult[end])
  end
  return f
end

function initialize_embedded_blowup_sequence(phi::BlowupMorphism, inc::CoveredClosedEmbedding)
  f = BlowUpSequence([phi])
  f.ex_div = [Oscar.exceptional_divisor(phi)]
  f.is_embedded = true
  f.transform_type = :strict
  if !is_one(center(phi))
    f.is_trivial = false
    X_strict,inc_strict,_ = strict_transform(phi,inc)
    f.embeddings = [f,inc_strict]
    f.resolves_sing = false                              # we have no information, whether we are done
                                                         # without further computation
    end
  else
    f.is_trivial = true
    f.embeddings = [inc,inc]
    f.resolves_sing = false
  end
  return f
end

function initialize_embedded_blowup_sequence(phi::BlowupMorphism, I::IdealSheaf, b::Int)
  f = BlowUpSequence([phi])
  f.ex_div = [Oscar.exceptional_divisor(phi)]
  f.is_embedded = true
  if !is_one(center(phi))
    f.is_trivial = false
    if b == 0
      I_trans,b = weak_transform_with_multiplicity(phi,I)
      f.transform_type = :weak
    ifelse b > 0
      I_trans = controlled_transform(phi,I,b)
      f.transform_type = :controlled
    end
    f.controlled_transform = I_trans                     # CAUTION: b is considered set once and for all
    f.ex_mult = [b]
    f.resolves_sing = false                              # we have no information, whether we are done
                                                         # without further computation
    end
  else
    f.is_trivial = true
    f.controlled_transform = I
    f.transform_type = :weak
    f.ex_mult = [0]
    f.resolves_sing = false
  end
  return f
end


##################################################################################################
# desingularization workers
##################################################################################################
function embedded_desingularization(f::Oscar.CoveredClosedEmbedding; algorithm::Symbol=:BEV)
  I_sl = Oscar.ideal_sheaf_of_singular_locus(domain(f))

  ## trivial case: domain(f) was already smooth
  if is_one(I_sl)
    id_W = identity_blow_up(codomain(f))
    phi = initialize_embedded_blowup_sequence(id_W)
    phi.resolves_sing = true
    return phi
  end

  ## I_sl non-empty, we need to do something
  dimX = dim(domain(f))
  if dimX == 1
@show "overriding algorithm for curve case"
    return _desing_emb_curve(f)
#  elseif ((dimX == 2) && (algorithm == :CJS))
#    return _desing_CJS(f)
#  elseif (algorithm == :BEV)
#    return _desing_BEV(f)
  end
# here the keyword algorithm ensures that the desired method is called
  error("not implemented yet")
end

function desingularization(X::AbsCoveredScheme; algorithm::Symbol=:Lipman)
  I_sl = Oscar.ideal_sheaf_of_singular_locus(X)
  
  ## trivial case: X is already smooth
  if is_one(I_sl)
    id_X = identity_blow_up(X)
    maps = [id_X] 
    return_value = BlowUpSequence(maps)
    return_value.resolves_sing = true
    return_value.is_trivial = true
    return return_value
  end

  ## I_sl non-empty, we need to do something 
# here the keyword algorithm ensures that the desired method is called
  dimX = dim(X)
  if dimX == 1
@show "overriding specified method for curves: use naive method"
    return_value = _desing_curve(X, I_sl)
  end
#  if ((dimX == 2) && (algorithm==:Lipman))
#    error("not implemented yet")
#    return_value = _desing_lipman(X, I_sl)
#    return return_value
#  end
#  if ((dimX == 2) && (algorithm==:Jung))
#    error("not implemented yet")
#    return_value = _desing_jung(X)
#   end       
  error("not implemented yet")    
end

function _desing_curve(X::AbsCoveredScheme, I_sl::IdealSheaf)
  ## note: I_sl not unit_ideal_sheaf, because this has been caught before in desingularization(X) 
  decomp = Oscar.maximal_associated_points(I_sl)
  I = pop!(decomp)
  current_blow_up = blow_up(I)
  f = initialize_blow_up_sequence(current_blow_up)
  decomp = [strict_transform(current_blow_up,J) for J in decomp]
  
  while !is_one(I_sl) 
    while length(decomp) > 0
      I = pop!(decomp)
      f = _do_blow_up(f,I)
      if length(decomp)>0 
        decomp = [strict_transform(last_map(f),J) for J in decomp]
      end
    end
    I_sl = Oscar.ideal_sheaf_of_singular_locus(domain(last_map(f)))
    decomp = Oscar.maximal_associated_points(I_sl)
  end

  f.resolves_sing = true
  return(f)
end

function _do_blow_up(f::AbsDesingMor, cent::IdealSheaf)
  old_sequence = maps(f)
  X = domain(old_sequence[end])
  X === scheme(cent) || error("center needs to be defined on same scheme")
  current_blow_up = blow_up(cent,var_name=string("v", length(old_sequence), "_"))
  add_map!(f, current_blow_up)
  return(f)
end


###################################################################################################
# Should go to IdealSheaf.jl, when PR is ready to merge
###################################################################################################

function unit_ideal_sheaf(X::AbsCoveredScheme)
  dd = IdDict{AbsSpec, Ideal}(U=>ideal(OO(U), [one(OO(U))]) for U in affine_patches(X))
  return IdealSheaf(X, dd, check=false)
end

function zero_ideal_sheaf(X::AbsCoveredScheme)
  dd = IdDict{AbsSpec, Ideal}(U=>ideal(OO(U), elem_type(OO(U))[]) for U in affine_patches(X))
  return IdealSheaf(X, dd, check=false)
end

function identity_blow_up(X::AbsCoveredScheme)
  f = BlowupMorphism(X,unit_ideal_sheaf(X))
  return f
end

