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
  resolves_sing::Bool                            # domain(maps[end]) smooth?
  is_trivial::Bool                               # codomain already smooth?

  # fields for caching, may be filled during computation
  ex_div::Vector{<:EffectiveCartierDivisor}               # list of exc. divisors arising from individual steps
                                                 # lives in domain(maps[end])

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
#  
#  embeddings::Vector{<:AbsCoveredSchemeMorphism} # assert codomain(maps[i])===codomain(embeddings[i]) 
#
#  # boolean flags
#  resolves_sing::Bool                            # domain not smooth yet?
#  is_trivial::Bool                               # codomain already smooth?
#
#  # fields for caching, may be filled during computation
#  ex_div::Vector{<:EffectiveCartierDivisor}               # list of exc. divisors arising from individual steps
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
  f.ex_div = push!(ex_div,Oscar.exceptional_divisor(phi))
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
  f.resolves_sing = false
  return f
end

##################################################################################################
# desingularization workers
##################################################################################################
function embedded_desingularization(f::Oscar.CoveredClosedEmbedding; algorithm::Symbol=:BEV)
  I_sl = Oscar.ideal_sheaf_of_singular_locus(domain(f))

  ## trivial case: domain(f) was already smooth
  if is_one(I_sl)
    id_X = identity_blow_up(codomain(f))
    maps = [id_X] 
    return_value = BlowUpSequence(maps)
    return_value.embeddings = [f,f]
    return_value.resolves_sing = true
    return_value.is_trivial = true
    return return_value
  end

  ## I_sl non-empty, we need to do something
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
  if dim(X) == 1
@show "overriding specified method for curves: use naive method"
    return_value = _desing_curve(X, I_sl)
    return return_value
  end
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

