mutable struct SingularStandardBasis{
    RingElemType<:MPolyElem,
    FreeModType<:FreeMod,
    SubModType<:ModuleFP,
    ModElemType<:FreeModElem, 
    OrderingType<:Oscar.Ordering.ModOrdering
   } <: AbsModuleStdBasis{MPolyElem}
  F::FreeModType
  M::SubModType
  o::OrderingType
  gens::Vector{ModElemType}

  # fields for caching
  elements::Vector{ModElemType}
  quotient_module::ModuleFP
  quotient_projection::Hecke.Map

  g::Oscar.ModuleGens

  function SingularStandardBasis(
      F::FreeMod,
      v::Vector{ModElemType},
      o::OrderingType; 
      generated_submodule::Oscar.SubModuleOfFreeModule=Oscar.SubModuleFreeModule(F, v),
      check::Bool=true
    ) where {
      ModElemType<:FreeModElem, 
      OrderingType<:Oscar.Ordering.ModOrdering
    }
    all(x->(parent(x) == F), v) || error("module elements do not have the correct parent")
    ambient_free_module(generated_submodule) == F || error("submodule does not have the correct ambient module")
    if check
      # Full membership does not make sense here, because its implementation will eventually 
      # be based on this code.
      all(x->(v in gens(generated_submodule)), v) || error("can not check whether the given module elements belong to the given submodule") 
    end
    R = base_ring(F)
    return new{elem_type(R), typeof(F), typeof(generated_submodule), ModElemType, typeof(o)}(F, generated_submodule, o, v)
  end
end

ordering(G::SingularStandardBasis) = G.o
ambient_free_module(G::SingularStandardBasis) = G.F

function sub(G::SingularStandardBasis; cached::Bool=true) 
  cached && return G.M, inclusion_map(M)
  return sub(G.M, G.F)
end
    
function sub(G::SingularStandardBasis; cached::Bool=true) 
  if cached 
    if isdefined(G, :quotient_module) 
      return G.quotient_module, G.quotient_projection
    end
    N, p = quo(ambient_free_module(G), G.M)
    G.quotient_module = N
    G.quotient_projection = p
    return N, p
  end
  return sub(G.M, G.F)
end

function elements(G::SingularStandardBasis)
  if !isdefined(G. :elements) 
    mg = Oscar.ModuleGens(gens(G))
    singular_assure(mg) # How can we pass on the ordering???
    mgstd = Oscar.ModuleGens(ambient_free_module(G), Singular.std(mg.S))
    G.g = mgstd
    G.elements = mgstd.O
  end
  return G.elements
end

