export CoherentSheaf
export free_module

mutable struct CoherentSheaf{SchemeType<:CoveredScheme, 
                             CoveringType<:Covering, 
                             SpecType<:Spec,
                             SubQuoType<:SubQuo
                            }
  X::SchemeType
  C::CoveringType
  M::Dict{SpecType, SubQuoType}
  T::Dict{Tuple{SpecType, SpecType}, Tuple{SubQuoHom, SubQuoHom, SubQuoHom, SubQuoHom}}

  function CoherentSheaf(
      X::CoveredScheme, 
      C::Covering, 
      M::Dict{SpecType, SubQuoType},
      T::Dict{Tuple{SpecType, SpecType}, Tuple{SubQuoHom, SubQuoHom, SubQuoHom, SubQuoHom}};
      check::Bool=true
    ) where {
             SpecType<:Spec,
             SubQuoType<:SubQuo,
            }
    C in coverings(X) || error("covering not found")
    for U in keys(M)
      U in C || error("patch $U not found in covering")
    end
    is_complete_covering(keys(M), C) || error("local module representations do not cover the scheme")
    if check
      # See whether the local compatibility is given on every triple overlap
    end
    for U in keys(M)
      for V in keys(M)
        (U, V) in keys(T) || error("not all transitions are given")
      end
    end
    return new{typeof(X), typeof(C), SpecType, SubQuoType}(X, C, M, T)
  end
end
    
function is_complete_covering(W::Vector{SpecType}, C::Covering) where {SpecType<:Spec}
  for U in basic_patches(C)
    if !(U in W)
      haskey(affine_refinements(C), U) || return false
      Uref = affine_refinements(C)[U]
      found = false
      for (V, a) in Uref
        if all(x->(x in W), affine_patches(V)) 
          found = true
          break
        end
      end
      !found && return false
    end
  end
  return true
end

function free_module(X::CoveredScheme, r::Integer)
  C = default_covering(X)
  SpecType = affine_patch_type(X)
  RingType = ring_type(SpecType)
  SubQuoType = subquo_type(RingType)
  module_dict = Dict{SpecType, SubQuoType}()
  for U in basic_patches(C)
    F = FreeMod(OO(U), r)
    M, _ = sub(F, gens(F))
    module_dict[U] = M
  end
  glueing_dict = Dict{Tuple{SpecType, SpecType}, Tuple{SubQuoHom, SubQuoHom, SubQuoHom, SubQuoHom}}()
  for U in basic_patches(C)
    for V in basic_patches(C)
      i, j, k, l = intersect(U, V, C)
      UnV = domain(k)
      VnU = domain(l)
      MUres, iM = change_base_ring(OO(UnV), module_dict[U])
      MVres, jM = change_base_ring(OO(VnU), module_dict[V])
      kM = hom(MVres, MUres, pullback(k))
      lM = hom(MUres, MVres, pullback(l))
      glueing_dict[(U, V)] = (iM, jM, kM, lM)
    end
  end
  return CoherentSheaf(X, C, module_dict, glueing_dict)
end

