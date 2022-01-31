export IdealSheaf

export scheme, covering, getindex, subscheme

mutable struct IdealSheaf{
    CoveredSchemeType<:CoveredScheme,
    CoveringType<:Covering,
    SpecType<:Spec,
    RingElemType<:MPolyElem
  }
  X::CoveredSchemeType # Der parent
  C::CoveringType
  ideal_gens::Dict{SpecType, Vector{RingElemType}}

  # fields for caching
  cached_ideals::Dict{SpecType, Any}

  function IdealSheaf(
      X::CoveredSchemeType, 
      C::CoveringType, 
      I::Dict{SpecType, Vector{RingElemType}};
      check::Bool=true
    ) where {
      CoveredSchemeType<:CoveredScheme,
      CoveringType<:Covering,
      SpecType<:Spec,
      RingElemType<:MPolyElem
    }
    if check
      C in coverings(X) || error("reference covering not found")
      for X in keys(I)
        X in patches(C) || error("affine patch not found")
        for f in I[X]
          parent(f) == base_ring(OO(X)) || error("element does not belong to the correct ring")
        end
      end
    end
    return new{CoveredSchemeType, CoveringType, SpecType, RingElemType}(X, C, I)
  end
end

scheme(I::IdealSheaf) = I.X
covering(I::IdealSheaf) = I.C
getindex(I::IdealSheaf, U::Spec) = I.ideal_gens[U]

function IdealSheaf(X::ProjectiveScheme, g::Vector{RingElemType}) where {RingElemType<:MPolyElem_dec}
  X_covered = as_covered_scheme(X)
  C = standard_covering(X)
  r = fiber_dimension(X)
  #TODO: Make this independent of the concrete instance of the covering
  I = Dict{affine_patch_type(X), Vector{elem_type(base_ring(OO(C[1])))}}()
  for i in 0:r
    @show dehomogenize(X, g, i)
    @show typeof(dehomogenize(X, g, i))

    I[C[i+1]] = lifted_numerator.(dehomogenize(X, g, i))
  end
  return IdealSheaf(X_covered, C, I, check=false)
end
  

### Given an ideal sheaf I, return the associated 
# subscheme
#
# **Note:** This must be cached!
function subscheme(I::IdealSheaf) 
  X = scheme(I)
  C = covering(I)
  new_patches = [subscheme(C[i], I[C[i]]) for i in 1:npatches(C)]
  new_glueings = Dict{Tuple{affine_patch_type(C), affine_patch_type(C)}, glueing_type(C)}()
  for (U, V) in keys(glueings(C))
    i = C[U]
    j = C[V]
    Unew = new_patches[i]
    Vnew = new_patches[j]
    new_glueings[(Unew, Vnew)] = restriction(C[U, V], Unew, Vnew)
    new_glueings[(Vnew, Unew)] = inverse(new_glueings[(Unew, Vnew)])
  end
  Cnew = Covering(new_patches, new_glueings)
  return CoveredScheme(Cnew)
end
