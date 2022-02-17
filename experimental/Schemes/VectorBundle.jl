export VectorBundle
export scheme, covering, transitions, set_name!, has_name, name_of, vector_bundle_type

@attributes mutable struct VectorBundle{
    CoveredSchemeType<:CoveredScheme,
    CoveringType<:Covering,
    SpecType<:Spec,
    MatrixType
  }
  X::CoveredSchemeType
  C::CoveringType
  rank::Int
  transitions::Dict{Tuple{SpecType, SpecType}, MatrixType}

  function VectorBundle(
      X::CoveredScheme, 
      C::Covering, 
      r::Int,
      g::Dict{Tuple{SpecType, SpecType}, MatrixType}; 
      check::Bool=true
    ) where {SpecType<:Spec, MatrixType}
    C in coverings(X) || error("covering is not listed")
    for (U, V) in keys(g)
      (U in patches(C) && V in patches(C)) || error("affine patches not listed")
      A = g[(U, V)]
      ncols(A) == nrows(A) == r || error("matrix does not have the correct size")
      if check 
        G = C[U, V]
        W1, W2 = glueing_domains(G)
        h1, h2 = glueing_morphisms(G)
        ### the following method is not yet implemented
        #B = map_entries(tmp -> pullback(h1, tmp), g[(V, U)])
        #inv(A) == B || error("transition matrices are not inverse to each other")
      end
    end
    if check
      #TODO: Do some more sophisticated checks
    end
    return new{typeof(X), typeof(C), affine_patch_type(X), MatrixType}(X, C, r, g)
  end
end

vector_bundle_type(::Type{T}) where {T<:CoveredScheme} = VectorBundle{
    T, 
    covering_type(T), 
    affine_patch_type(T), 
    AbstractAlgebra.Generic.MatSpaceElem{
        elem_type(SpecOpenRing{affine_patch_type(T), open_subset_type(affine_patch_type(T))})}
   }
vector_bundle_type(C::CoveredScheme) = vector_bundle_type(typeof(C))
vector_bundle_type(::Type{T}) where {T<:ProjectiveScheme} = vector_bundle_type(covered_scheme_type(T))
vector_bundle_type(P::ProjectiveScheme) = vector_bundle_type(typeof(P))

scheme(E::VectorBundle) = E.X
covering(E::VectorBundle) = E.C
transitions(E::VectorBundle) = E.transitions
getindex(E::VectorBundle, U::SpecType, V::SpecType) where {SpecType<:Spec} = E.transitions[(U, V)]

set_name!(E::VectorBundle, name::String) = set_attribute!(E, :name, name)
has_name(E::VectorBundle) = has_attribute(E, :name)
name_of(E::VectorBundle) = get_attribute!(E, :name)::String

function OO(P::ProjectiveScheme, d::Int)
  PC = as_covered_scheme(P)
  if !has_attribute(PC, :tautological_bundles)
    bundle_dict = Dict{Int, vector_bundle_type(PC)}()
    set_attribute!(PC, :tautological_bundles, bundle_dict)
  end
  bundle_dict = get_attribute(PC, :tautological_bundles)::Dict{Int, vector_bundle_type(PC)}
  if !haskey(bundle_dict, d) 
    C = coverings(PC)[1]
    SpecType = affine_patch_type(PC)
    MatrixType = AbstractAlgebra.Generic.MatSpaceElem{elem_type(SpecOpenRing{SpecType, open_subset_type(SpecType)})}
    trans = Dict{Tuple{affine_patch_type(PC), affine_patch_type(PC)}, MatrixType}()
    for i in 1:npatches(C)-1
      for j in i+1:npatches(C)
        X = C[i]
        Y = C[j]
        G = C[i, j]
        U, V = glueing_domains(G)
        MU = MatrixSpace(OO(U), 1, 1)
        MV = MatrixSpace(OO(V), 1, 1)
        d >= 0 ? trans[(X, Y)] = MU(OO(U)((gens(OO(X))[j-1])^d)) : trans[(X, Y)] = MU(inv(OO(U)(gens(OO(X))[j-1]))^(-d))
        d >= 0 ? trans[(Y, X)] = MV(OO(V)((gens(OO(Y))[i])^d)) : trans[(Y, X)] = MV(inv(OO(V)(gens(OO(Y))[i]))^(-d))
      end
    end
    bundle_dict[d] = VectorBundle(PC, C, 1, trans)
  end
  return bundle_dict[d]
end
      
