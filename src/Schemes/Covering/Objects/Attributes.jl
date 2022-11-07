export base_ring_type
export patches, basic_patches, npatches, glueings, glueing_graph

########################################################################
# Type getters                                                         #
########################################################################
base_ring_type(::Type{Covering{T}}) where {T} = T
base_ring_type(C::Covering) = base_ring_type(typeof(C))

### type constructors
#covering_type(::Type{T}) where {T<:Spec} = Covering{T, glueing_type(T)}
#covering_type(X::Spec) = covering_type(typeof(X))

########################################################################
# Basic getters                                                        #
########################################################################
patches(C::Covering) = C.patches
basic_patches(C::Covering) = C.patches
npatches(C::Covering) = length(C.patches)
glueings(C::Covering) = C.glueings
getindex(C::Covering, i::Int) = C.patches[i]
getindex(C::Covering, i::Int, j::Int) = glueings(C)[(patches(C)[i], patches(C)[j])]
getindex(C::Covering, X::AbsSpec, Y::AbsSpec) = glueings(C)[(X, Y)]
#edge_dict(C::Covering) = C.edge_dict

function glueing_graph(C::Covering) 
  if !isdefined(C, :glueing_graph)
    update_glueing_graph(C)
  end
  return C.glueing_graph
end

