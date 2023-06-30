
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
@doc raw"""
    affine_charts(C::Covering)

Return the list of affine charts that make up the `Covering` `C`.
"""
affine_charts(C::Covering) = C.patches
npatches(C::Covering) = length(C.patches)
@doc raw"""
    glueings(C::Covering)

Return a dictionary of glueings of the `affine_chart`s of `C`.

The keys are pairs `(U, V)` of `affine_chart`s. 
One can also use `C[U, V]` to obtain the respective glueing.

**Note:** Glueings are lazy in the sense that they are in general 
only computed when asked for. This method only returns the internal 
cache, but does not try to compute new glueings.
"""
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

@doc raw"""
    decomposition_info(C::Covering)

Return an `IdDict` `D` with the `patches` of `C` as keys and values `D[U]` a list 
of elements ``f₁,…,fᵣ ∈ 𝒪(U)``. These elements are chosen so that for 
``Xᵤ = V(f₁,…,fᵣ) ⊂ U`` the scheme ``X`` covered by ``C`` decomposes as a 
disjoint union ``X = ∪ Xᵤ``.

!!! note This attribute might not be defined!
"""
function decomposition_info(C::Covering)
  return C.decomp_info
end

function set_decomposition_info!(C::Covering, U::AbsSpec, f::Vector{<:RingElem})
  if !isdefined(C, :decomp_info)
    C.decomp_info = IdDict{AbsSpec, Vector{RingElem}}()
  end
  all(x->parent(x) === OO(U), f) || error("elements do not belong to the correct ring")
  decomp_info(C)[U] = f
end

function set_decomposition_info!(C::Covering, D::IdDict{<:AbsSpec, <:Vector{<:RingElem}})
  C.decomp_info = D
end

has_decomposition_info(C::Covering) = isdefined(C, :decomp_info)
