
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

function glueing_graph(C::Covering; all_dense::Bool=false)
  if !isdefined(C, :glueing_graph)
    update_glueing_graph(C, all_dense=all_dense)
  end
  return C.glueing_graph
end

function glueing_tree(C::Covering; all_dense::Bool=false)
  if !isdefined(C, :glueing_tree)
    U = copy(patches(C))
    n = npatches(C)
    gg = Graph{Undirected}(n)
    heap = AbsSpec[pop!(U)]
    while !isempty(U)
      W = last(heap) 
      next = [B for (A, B) in keys(glueings(C)) if A === W && any(x->x===B, U)]
      if isempty(next)
        pop!(W)
        isempty(W) && error("scheme is not connected")
        continue
      end
      for B in next
        if all_dense
          add_edge!(gg, C[W], C[B])
          add_edge!(gg, C[B], C[W])
        else
          WB, BW = glueing_domains(C[W, B])
          is_dense(WB) && add_edge!(gg, C[W], C[B])
          is_dense(BW) && add_edge!(gg, C[B], C[W])
        end
      end
      U = [V for V in U if !any(x->x===V, next)]
      heap = vcat(heap, next)
    end
    C.glueing_tree = gg
  end
  return C.glueing_tree
end

@doc raw"""
    decomposition_info(C::Covering)

Return an `IdDict` `D` with the `patches` ``Uᵢ`` of `C` as keys and values a list 
of elements ``fᵢ₁,…,fᵢᵣ ∈ 𝒪(Uᵢ)``. These elements are chosen so that for every 
affine patch `Uᵢ` of ``X`` in the covering `C` the closed subvarieties ``Zᵢ ⊂ Uᵢ`` 
defined by the ``fᵢⱼ`` give rise to a decomposition of ``X`` as a **disjoint** union 
``X = \bigcup_{i} Z_i`` of locally closed subvarieties. 

This information can be used for local computations in any chart ``Uᵢ`` of ``X`` 
as above to focus on phenomena occuring exclusively along ``Zᵢ`` and assuming 
that other cases have been handled by computations in other charts. A key 
application is counting points of zero-dimensional subschemes: To avoid overcounting, 
we need to only consider points in ``Uᵢ`` which are located within ``Zᵢ`` and 
then sum these up to all points in ``X``.

!!! note This attribute might not be defined! Use `has_decomposition_info(C)` to check whether this information is available for the given covering.
"""
function decomposition_info(C::Covering)
  return C.decomp_info
end

function set_decomposition_info!(C::Covering, U::AbsSpec, f::Vector{<:RingElem})
  if !isdefined(C, :decomp_info)
    C.decomp_info = IdDict{AbsSpec, Vector{RingElem}}()
  end
  all(x->parent(x) === OO(U), f) || error("elements do not belong to the correct ring")
  decomposition_info(C)[U] = f
end

function set_decomposition_info!(C::Covering, D::IdDict{<:AbsSpec, <:Vector{<:RingElem}})
  C.decomp_info = D
end

function has_decomposition_info(C::Covering) 
  isdefined(C, :decomp_info) || return false
  all(x->haskey(C.decomp_info, x), patches(C)) || return false
  return true
end
