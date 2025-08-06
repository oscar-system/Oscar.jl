
########################################################################
# Type getters                                                         #
########################################################################
base_ring_type(::Type{Covering{T}}) where {T} = T

### type constructors
#covering_type(::Type{T}) where {T<:AffineScheme} = Covering{T, gluing_type(T)}
#covering_type(X::AffineScheme) = covering_type(typeof(X))

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
number_of_patches(C::Covering) = length(C.patches)
@doc raw"""
    gluings(C::Covering)

Return a dictionary of gluings of the `affine_chart`s of `C`.

The keys are pairs `(U, V)` of `affine_chart`s.
One can also use `C[U, V]` to obtain the respective gluing.

**Note:** Gluings are lazy in the sense that they are in general
only computed when asked for. This method only returns the internal
cache, but does not try to compute new gluings.
"""
gluings(C::Covering) = C.gluings
getindex(C::Covering, i::Int) = C.patches[i]
getindex(C::Covering, i::Int, j::Int) = gluings(C)[(patches(C)[i], patches(C)[j])]
function getindex(C::Covering, X::AbsAffineScheme, Y::AbsAffineScheme) 
  @vprintln :Gluing 1 "Requesting gluing ($(indexin(X,C)[1]),$(indexin(Y,C)[1]))"
  return gluings(C)[(X, Y)]
end
#edge_dict(C::Covering) = C.edge_dict

function gluing_graph(C::Covering; all_dense::Bool=false)
  if !isdefined(C, :gluing_graph)
    update_gluing_graph(C, all_dense=all_dense)
  end
  return C.gluing_graph
end

@doc raw"""
    decomposition_info(C::Covering)

Return an `IdDict` `D` with the `patches` ``Uᵢ`` of `C` as keys and values a list
of elements ``fᵢ₁,…,fᵢᵣ ∈ 𝒪(Uᵢ)``. These elements are chosen so that for every
affine patch `Uᵢ` of ``X`` in the covering `C` the closed subvarieties ``Zᵢ ⊂ Uᵢ``
defined by the ``fᵢⱼ`` give rise to a decomposition of ``X`` as a **disjoint** union
``X = \bigcup_{i} Z_i`` of locally closed subvarieties.

This information can be used for local computations in any chart ``Uᵢ`` of ``X``
as above to focus on phenomena occurring exclusively along ``Zᵢ`` and assuming
that other cases have been handled by computations in other charts. A key
application is counting points of zero-dimensional subschemes: To avoid overcounting,
we need to only consider points in ``Uᵢ`` which are located within ``Zᵢ`` and
then sum these up to all points in ``X``.

!!! note This attribute might not be defined! Use `has_decomposition_info(C)` to check whether this information is available for the given covering.
"""
function decomposition_info(C::Covering)
  return C.decomp_info
end

function set_decomposition_info!(C::Covering, U::AbsAffineScheme, f::Vector{<:RingElem})
  if !isdefined(C, :decomp_info)
    C.decomp_info = IdDict{AbsAffineScheme, Vector{RingElem}}()
  end
  all(x->parent(x) === OO(U), f) || error("elements do not belong to the correct ring")
  decomposition_info(C)[U] = f
end

function set_decomposition_info!(C::Covering, D::IdDict{<:AbsAffineScheme, <:Vector{<:RingElem}})
  C.decomp_info = D
end

function has_decomposition_info(C::Covering)
  isdefined(C, :decomp_info) || return false
  all(x->haskey(C.decomp_info, x), patches(C)) || return false
  return true
end

function inherit_decomposition_info!(
    X::AbsCoveredScheme, ref_cov::Covering;
    orig_cov::Covering=default_covering(X)
  )
  !has_decomposition_info(orig_cov) && return ref_cov
  OX = OO(X)

  # For every chart U of `ref_cov`, we find its patch `V` in `orig_cov`, 
  # together with a list of equations h₁,…,hₙ ∈ 𝒪 (V) such that 
  # D(h₁⋅…⋅hₙ) ≅ U.
  decomp_dict = IdDict{AbsAffineScheme, Tuple{AbsAffineScheme, Vector{RingElem}}}()
  for U in patches(ref_cov)
    inc_U, d_U = _find_chart(U, orig_cov)
    decomp_dict[U] = (codomain(inc_U), d_U)
  end

  for V in patches(orig_cov)
    # Collect the patches in `ref_cov` refining V
    V_ref = [(U, h) for (U, (UV, h)) in decomp_dict if UV === V]
    _compl(a) = sum(length(lifted_numerator(b)) + length(lifted_denominator(b)) for b in a[2]; init=0)
    V_ref = sort!(V_ref, by=_compl)
    # Start out from the original decomposition info
    dec_inf = copy(decomposition_info(orig_cov)[V])
    for (i, (U, h)) in enumerate(V_ref)
      # Cast the already made decomposition info down to U
      tmp = elem_type(OO(U))[OX(V, U)(i) for i in dec_inf] # help the compiler
      # Append the equations for the previously covered patches
      for j in 1:i-1
        _, hh = V_ref[j]
        if is_empty(hh) # The patch for hh already covered everything; this one can hence be discarded.
          push!(tmp, one(OO(U)))
          break
        end
        surplus = OX(V, U).(hh)
        push!(tmp, prod(surplus; init=one(OO(U))))
      end
      set_decomposition_info!(ref_cov, U, tmp)
    end
  end
  return ref_cov
end


