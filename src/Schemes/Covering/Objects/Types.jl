export Covering

########################################################################
# Coverings for covered schemes                                        #
########################################################################
@Markdown.doc """
    Covering

A covering of a scheme ``X`` by affine charts ``Uᵢ`` which are glued
along isomorphisms ``gᵢⱼ : Uᵢ⊃ Vᵢⱼ →  Vⱼᵢ ⊂ Uⱼ``.

**Note:** The distinction between the different affine charts of the scheme
is made from their hashes. Thus, an affine scheme must not appear more than once
in any covering!
"""
mutable struct Covering{BaseRingType}
  patches::Vector{<:AbsSpec} # the basic affine patches of X
  glueings::IdDict{Tuple{<:AbsSpec, <:AbsSpec}, <:AbsGlueing} # the glueings of the basic affine patches
  affine_refinements::IdDict{<:AbsSpec, <:Vector{<:Tuple{<:SpecOpen, Vector{<:RingElem}}}} # optional lists of refinements
      # of the basic affine patches.
      # These are stored as pairs (U, a) where U is a 'trivial' SpecOpen,
      # meaning that its list of hypersurface equation (f₁,…,fᵣ) has empty
      # intersection in the basic affine patch X and hence satisfies
      # some equality 1 ≡ a₁⋅f₁ + a₂⋅f₂ + … + aᵣ⋅fᵣ on X.
      # Since the coefficients aᵢ of this equality are crucial for computations,
      # we store them in an extra tuple.

  # fields for caching
  glueing_graph::Graph{Undirected}
  transition_graph::Graph{Undirected}
  edge_dict::Dict{Tuple{Int, Int}, Int}

  function Covering(
      patches::Vector{<:AbsSpec},
      glueings::IdDict{Tuple{<:AbsSpec, <:AbsSpec}, <:AbsGlueing};
      check::Bool=true,
      affine_refinements::IdDict{
          <:AbsSpec,
          <:Vector{<:Tuple{<:SpecOpen, <:Vector{<:RingElem}}}
         }=IdDict{AbsSpec, Vector{Tuple{SpecOpen, Vector{RingElem}}}}()
    )
    n = length(patches)
    n > 0 || error("can not glue the empty scheme")
    kk = coefficient_ring(ambient_coordinate_ring(patches[1]))
    for i in 2:n
      kk == coefficient_ring(base_ring(OO(patches[i]))) || error("schemes are not defined over the same base ring")
    end
    # Check that no patch appears twice
    for i in 1:n-1
      for j in i+1:n
        patches[i] === patches[j] && error("affine schemes must not appear twice among the patches")
      end
    end
    for (X, Y) in keys(glueings)
      X in patches || error("glueings are not compatible with the patches")
      Y in patches || error("glueings are not compatible with the patches")
      if haskey(glueings, (Y, X))
        if check
          inverse(glueings[(X, Y)]) == glueings[(Y, X)] || error("glueings are not inverse of each other")
        end
      else
        glueings[(Y, X)] = inverse(glueings[(X, Y)])
      end
    end

    # check the affine refinements
    for U in keys(affine_refinements)
      for (V, a) in affine_refinements[U]
        ambient(V) == U && error("the ambient scheme of the refinement of X must be X")
        U in patches && error("the ambient scheme of the refinement can not be found in the affine patches")
        if check
          isone(OO(U)(sum([c*g for (c, g) in zip(a, gens(U))]))) || error("the patch $V does not cover $U")
        end
      end
    end
    return new{base_ring_type(patches[1])}(patches, glueings, affine_refinements)
  end

  ### the empty covering
  function Covering(kk::Ring)
    return new{typeof(kk)}(Vector{AbsSpec}(), IdDict{Tuple{AbsSpec, AbsSpec}, AbsGlueing}(),
                           IdDict{AbsSpec, Vector{Tuple{SpecOpen, Vector{RingElem}}}}())
  end
end

