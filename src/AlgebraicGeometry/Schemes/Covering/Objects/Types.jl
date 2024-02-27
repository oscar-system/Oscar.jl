
########################################################################
# Coverings for covered schemes                                        #
########################################################################
@doc raw"""
    Covering

A covering of a scheme ``X`` by affine charts ``Uᵢ`` which are glued
along isomorphisms ``gᵢⱼ : Uᵢ⊃ Vᵢⱼ →  Vⱼᵢ ⊂ Uⱼ``.

**Note:** The distinction between the different affine charts of the scheme
is made from their hashes. Thus, an affine scheme must not appear more than once
in any covering!
"""
mutable struct Covering{BaseRingType}
  patches::Vector{<:AbsAffineScheme} # the basic affine patches of X
  gluings::IdDict{Tuple{<:AbsAffineScheme, <:AbsAffineScheme}, <:AbsGluing} # the gluings of the basic affine patches
  affine_refinements::IdDict{<:AbsAffineScheme, <:Vector{<:Tuple{<:AffineSchemeOpenSubscheme, Vector{<:RingElem}}}} # optional lists of refinements
      # of the basic affine patches.
      # These are stored as pairs (U, a) where U is a 'trivial' AffineSchemeOpenSubscheme,
      # meaning that its list of hypersurface equation (f₁,…,fᵣ) has empty
      # intersection in the basic affine patch X and hence satisfies
      # some equality 1 ≡ a₁⋅f₁ + a₂⋅f₂ + … + aᵣ⋅fᵣ on X.
      # Since the coefficients aᵢ of this equality are crucial for computations,
      # we store them in an extra tuple.

  # fields for caching
  gluing_graph::Graph{Undirected}
  transition_graph::Graph{Undirected}
  edge_dict::Dict{Tuple{Int, Int}, Int}
  decomp_info::IdDict{<:AbsAffineScheme, <:Vector{<:RingElem}}

  function Covering(
      patches::Vector{<:AbsAffineScheme},
      gluings::IdDict{<:Tuple{<:AbsAffineScheme, <:AbsAffineScheme}, <:AbsGluing};
      check::Bool=true,
      affine_refinements::IdDict{
          <:AbsAffineScheme,
          <:Vector{<:Tuple{<:AffineSchemeOpenSubscheme, <:Vector{<:RingElem}}}
         }=IdDict{AbsAffineScheme, Vector{Tuple{AffineSchemeOpenSubscheme, Vector{RingElem}}}}()
    )
    n = length(patches)
    n > 0 || error("can not glue the empty scheme")
    kk = coefficient_ring(ambient_coordinate_ring(patches[1]))
    for i in 2:n
      kk == coefficient_ring(ambient_coordinate_ring(patches[i])) || error("schemes are not defined over the same base ring")
    end
    # Check that no patch appears twice
    for i in 1:n-1
      for j in i+1:n
        patches[i] === patches[j] && error("affine schemes must not appear twice among the patches")
      end
    end
    for (X, Y) in keys(gluings)
      any(x->x===X, patches) || error("gluings are not compatible with the patches")
      any(y->y===Y, patches) || error("gluings are not compatible with the patches")
      if haskey(gluings, (Y, X))
        @check inverse(gluings[(X, Y)]) == gluings[(Y, X)] "gluings are not inverse of each other"
      else
        gluings[(Y, X)] = inverse(gluings[(X, Y)])
      end
    end

    # check the affine refinements
    for U in keys(affine_refinements)
      for (V, a) in affine_refinements[U]
        ambient(V) == U && error("the ambient scheme of the refinement of X must be X")
        any(u->u===U, patches) && error("the ambient scheme of the refinement can not be found in the affine patches")
        @check isone(OO(U)(sum([c*g for (c, g) in zip(a, gens(U))]))) "the patch $V does not cover $U"
      end
    end
    return new{base_ring_type(patches[1])}(patches, gluings, affine_refinements)
  end

  ### the empty covering
  function Covering(kk::Ring)
    return new{typeof(kk)}(Vector{AbsAffineScheme}(), IdDict{Tuple{AbsAffineScheme, AbsAffineScheme}, AbsGluing}(),
                           IdDict{AbsAffineScheme, Vector{Tuple{AffineSchemeOpenSubscheme, Vector{RingElem}}}}())
  end
end

