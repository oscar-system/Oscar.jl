@doc raw"""
    MaxContactObject{BaseRingType}

Internal data type of desingularization! Do not use outside!

This object comprises an ambient Covered scheme ``W `` covered by
a covering ``C``, which is refined in such a way that on each chart
a hypersurface of maximal contact may be chosen uniformly, and
internal data on the maximal contact.
"""
@attributes mutable struct MaxContactObject{BaseRingType}
  W_orig::AbsCoveredScheme                  ## the true ambient space
  C::Covering                               ## refined covering of W allowing max contact hypersurfaces
  max_contact_data::IdDict{<:AbsAffineScheme,
                       <:Tuple{Vector{MPolyRingElem},Vector{Int}}}
                                            ## first vector: generators of maximal contact ambient space
                                            ## second vector: order in which maximal contact hypersurface
                                            ##        max_contact_data(U)[1][i] occurred
                                            ## note: original ambient equations are marked by order 0
  ambient_param_data::IdDict{<:AbsAffineScheme,<:Tuple{Vector{Int},Vector{Tuple{RingElem,Int}},
                       Vector{<:MatrixElem}}}
                                            ## vector: indices of rows in minor i.e.dependent variables
                                            ## ring element: data on minor selected by the row indices
                                            ##      (minor of a single max contact step,blocksize of step)
                                            ##      columns of minor follow order of max_contact_data(U)[1]
                                            ## matrix: jacobian matrix * pseudoinverse
                                            ##      for jacobian matrix of max_contact_data(U)[1]
                                            ##      again in blocks corresponding to the steps

  function MaxContactObject(
               W::AbsCoveredScheme,
               Cov::Covering,
               max_contact_data::IdDict{<:AbsAffineScheme,<:Tuple{Vector{MPolyRingElem},Vector{Int}}},
               ambient_param_data::IdDict{<:AbsAffineScheme,
                                          <:Tuple{Vector{Int},Vector{Tuple{RingElem,Int}}, Vector{<:MatrixElem}}})
    return new{base_ring_type(W)}(W,Cov,max_contact_data, ambient_param_data)
  end
end