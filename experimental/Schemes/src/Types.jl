@doc raw"""
    MaxContactChart{BaseRingType}

Internal data type of desingularization! Do not use outside!

This object comprises the local data of a chart, which arose
by hypersurfaces of maximal contact.
"""
@attributes mutable struct MaxContactChart{BaseRingType}
  U::AbsAffineScheme
  ambient_gens::Vector{<:MPolyRingElem}                # generators chosen in sequence of max.contact steps
  ambient_orders::Vector{Int}                          # maximal orders at which hypersurfaces have been chosen
  dependent_vars::Vector{Int}                          # variables not used for system of parameters
  minor_data::Vector{<:Tuple{<:RingElem, Int}}           # chart is complement of product of these minors
                                                       #   (m,i) stands for the minor m and its size i
                                                       #   sum over sizes is lenght(dependent_vars)
  ambient_jacobi::Vector{<:MatrixElem}                 # jacobi-matrix*pseudoinverse for reduction

  function MaxContactChart(U::AbsAffineScheme, ambient_gens::Vector{<:MPolyRingElem},
                         ambient_orders::Vector{Int},
                         dependent_vars::Vector{Int},
                         minor_data::Vector{<:Tuple{<:RingElem, Int}},
                         ambient_jacob::Vector{<:MatrixElem})
    return new{base_ring_type(U)}(U,ambient_gens,ambient_orders,dependent_vars,minor_data,ambient_jacob)
  end
end

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
  max_contact_data::IdDict{<:AbsAffineScheme,MaxContactChart}

  function MaxContactObject(
               W::AbsCoveredScheme,
               Cov::Covering,
               max_contact_data::IdDict{AbsAffineScheme,MaxContactChart})
    return new{base_ring_type(W)}(W,Cov,max_contact_data)
  end
end