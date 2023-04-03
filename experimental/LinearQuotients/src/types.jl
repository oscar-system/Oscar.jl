mutable struct LinearQuotient{S, T}
  group::MatrixGroup{S, T}

  # We fix a primitive e-th root of unity, where e is the exponent of the group
  # for consistency
  root_of_unity::Tuple{S, Int}

  class_group::Tuple{GrpAbFinGen, Generic.CompositeMap{MatrixGroup{S, T}, GrpAbFinGen}}

  function LinearQuotient(G::MatrixGroup{S, T}) where {S, T}
    L = new{S, T}()
    L.group = G
    return L
  end
end
