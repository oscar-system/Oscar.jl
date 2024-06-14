###############################################################################
# A place to accumulate code that should eventually be moved to Hecke.jl
###############################################################################

function getindex(r::Hecke.SRow, u::AbstractUnitRange)
  s = sparse_row(base_ring(r))
  shift = 1-first(u)
  for (p,v) = r
    if p in u
      push!(s.pos, p+shift)
      push!(s.values, v)
    end
  end
  return s
end

Oscar.canonical_unit(x::AbsSimpleNumFieldOrderQuoRingElem) = one(parent(x))
