### Iteration over monomials in modules of a certain degree
struct AllModuleMonomials{ModuleType<:FreeMod}
  F::ModuleType
  d::Int

  function AllModuleMonomials(F::FreeMod{T}, d::Int) where {T <: MPolyDecRingElem}
    is_graded(F) || error("module must be graded")
    S = base_ring(F)
    is_standard_graded(S) || error("iterator implemented only for the standard graded case")
    return new{typeof(F)}(F, d)
  end
end

underlying_module(amm::AllModuleMonomials) = amm.F
degree(amm::AllModuleMonomials) = amm.d

function all_monomials(F::FreeMod{T}, d::Int) where {T<:MPolyDecRingElem}
  return AllModuleMonomials(F, d)
end

Base.eltype(amm::AllModuleMonomials{T}) where {T} = elem_type(T)

function Base.length(amm::AllModuleMonomials)
  F = underlying_module(amm)
  r = rank(F)
  R = base_ring(F)
  d = degree(amm)
  result = 0
  for i in 1:r
    d_loc = d - Int(degree(F[i])[1])
    d_loc < 0 && continue
    result = result + length(monomials_of_degree(R, d_loc))
  end
  return result
end
  
function Base.iterate(amm::AllModuleMonomials, state::Nothing = nothing)
  i = 1
  F = underlying_module(amm)
  d = degree(amm)
  R = base_ring(F)

  i = findfirst(i -> d - Int(degree(F[i])[1]) >= 0, 1:ngens(F))
  i === nothing && return nothing
  d_loc = d - Int(degree(F[i])[1])

  mon_it = monomials_of_degree(R, d_loc)
  res = iterate(mon_it, nothing)
  res === nothing && i == ngens(F) && return nothing

  x, s = res
  return x*F[1], (i, mon_it, s)
end

function Base.iterate(amm::AllModuleMonomials, state::Tuple{Int, AllMonomials, Vector{Int}})
  F = underlying_module(amm)
  d = degree(amm)
  R = base_ring(F)

  i, mon_it, s = state
  res = iterate(mon_it, s)
  if res === nothing
    i = findnext(i -> d - Int(degree(F[i])[1]) >= 0, 1:ngens(F), i + 1)
    i === nothing && return nothing
    d_loc = d - Int(degree(F[i])[1])

    mon_it = monomials_of_degree(R, d_loc)
    res_loc = iterate(mon_it, nothing)
    res_loc === nothing && i == ngens(F) && return nothing

    x, s = res_loc
    return x*F[i], (i, mon_it, s)
  end

  x, s = res
  return x*F[i], (i, mon_it, s)
end

