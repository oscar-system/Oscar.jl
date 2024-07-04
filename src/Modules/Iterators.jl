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
    d_loc = d - Int(degree(F[i]; check=false)[1])
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

  i = findfirst(i -> d - Int(degree(F[i]; check=false)[1]) >= 0, 1:ngens(F))
  i === nothing && return nothing
  d_loc = d - Int(degree(F[i]; check=false)[1])

  mon_it = monomials_of_degree(R, d_loc)
  res = iterate(mon_it, nothing)
  res === nothing && i == ngens(F) && return nothing

  x, s = res
  return x*F[i], (i, mon_it, s)
end

function Base.iterate(amm::AllModuleMonomials, state::Tuple{Int, AllMonomials, Vector{Int}})
  F = underlying_module(amm)
  d = degree(amm)
  R = base_ring(F)

  i, mon_it, s = state
  res = iterate(mon_it, s)
  if res === nothing
    i = findnext(i -> d - Int(degree(F[i]; check=false)[1]) >= 0, 1:ngens(F), i + 1)
    i === nothing && return nothing
    d_loc = d - Int(degree(F[i]; check=false)[1])

    mon_it = monomials_of_degree(R, d_loc)
    res_loc = iterate(mon_it, nothing)
    res_loc === nothing && i == ngens(F) && return nothing

    x, s = res_loc
    return x*F[i], (i, mon_it, s)
  end

  x, s = res
  return x*F[i], (i, mon_it, s)
end

### The same for module exponents
struct AllModuleExponents{ModuleType<:FreeMod}
  F::ModuleType
  d::Int

  function AllModuleExponents(F::FreeMod{T}, d::Int) where {T <: MPolyDecRingElem}
    is_graded(F) || error("module must be graded")
    S = base_ring(F)
    is_standard_graded(S) || error("iterator implemented only for the standard graded case")
    return new{typeof(F)}(F, d)
  end
end

underlying_module(amm::AllModuleExponents) = amm.F
degree(amm::AllModuleExponents) = amm.d

function all_exponents(F::FreeMod{T}, d::Int) where {T<:MPolyDecRingElem}
  return AllModuleExponents(F, d)
end

Base.eltype(amm::AllModuleExponents{T}) where {T} = Tuple{Vector{Int}, Int}

function Base.length(amm::AllModuleExponents)
  F = underlying_module(amm)
  r = rank(F)
  R = base_ring(F)
  n = ngens(R)
  d = degree(amm)
  result = 0
  for i in 1:r
    d_loc = d - Int(degree(F[i]; check=false)[1])
    d_loc < 0 && continue
    result = result + length(MultiIndicesOfDegree(n, d_loc))
  end
  return result
end
  
function Base.iterate(amm::AllModuleExponents, state::Nothing = nothing)
  i = 1
  F = underlying_module(amm)
  d = degree(amm)
  R = base_ring(F)
  n = ngens(R)

  i = findfirst(i -> d - Int(degree(F[i]; check=false)[1]) >= 0, 1:ngens(F))
  i === nothing && return nothing
  d_loc = d - Int(degree(F[i]; check=false)[1])

  exp_it = MultiIndicesOfDegree(n, d_loc)
  res = iterate(exp_it, nothing)
  res === nothing && i == ngens(F) && return nothing

  e, _ = res
  return (e, i), (i, exp_it, e)
end

function Base.iterate(amm::AllModuleExponents, state::Tuple{Int, MultiIndicesOfDegree, Vector{Int}})
  F = underlying_module(amm)
  d = degree(amm)
  R = base_ring(F)
  n = ngens(R)

  i, exp_it, e = state
  res = iterate(exp_it, e)
  if res === nothing
    i = findnext(i -> d - Int(degree(F[i]; check=false)[1]) >= 0, 1:ngens(F), i + 1)
    i === nothing && return nothing
    d_loc = d - Int(degree(F[i]; check=false)[1])

    exp_it = MultiIndicesOfDegree(n, d_loc)
    res_loc = iterate(exp_it, nothing)
    res_loc === nothing && i == ngens(F) && return nothing

    e, _ = res_loc
    return (e, i), (i, exp_it, e)
  end

  e, _ = res
  return (e, i), (i, exp_it, e)
end

### Iteration over monomials in Subquos
#
#= Disabled for the moment and continued soon (10.1.2024)
#
# We need a Groebner basis over a field for this to work; hence the 
# type restriction.
mutable struct AllSubquoMonomials{ModuleType<:SubquoModule{<:MPolyDecRingElem{<:FieldElem}}}
  M::ModuleType
  d::Int
  pres::ComplexOfMorphisms # A field for the presentation of M
  lead::SubModuleOfFreeModule # The leading module of the modulus in the presentation

  function AllSubquoMonomials(M::SubquoModule{T}, d::Int) where {U<:FieldElem, T <: MPolyDecRingElem{U}}
    is_graded(M) || error("module must be graded")
    S = base_ring(M)
    is_standard_graded(S) || error("iterator implemented only for the standard graded case")
    return new{typeof(M)}(M, d)
  end
end

underlying_module(amm::AllSubquoMonomials) = amm.M
degree(amm::AllSubquoMonomials) = amm.d

function presentation(amm::AllSubquoMonomials)
  !isdefined(amm, :pres) && (amm.pres = presentation(underlying_module(amm)))
  return amm.pres
end

function leading_module_of_presentation(amm::AllSubquoMonomials)
  if !isdefined(amm, :lead)
    I, _ = image(map(presentation(amm), 1))
    amm.lead = leading_module(I.sub)
  end
  return amm.lead
end

function all_monomials(M::SubquoModule{T}, d::Int) where {U<:FieldElem, T <: MPolyDecRingElem{U}}
  return AllModuleMonomials(M, d)
end

Base.eltype(amm::AllSubquoMonomials{T}) where {T} = elem_type(T)

function Base.length(amm::AllSubquoMonomials)
  p = presentation(amm)
  F = p[0]::FreeMod
  result = 0
  I = leading_module_of_presentation(amm)
  for x in all_monomials(F, degree(amm))
    x in I && continue
    result += 1
  end
  return result
end
  
function Base.iterate(amm::AllSubquoMonomials, state::Nothing = nothing)
  F = presentation(amm)[0]
  mon_it = all_monomials(F, degree(amm))
  I = leading_module_of_presentation(amm)
  next = zero(underlying_module(amm))
  res = iterate(mon_it)
  res === nothing && return nothing
  m, s = res
  while m in I
    res = iterate(mon_it, s)
    res === nothing && return nothing
    m, s = res
  end
  next = map(presentation(amm), 0)(m)
  return next, (mon_it, s)
end

function Base.iterate(amm::AllSubquoMonomials, state::Tuple)
  mon_it, s = state
  res = iterate(mon_it, s)
  res === nothing && return nothing
  m, s = res
  I = leading_module_of_presentation(amm)
  while m in I
    res = iterate(mon_it, s)
    res === nothing && return nothing
    m, s = res
  end
  next = map(presentation(amm), 0)(m)
  return next, (mon_it, s)
end

### Iteration over exponent vectors of SubquoModules
mutable struct AllSubquoExponents{ModuleType<:SubquoModule{<:MPolyDecRingElem{<:FieldElem}}}
  M::ModuleType
  d::Int
  pres::ComplexOfMorphisms # A field for the presentation of M
  lead::SubModuleOfFreeModule # The leading module of the modulus in the presentation

  function AllSubquoExponents(M::SubquoModule{T}, d::Int) where {U<:FieldElem, T <: MPolyDecRingElem{U}}
    is_graded(M) || error("module must be graded")
    S = base_ring(M)
    is_standard_graded(S) || error("iterator implemented only for the standard graded case")
    return new{typeof(M)}(M, d)
  end
end

underlying_module(amm::AllSubquoExponents) = amm.M
degree(amm::AllSubquoExponents) = amm.d

function presentation(amm::AllSubquoExponents)
  !isdefined(amm, :pres) && (amm.pres = presentation(underlying_module(amm)))
  return amm.pres
end

function leading_module_of_presentation(amm::AllSubquoExponents)
  if !isdefined(amm, :lead)
    I, _ = image(map(presentation(amm), 1))
    amm.lead = leading_module(I.sub)
  end
  return amm.lead
end

function all_exponents(M::SubquoModule{T}, d::Int) where {U<:FieldElem, T <: MPolyDecRingElem{U}}
  return AllSubquoExponents(M, d)
end

Base.eltype(amm::AllSubquoExponents) = Tuple{Vector{Int}, Int}

function Base.length(amm::AllSubquoExponents)
  p = presentation(amm)
  F = p[0]::FreeMod
  R = base_ring(F)
  v = gens(R)
  result = 0
  I = leading_module_of_presentation(amm)
  for e in all_exponents(F, degree(amm))
    m = prod(x^k for (x, k) in zip(v, e[1]); init=one(R))*F[e[2]]
    m in I && continue
    result += 1
  end
  return result
end
  
function Base.iterate(amm::AllSubquoExponents, state::Nothing = nothing)
  F = presentation(amm)[0]
  R = base_ring(F)
  v = gens(R)
  mon_it = all_exponents(F, degree(amm))
  I = leading_module_of_presentation(amm)
  next = zero(underlying_module(amm))
  res = iterate(mon_it)
  res === nothing && return nothing
  e, s = res
  m = prod(x^k for (x, k) in zip(v, e[1]); init=one(R))*F[e[2]]
  while m in I
    res = iterate(mon_it, s)
    res === nothing && return nothing
    e, s = res
    m = prod(x^k for (x, k) in zip(v, e[1]); init=one(R))*F[e[2]]
  end
  return e, (mon_it, s)
end

function Base.iterate(amm::AllSubquoExponents, state::Tuple)
  R = base_ring(underlying_module(amm))
  v = gens(R)
  mon_it, s = state
  res = iterate(mon_it, s)
  res === nothing && return nothing
  e, s = res
  F = presentation(amm)[0]
  m = prod(x^k for (x, k) in zip(v, e[1]); init=one(R))*F[e[2]]
  I = leading_module_of_presentation(amm)
  while m in I
    res = iterate(mon_it, s)
    res === nothing && return nothing
    e, s = res
    m = prod(x^k for (x, k) in zip(v, e[1]); init=one(R))*F[e[2]]
  end
  return e, (mon_it, s)
end
=#
