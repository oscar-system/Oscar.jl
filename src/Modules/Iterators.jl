### Iteration over monomials in modules of a certain degree
struct AllModuleMonomials{ModuleType<:FreeMod, DegreeType<:Union{Int, FinGenAbGroupElem}}
  F::ModuleType
  d::DegreeType
  # Cache for the exponent vectors of the `base_ring` of the module in degree `d`
  # For multigradings we're using polymake and this doesn't return iterators. 
  # Since the computation is expensive, we don't want to do it over and over again. 
  exp_cache::Dict{FinGenAbGroupElem, Vector{Vector{Int}}}

  function AllModuleMonomials(F::FreeMod{T}, d::Int) where {T <: MPolyDecRingElem}
    is_graded(F) || error("module must be graded")
    S = base_ring(F)
    is_standard_graded(S) || error("iterator implemented only for the standard graded case")
    return new{typeof(F), Int}(F, d)
  end
  function AllModuleMonomials(F::FreeMod{T}, d::FinGenAbGroupElem) where {T <: MPolyDecRingElem}
    is_graded(F) || error("module must be graded")
    S = base_ring(F)
    is_zm_graded(S) || error("iterator implemented only for the ZZ^m-graded case")
    return new{typeof(F), FinGenAbGroupElem}(F, d, Dict{FinGenAbGroupElem, Vector{Vector{Int}}}())
  end
end

underlying_module(amm::AllModuleMonomials) = amm.F
degree(amm::AllModuleMonomials) = amm.d

function all_monomials(F::FreeMod{T}, d::Union{Int, FinGenAbGroupElem}) where {T<:MPolyDecRingElem}
  return AllModuleMonomials(F, d)
end

Base.eltype(::Type{<:AllModuleMonomials{T}}) where {T} = elem_type(T)

function Base.length(amm::AllModuleMonomials{<:FreeMod, Int})
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
  
function Base.length(amm::AllModuleMonomials{<:FreeMod, FinGenAbGroupElem})
  F = underlying_module(amm)
  r = rank(F)
  R = base_ring(F)
  d = degree(amm)::FinGenAbGroupElem
  return sum(length(get!(amm.exp_cache, d - degree(g; check=false), all_exponents(R, d - degree(g)))) for g in gens(F); init=0)
end
  
function Base.iterate(amm::AllModuleMonomials{<:FreeMod, Int}, state::Nothing = nothing)
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

function Base.iterate(
    amm::AllModuleMonomials{<:FreeMod, Int}, 
    state::Tuple{Int, AllMonomials, Vector{Int}}
  )
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

function Base.iterate(amm::AllModuleMonomials{<:FreeMod, FinGenAbGroupElem}, state::Nothing = nothing)
  F = underlying_module(amm)
  d = degree(amm)::FinGenAbGroupElem
  R = base_ring(F)

  for (i, g) in enumerate(gens(F))
    d_loc = d - degree(g; check=false)

    exps = get!(amm.exp_cache, d_loc) do
      all_exponents(R, d_loc)
    end

    res = iterate(exps)
    res === nothing && continue

    e, s = res
    poly_ctx = MPolyBuildCtx(R)
    push_term!(poly_ctx, one(coefficient_ring(R)), e)
    return (finish(poly_ctx)*g)::elem_type(F), (i, exps, s)
  end
  return nothing
end

function Base.iterate(
    amm::AllModuleMonomials{<:FreeMod, FinGenAbGroupElem}, 
    state::Tuple{Int, Vector{Vector{Int}}, Int}
  )
  F = underlying_module(amm)
  d = degree(amm)
  R = base_ring(F)

  i, exps, s = state
  res = iterate(exps, s)
  if res === nothing
    for j in i+1:ngens(F)
      d_loc = d - degree(F[j]; check=false)
      exps = get!(amm.exp_cache, d_loc) do
        all_exponents(R, d_loc)
      end
      res_loc = iterate(exps)
      res_loc === nothing && continue
      e, s = res_loc
      poly_ctx = MPolyBuildCtx(R)
      push_term!(poly_ctx, one(coefficient_ring(R)), e)
      return finish(poly_ctx)*F[j], (j, exps, s)
    end
    return nothing
  end
  e, s = res
  poly_ctx = MPolyBuildCtx(R)
  push_term!(poly_ctx, one(coefficient_ring(R)), e)
  return (finish(poly_ctx)*F[i])::elem_type(F), (i, exps, s)
end

### The same for module exponents
struct AllModuleExponents{ModuleType<:FreeMod, DegreeType<:Union{Int, FinGenAbGroupElem}}
  F::ModuleType
  d::DegreeType
  exp_cache::Dict{FinGenAbGroupElem, Vector{Vector{Int}}}
  check::Bool

  function AllModuleExponents(F::FreeMod{T}, d::Int; check::Bool=true) where {T <: MPolyDecRingElem}
    is_graded(F) || error("module must be graded")
    S = base_ring(F)
    is_standard_graded(S) || error("iterator implemented only for the standard graded case")
    return new{typeof(F), Int}(F, d, check)
  end
  function AllModuleExponents(F::FreeMod{T}, d::FinGenAbGroupElem; check::Bool=true) where {T <: MPolyDecRingElem}
    is_graded(F) || error("module must be graded")
    S = base_ring(F)
    is_zm_graded(S) || error("iterator implemented only for the ZZ^m-graded case")
    return new{typeof(F), FinGenAbGroupElem}(F, d, Dict{FinGenAbGroupElem, Vector{Vector{Int}}}(), check)
  end
end

underlying_module(amm::AllModuleExponents) = amm.F
degree(amm::AllModuleExponents) = amm.d

function all_exponents(F::FreeMod{T}, d::Union{Int, FinGenAbGroupElem}; check::Bool=true) where {T<:MPolyDecRingElem}
  return AllModuleExponents(F, d; check)
end

Base.eltype(::Type{<:AllModuleExponents{T}}) where {T} = Tuple{Vector{Int}, Int}

function Base.length(amm::AllModuleExponents{<:FreeMod, Int})
  F = underlying_module(amm)
  r = rank(F)
  R = base_ring(F)
  n = ngens(R)
  d = degree(amm)
  result = 0
  for i in 1:r
    d_loc = d - Int(degree(F[i]; check=false)[1])
    d_loc < 0 && continue
    result = result + Int(number_of_weak_compositions(d_loc, n))
  end
  return result
end
  
function Base.iterate(amm::AllModuleExponents{<:FreeMod, Int}, state::Nothing = nothing)
  i = 1
  F = underlying_module(amm)
  d = degree(amm)
  R = base_ring(F)
  n = ngens(R)

  i = findfirst(i -> d - Int(degree(F[i]; check=false)[1]) >= 0, 1:ngens(F))
  i === nothing && return nothing
  d_loc = d - Int(degree(F[i]; check=false)[1])

  exp_it = weak_compositions(d_loc, n)
  res = iterate(exp_it, nothing)
  res === nothing && i == ngens(F) && return nothing

  e, s = res
  return (data(e), i), (i, exp_it, s)
end

function Base.iterate(amm::AllModuleExponents{<:FreeMod, Int}, state::Tuple{Int, WeakCompositions{Int}, Vector{Int}})
  F = underlying_module(amm)
  d = degree(amm)
  R = base_ring(F)
  n = ngens(R)

  i, exp_it, s = state
  res = iterate(exp_it, s)
  if res === nothing
    i = findnext(i -> d - Int(degree(F[i]; check=false)[1]) >= 0, 1:ngens(F), i + 1)
    i === nothing && return nothing
    d_loc = d - Int(degree(F[i]; check=false)[1])

    exp_it = weak_compositions(d_loc, n)
    res_loc = iterate(exp_it, nothing)
    res_loc === nothing && i == ngens(F) && return nothing

    e, s = res_loc
    return (data(e), i), (i, exp_it, s)
  end

  e, s = res
  return (data(e), i), (i, exp_it, s)
end

function Base.length(amm::AllModuleExponents{<:FreeMod, FinGenAbGroupElem})
  F = underlying_module(amm)
  r = rank(F)
  R = base_ring(F)
  d = degree(amm)::FinGenAbGroupElem
  return sum(length(get!(amm.exp_cache, d - degree(g; check=false), all_exponents(R, d - degree(g); check=amm.check))) for g in gens(F); init=0)
  result = 0
  for i in 1:r
    d_loc = d - degree(F[i]; check=false)
    any(Int(d_loc[i]) < 0 for i in 1:ngens(parent(d)))&& continue
    exps = get!(amm.exp_cache, d_loc) do
      return all_exponents(R, d_loc; check=amm.check)
    end
    result += length(exps)
  end
  return result
end
  
function Base.iterate(amm::AllModuleExponents{<:FreeMod, FinGenAbGroupElem}, state::Nothing = nothing)
  F = underlying_module(amm)
  d = degree(amm)::FinGenAbGroupElem
  R = base_ring(F)

  for (i, g) in enumerate(gens(F))
    d_loc = d - degree(g; check=false)

    exps = get!(amm.exp_cache, d_loc) do
      all_exponents(R, d_loc; check=amm.check)
    end

    res = iterate(exps)
    res === nothing && continue

    e, s = res
    return (e, i), (i, exps, s)
  end
  return nothing
end

function Base.iterate(
    amm::AllModuleExponents{<:FreeMod, FinGenAbGroupElem}, 
    state::Tuple{Int, Vector{Vector{Int}}, Int}
  )
  F = underlying_module(amm)
  d = degree(amm)
  R = base_ring(F)

  i, exps, s = state
  res = iterate(exps, s)
  if res === nothing
    for j in i+1:ngens(F)
      d_loc = d - degree(F[j]; check=false)
      exps = get!(amm.exp_cache, d_loc) do
        all_exponents(R, d_loc; check=amm.check)
      end
      res_loc = iterate(exps)
      res_loc === nothing && continue
      e, s = res_loc
      return (e, j), (j, exps, s)
    end
    return nothing
  end
  e, s = res
  return (e, i), (i, exps, s)
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
