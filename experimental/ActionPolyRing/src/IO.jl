###############################################################################
#
#  Action polynomial rings
#
###############################################################################

### Difference ###
function Base.show(io::IO, ::MIME"text/plain", dpr::DifferencePolyRing)
  io = pretty(io)
  n = n_elementary_symbols(dpr)
  print(io, "Difference polynomial ring in $n elementary symbols ")
  join(io, elementary_symbols(dpr), ", ")
  print(io, "\n")
  print(io, "with $(n_action_maps(dpr)) commuting endomorphisms\n")
  print(io, Indent())
  print(io, "over ", Lowercase(), coefficient_ring(dpr))
  print(io, Dedent())
end

function Base.show(io::IO, dpr::DifferencePolyRing)
  io = pretty(io)
  if is_terse(io)
    print(io, "Difference polynomial ring")
  else
    print(io, "Difference polynomial ring in $(n_elementary_symbols(dpr)) elementary symbols over ")
    print(terse(io), Lowercase(), coefficient_ring(dpr))
  end
end

### Difference ###
function Base.show(io::IO, ::MIME"text/plain", dpr::DifferentialPolyRing)
  io = pretty(io)
  n = n_elementary_symbols(dpr)
  print(io, "Differential polynomial ring in $n elementary symbols ")
  join(io, elementary_symbols(dpr), ", ")
  print(io, "\n")
  print(io, "with $(n_action_maps(dpr)) commuting derivations\n")
  print(io, Indent())
  print(io, "over ", Lowercase(), coefficient_ring(dpr))
  print(io, Dedent())
end

function Base.show(io::IO, dpr::DifferentialPolyRing)
  io = pretty(io)
  if is_terse(io)
    print(io, "Differential polynomial ring")
  else
    print(io, "Differential polynomial ring in $(n_elementary_symbols(dpr)) elementary symbols over ")
    print(terse(io), Lowercase(), coefficient_ring(dpr))
  end
end

###############################################################################
#
#  Expressify
#
###############################################################################

@enable_all_show_via_expressify ActionPolyRingElem

function __expressify_coeff_monomial!(coeff_prod::Expr, x, e, ld_ind)
  @inbounds for i in (ld_ind + 1):length(e)
    if e[i] > 1 
      push!(coeff_prod.args, Expr(:call, :^, x[i], e[i]))
    elseif e[i] == 1
      push!(coeff_prod.args, x[i])
    end
  end
end

function expressify(a::ActionPolyRingElem, x = symbols(parent(a)); context = nothing)
  es = exponents(a)
  if length(es) == 0
    return Expr(:call, :+)
  end

  # Find index and highest exponent of the leader of a
  ld_ind, cur_exp = (0, 0)
  for (i, v) in enumerate(first(es))
    if v != 0
      ld_ind, cur_exp = i, v
      break
    end
  end

  if ld_ind == 0 # a is an element of the base ring
    return Expr(:call, :+, expressify(coeff(a, 1), context = context))  
  end

  sum = Expr(:call, :+)
  ld = x[ld_ind]
  prev_exp = cur_exp
  coeff_sum = Expr(:call, :+)
  for (c, e) in zip(coefficients(a), es)
    cur_exp = e[ld_ind]
    coeff_prod = Expr(:call, :*) 
    push!(coeff_prod.args, expressify(c, context = context))
    __expressify_coeff_monomial!(coeff_prod, x, e, ld_ind)

    if cur_exp == prev_exp # still part of the coefficient sum
      push!(coeff_sum.args, coeff_prod)
      continue
    end
      
    # At this point, coeff_prod is part of the new coeff_sum, so we push the old (with the old exponent) and reset afterwards
    coeff_ld_exp_prod = Expr(:call, :*, coeff_sum) # Multiply coefficient term and exponent of the leader
    if prev_exp > 1
      push!(coeff_ld_exp_prod.args, Expr(:call, :^, ld, prev_exp))
    elseif prev_exp == 1
      push!(coeff_ld_exp_prod.args, ld)
    end
    push!(sum.args, coeff_ld_exp_prod)

    # Reset
    prev_exp = cur_exp
    coeff_sum = Expr(:call, :+, coeff_prod)
  end 

  # At this point, the last summand is not yet pushed, so we do it manually
  coeff_ld_exp_prod = Expr(:call, :*, coeff_sum) # Multiply coefficient term and exponent of the leader
  if cur_exp > 1
    push!(coeff_ld_exp_prod.args, Expr(:call, :^, ld, cur_exp))
  elseif cur_exp == 1
    push!(coeff_ld_exp_prod.args, ld)
  end
  push!(sum.args, coeff_ld_exp_prod)

  return sum 
end

###############################################################################
#
#  Iterators 
#
###############################################################################

### Coefficients ###
function Base.show(io::IO, ci::ActionPolyCoeffs)
  io = pretty(io)
  if is_terse(io)
    print(io, "Coefficients iterator")
  else
    print(io, "Coefficients iterator of ")
    print(terse(io), ci.poly)
  end
end

### Exponents ###
function Base.show(io::IO, ei::ActionPolyExponentVectors)
  io = pretty(io)
  if is_terse(io)
    print(io, "Exponents iterator")
  else
    print(io, "Exponents iterator of ")
    print(terse(io), ei.poly)
  end
end

### Monomials ###
function Base.show(io::IO, mi::ActionPolyMonomials)
  io = pretty(io)
  if is_terse(io)
    print(io, "Monomials iterator")
  else
    print(io, "Monomials iterator of ")
    print(terse(io), mi.poly)
  end
end

### Terms ###
function Base.show(io::IO, ti::ActionPolyTerms)
  io = pretty(io)
  if is_terse(io)
    print(io, "Terms iterator")
  else
    print(io, "Terms iterator of ")
    print(terse(io), ti.poly)
  end
end

###############################################################################
#
#  Rankings 
#
###############################################################################

### Difference ###
function Base.show(io::IO, ::MIME"text/plain", ran::ActionPolyRingRanking)
  io = pretty(io)
  print(io, "Ranking of ", Lowercase(), parent(ran))
  print(io, "\n")
  print(io, "with elementary symbols partitioned by\n")
  print(io, Indent())
  print(io, partition(ran))
  print(io, Dedent())
  print(io, "\nand ordering of the indices defined by\n")
  print(io, Indent())
  show(io, "text/plain", index_ordering_matrix(ran))
  print(io, Dedent())
end

function Base.show(io::IO, ran::ActionPolyRingRanking)
  io = pretty(io)
  if is_terse(io)
    print(io, "Ranking")
  else
    if parent(ran) isa DifferencePolyRing
      print(terse(io), "Ranking of difference polynomial ring")
    end
    if parent(ran) isa DifferentialPolyRing
      print(terse(io), "Ranking of differential polynomial ring")
    end
  end
end

