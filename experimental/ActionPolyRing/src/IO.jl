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
  print(io, "with $(ndiffs(dpr)) commuting endomorphisms\n")
  print(io, Indent())
  print(io, "over ", Lowercase(), base_ring(dpr))
  print(io, Dedent())
end

function Base.show(io::IO, dpr::DifferencePolyRing)
  io = pretty(io)
  if is_terse(io)
    print(io, "Difference polynomial ring")
  else
    print(io, "Difference polynomial ring in $(n_elementary_symbols(dpr)) elementary symbols over ")
    print(terse(io), Lowercase(), base_ring(dpr))
  end
end

### Difference ###
function Base.show(io::IO, ::MIME"text/plain", dpr::DifferentialPolyRing)
  io = pretty(io)
  n = n_elementary_symbols(dpr)
  print(io, "Differential polynomial ring in $n elementary symbols ")
  join(io, elementary_symbols(dpr), ", ")
  print(io, "\n")
  print(io, "with $(ndiffs(dpr)) commuting derivations\n")
  print(io, Indent())
  print(io, "over ", Lowercase(), base_ring(dpr))
  print(io, Dedent())
end

function Base.show(io::IO, dpr::DifferentialPolyRing)
  io = pretty(io)
  if is_terse(io)
    print(io, "Differential polynomial ring")
  else
    print(io, "Differential polynomial ring in $(n_elementary_symbols(dpr)) elementary symbols over ")
    print(terse(io), Lowercase(), base_ring(dpr))
  end
end

###############################################################################
#
#  Expressify
#
###############################################################################

@enable_all_show_via_expressify ActionPolyRingElem

function _expressify_monomial!(prod::Expr, x, e)
  @inbounds for i in 1:length(e)
    if e[i] > 1
      push!(prod.args, Expr(:call, :^, x[i], e[i]))  # deepcopy not needed for immutables
    elseif e[i] == 1
      push!(prod.args, x[i])
    end
  end
end

function expressify(a::ActionPolyRingElem, x = symbols(parent(a)); context = nothing)
  sum = Expr(:call, :+)

  for (c, e) in zip(coefficients(a), exponents(a))
    if all(zero(eltype(e)) == ei for ei in e)  # constant term
      push!(sum.args, expressify(c, context = context))
    else
      prod = Expr(:call, :*)
      push!(prod.args, expressify(c, context = context))
      _expressify_monomial!(prod, x, e)
      push!(sum.args, prod)
    end
  end

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
  end
  print(io, "Coefficients iterator of ")
  print(terse(io), ci.poly)
end

### Exponents ###
function Base.show(io::IO, ei::ActionPolyExponentVectors)
  io = pretty(io)
  if is_terse(io)
    print(io, "Exponents iterator")
  end
  print(io, "Exponents iterator of ")
  print(terse(io), ei.poly)
end

### Monomials ###
function Base.show(io::IO, mi::ActionPolyMonomials)
  io = pretty(io)
  if is_terse(io)
    print(io, "Monomials iterator")
  end
  print(io, "Monomials iterator of ")
  print(terse(io), mi.poly)
end

### Terms ###
function Base.show(io::IO, ti::ActionPolyTerms)
  io = pretty(io)
  if is_terse(io)
    print(io, "Terms iterator")
  end
  print(io, "Terms iterator of ")
  print(terse(io), ti.poly)
end

###############################################################################
#
#  Rankings 
#
###############################################################################

### Difference ###
function Base.show(io::IO, ::MIME"text/plain", ran::ActionPolyRingRanking)
  io = pretty(io)
  print(io, "Ranking of ", Lowercase(), base_ring(ran))
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
    if base_ring(ran) isa DifferencePolyRing
      print(terse(io), "Ranking of difference polynomial ring")
    end
    if base_ring(ran) isa DifferentialPolyRing
      print(terse(io), "Ranking of differential polynomial ring")
    end
  end
end

