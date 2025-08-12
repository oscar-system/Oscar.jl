###############################################################################
#
#  Action polynomial rings
#
###############################################################################

### Difference ###
function Base.show(io::IO, ::MIME"text/plain", dpr::DifferencePolyRing)
  io = pretty(io)
  n = nelementary_symbols(dpr)
  print(io, "Difference polynomial ring in $n elementary symbols ")
  for i in 1:n-1
    print(io, string(elementary_symbols(dpr)[i]) * ", ")
  end
  print(io, string(elementary_symbols(dpr)[n])*"\n")
  print(io, "with $(ndiffs(dpr)) commuting endomorphisms\n")
  print(io, Indent())
  print(io, "over $(base_ring(dpr))")
  print(io, Dedent())
end

function Base.show(io::IO, dpr::DifferencePolyRing)
  io = pretty(io)
  if is_terse(io)
    print(io, "Difference polynomial ring")
  else
    print(terse(io), "Difference polynomial ring in $(nelementary_symbols(dpr)) elementary symbols over ")
    print(terse(io), base_ring(dpr))
  end
end

### Difference ###
function Base.show(io::IO, ::MIME"text/plain", dpr::DifferentialPolyRing)
  io = pretty(io)
  n = nelementary_symbols(dpr)
  print(io, "Differential polynomial ring in $n elementary symbols ")
  for i in 1:n-1
    print(io, string(elementary_symbols(dpr)[i]) * ", ")
  end
  print(io, string(elementary_symbols(dpr)[n])*"\n")
  print(io, "with $(ndiffs(dpr)) commuting derivations\n")
  print(io, Indent())
  print(io, "over $(base_ring(dpr))")
  print(io, Dedent())
end

function Base.show(io::IO, dpr::DifferentialPolyRing)
  io = pretty(io)
  if is_terse(io)
    print(io, "Differential polynomial ring")
  else
    print(terse(io), "Differential polynomial ring in $(nelementary_symbols(dpr)) elementary symbols over ")
    print(terse(io), base_ring(dpr))
  end
end

###############################################################################
#
# Expressify
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
            if !(isone(c))
                push!(prod.args, expressify(c, context = context))
            end
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
function Base.show(io::IO, ::MIME"text/plain", ci::Union{DifferencePolyCoeffs, DifferentialPolyCoeffs})
  io = pretty(io)
  print(io, "coefficients iterator of ")
  show(io, ci.poly)
end

function Base.show(io::IO, ci::Union{DifferencePolyCoeffs, DifferentialPolyCoeffs})
  io = pretty(io)
  if is_terse(io)
    print(io, "iterator")
  end
  print(terse(io), "coefficients iterator")
end

### Exponents ###
function Base.show(io::IO, ::MIME"text/plain", ei::Union{DifferencePolyExponentVectors, DifferentialPolyExponentVectors})
  io = pretty(io)
  print(io, "exponents iterator of ")
  show(io, ei.poly)
end

function Base.show(io::IO, ei::Union{DifferencePolyExponentVectors, DifferentialPolyExponentVectors})
  io = pretty(io)
  if is_terse(io)
    print(io, "iterator")
  end
  print(terse(io), "exponents iterator")
end

### Monomials ###
function Base.show(io::IO, ::MIME"text/plain", mi::Union{DifferencePolyMonomials, DifferentialPolyMonomials})
  io = pretty(io)
  print(io, "monomials iterator of ")
  show(io, mi.poly)
end

function Base.show(io::IO, mi::Union{DifferencePolyMonomials, DifferentialPolyMonomials})
  io = pretty(io)
  if is_terse(io)
    print(io, "iterator")
  end
  print(terse(io), "monomials iterator")
end

### Terms ###
function Base.show(io::IO, ::MIME"text/plain", ti::Union{DifferencePolyTerms, DifferentialPolyTerms})
  io = pretty(io)
  print(io, "terms iterator of ")
  show(io, ti.poly)
end

function Base.show(io::IO, ti::Union{DifferencePolyTerms, DifferentialPolyTerms})
  io = pretty(io)
  if is_terse(io)
    print(io, "iterator")
  end
  print(terse(io), "coefficients iterator")
end

###############################################################################
#
#  Rankings 
#
###############################################################################

### Difference ###
function Base.show(io::IO, ::MIME"text/plain", ran::DifferenceRanking)
  io = pretty(io)
  print(io, "Ranking of $(base_ring(ran))\n")
  print(io, "with elementary symbols partitioned by\n")
  print(io, Indent())
  print(io, partition(ran))
  print(io, Dedent())
  print(io, "\nand ordering of the indices defined by\n")
  print(io, Indent())
  show(io, "text/plain", index_ordering_matrix(ran))
  print(io, Dedent())
end

function Base.show(io::IO, ran::DifferenceRanking)
  io = pretty(io)
  if is_terse(io)
    print(io, "Ranking")
  end
  print(terse(io), "Ranking of difference polynomial ring")
end

### Differential ###
function Base.show(io::IO, ::MIME"text/plain", ran::DifferentialRanking)
  io = pretty(io)
  print(io, "Ranking of $(base_ring(ran))\n")
  print(io, "with elementary symbols partitioned by\n")
  print(io, Indent())
  print(io, partition(ran))
  print(io, Dedent())
  print(io, "\nand ordering of the indices defined by\n")
  print(io, Indent())
  show(io, "text/plain", index_ordering_matrix(ran))
  print(io, Dedent())
end

function Base.show(io::IO, ran::DifferentialRanking)
  io = pretty(io)
  if is_terse(io)
    print(io, "Ranking")
  end
  print(terse(io), "Ranking of differential polynomial ring")
end

