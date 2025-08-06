#######################################
#
# String I/O
#
#######################################

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

function Base.show(io::IO, apre::ActionPolyRingElem)
  io = pretty(io)
  show(io, data(apre))
end

### Rankings ###
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
    print(io, "Ranking of difference polynomial ring")
  end
  print(terse(io), "Ranking of difference polynomial ring")
end

