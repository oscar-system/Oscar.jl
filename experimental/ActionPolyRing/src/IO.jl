#######################################
#
# String I/O
#
#######################################

### Difference ###
function Base.show(io::IO, ::MIME"text/plain", dpr::DifferencePolyRing)
  io = pretty(io)
  n = nelemvars(dpr)
  print(io, "Difference polynomial ring in $n elementary variables ")
  for i in 1:n-1
    print(io, string(symbols(dpr)[i]) * ", ")
  end
  print(io, string(symbols(dpr)[n])*"\n")
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
    print(terse(io), "Difference polynomial ring in $(nelemvars(dpr)) elementary variables over ")
    print(terse(io), base_ring(dpr))
  end
end

