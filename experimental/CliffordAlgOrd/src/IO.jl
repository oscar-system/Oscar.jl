##################################################
#
#  Clifford algebra and elements
#
##################################################
 
##### Algebra #####
function Base.show(io::IO, ::MIME"text/plain", C::CliffordAlgebra)
  io = pretty(io)
  print(io, LowercaseOff(), "Clifford algebra of quadratic space with Gram matrix\n")
  print(io, Indent())
  show(io, "text/plain", gram_matrix(C))
  print(io, Dedent(), "\ndefined over ", Lowercase(), base_ring(C))
end
 
function Base.show(io::IO, C::CliffordAlgebra)
  io = pretty(io)
  if is_terse(io)
    print(io, LowercaseOff(), "Clifford algebra")
  else
    print(terse(io), LowercaseOff(), "Clifford algebra over ")
    print(terse(io), Lowercase(), base_ring(C))
  end
end

##### Elements #####
function Base.show(io::IO, x::CliffordAlgebraElem)
  cf = coefficients(x)
  print(IOContext(io, :typeinfo => typeof(cf)), cf)
end

##################################################
#   
#  Clifford Order and elements
#
##################################################
  
##### Order #####
function Base.show(io::IO, ::MIME"text/plain", C::CliffordOrder)
  io = pretty(io)
  print(io, LowercaseOff(), "Clifford order of even lattice over ")
  print(io, Lowercase(), base_ring(C))
  print(io, " with Gram matrix\n")
  print(io, Indent())
  show(io, "text/plain", gram_matrix(C))
  print(io, Dedent(), "\nand coefficient ideals of the lattice\n")
  print(io, Indent())
  show(io, "text/plain", _coefficient_ideals_of_lattice(lattice(C)))
  print(io, Dedent())
end

function Base.show(io::IO, ::MIME"text/plain", C::ZZCliffordOrder)
  io = pretty(io)
  print(io, LowercaseOff(), "Clifford order of even integer lattice with Gram matrix\n")
  print(io, Indent())
  show(io, "text/plain", gram_matrix(C))
  print(io, Dedent())
end

function Base.show(io::IO, C::Union{ZZCliffordOrder, CliffordOrder})
  io = pretty(io)
  if is_terse(io)
    print(io, LowercaseOff(), "Clifford order")
  else
    print(terse(io), LowercaseOff(), "Clifford order over ")
    print(terse(io), Lowercase(), base_ring(C))
  end
end

##### Elements #####
function Base.show(io::IO, x::Union{ZZCliffordOrderElem, CliffordOrderElem})
  cf = coefficients(x)
  print(IOContext(io, :typeinfo => typeof(cf)), cf)
end

