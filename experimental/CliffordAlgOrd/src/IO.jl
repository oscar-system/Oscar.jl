##################################################
#
#  Clifford algebra and elements
#
##################################################
 
##### Algebra #####
function Base.show(io::IO, ::MIME"text/plain", C::CliffordAlgebra)
  io = pretty(io)
  print(io, "Clifford algebra of quadratic space with Gram matrix\n")
  print(io, Indent())
  show(io, "text/plain", gram_matrix(C))
  print(io, Dedent(), "\ndefined over $(base_ring(C))")
end
 
function Base.show(io::IO, C::CliffordAlgebra)
  if is_terse(io)
    print(io,"Clifford algebra")
  else
    print(terse(io), "Clifford algebra over ")
    print(terse(io), base_ring(C))
  end
end

##### Elements #####
function Base.show(io::IO, x::CliffordAlgebraElem)
  print(io, "[")
  foreach(y -> print(io, "$y "), coefficients(x)[1:(end - 1)])
  print(io, "$(coefficients(x)[end])]")
end

##################################################
#   
#  Clifford Order and elements
#
##################################################
  
##### Order #####
function Base.show(io::IO, ::MIME"text/plain", C::CliffordOrder)
  io = pretty(io)
  print(io, "Clifford order of even lattice over $(base_ring(C)) with Gram matrix\n")
  print(io, Indent())
  show(io, "text/plain", gram_matrix(C))
  print(io, Dedent(), "\nand coefficient ideals of the lattice\n")
  print(io, Indent())
  show(io, "text/plain", _coefficient_ideals_of_lattice(lattice(C)))
  print(io, Dedent())
end

function Base.show(io::IO, ::MIME"text/plain", C::ZZCliffordOrder)
  io = pretty(io)
  print(io, "Clifford order of even integer lattice with Gram matrix\n")
  print(io, Indent())
  show(io, "text/plain", gram_matrix(C))
  print(io, Dedent())
end

function Base.show(io::IO, C::Union{ZZCliffordOrder, CliffordOrder})
  if is_terse(io)
    print(io, "Clifford order")
  else
    print(terse(io), "Clifford order over ")
    print(terse(io), base_ring(C))
  end
end

##### Elements #####
function Base.show(io::IO, x::Union{ZZCliffordOrderElem, CliffordOrderElem})
  print(io, "[")
  foreach(y -> print(io, "$y "), coefficients(x)[1:(end - 1)])
  print(io, "$(coefficients(x)[end])]")
end

