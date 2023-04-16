########################################################
# (1) Display
########################################################

# detailed printing
function Base.show(io::IO, ::MIME"text/plain", X::AffineVariety{<:Field,<:MPolyQuoRing})
  println(io, "Affine variety")
  print(io, " in ")
  println(io, ambient_space(X))
  print(io, "defined by ")
  print(io, vanishing_ideal(X))
end

function Base.show(io::IO, X::AffineVariety{<:Field,<:MPolyQuoRing})
# one line printing
  if get(io, :supercompact, false)
    print(io, "Affine variety")
  else
    print(io, "Affine variety defined by $(vanishing_ideal(X))")
  end
end


# The case that X is not Zariski closed in its ambient space.
function Base.show(io::IO, ::MIME"text/plain", X::AffineVariety)
  println(io, "Affine variety")
  print(io, "  in ")
  println(io, ambient_space(X))
  print(io, "with coordinate ring $(OO(X))")
end

function Base.show(io::IO, X::AffineVariety)
  if get(io, :supercompact, false)
    print(io, "Affine variety")
  else
    print(io, "Affine variety contained in the vanishing locus of $(vanishing_ideal(X))")
  end
end


# For affine space
function Base.show(io::IO, ::MIME"text/plain", X::AffineVariety{<:Field,<:MPolyRing})
  println(io, "Affine $(dim(X))-space")
  print(io, "  over ")
  println(io, base_ring(X))
  println(io, "with coordinates")
  show(io, "text/plain", coordinates(X))
end

function Base.show(io::IO, X::AffineVariety{<:Field,<:MPolyRing})
  if get(io, :supercompact, false)
    print(io, "AA^$(dim(X))")
  else
    print(io, "Affine $(dim(X))-space over ")
    print(IOContext(io, :supercompact => true), base_ring(X))
  end
end
