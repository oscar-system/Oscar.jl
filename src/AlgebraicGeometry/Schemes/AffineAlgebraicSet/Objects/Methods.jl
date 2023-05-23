########################################################
# (1) Display
########################################################

function Base.show(io::IO, ::MIME"text/plain", X::AffineAlgebraicSet{<:Field,<:MPolyQuoRing})
  println(io, "Vanishing locus")  # at least one new line is needed
  println(io, "  in $(ambient_space(X))")
  print(io, "  of $(vanishing_ideal(X))")
  # the last print statement must not add a new line
end

function Base.show(io::IO, X::AffineAlgebraicSet{<:Field,<:MPolyQuoRing})
  if get(io, :supercompact, false)
    # no nested printing
    print(io, "Affine algebraic set")
  else
    # nested printing allowed, preferably supercompact
    print(io, "Vanishing locus of $(vanishing_ideal(X))")
  end
end

# special case for Zariski opens
function Base.show(io::IO, ::MIME"text/plain", X::AffineAlgebraicSet)
  println(io, "Affine algebraic set")  # at least one new line is needed
  print(io, "  with coordinate ring $(OO(X))")
  # the last print statement must not add a new line
end


