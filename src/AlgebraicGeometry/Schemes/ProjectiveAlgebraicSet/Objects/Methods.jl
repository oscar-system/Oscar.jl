########################################################
# (1) Display
########################################################


function Base.show(io::IO, ::MIME"text/plain",
                   X::AbsProjectiveAlgebraicSet{<:Field,<:MPolyQuoRing})
  println(io, "Vanishing locus ")
  println(io, "  in ", ambient_space(X))
  print(io, "  of $(defining_ideal(X))")
end


function Base.show(io::IO, X::AbsProjectiveAlgebraicSet{<:Field, <:MPolyQuoRing})
  if get(io, :supercompact, false)
    print(io, "Projective algebraic set")
  else
    print(io, "Vanishing locus in ")
    print(IOContext(io, :supercompact => true), ambient_space(X))
    print(io, " of ", vanishing_ideal(X))
  end
end

