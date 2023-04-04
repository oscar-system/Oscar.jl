########################################################
# (1) Display
########################################################


function Base.show(io::IO, ::MIME"text/plain", X::AbsProjectiveAlgebraicSet{<:Field,<:MPolyQuoRing})
    println(io, "Vanishing locus ")
    println(io, "  in ", ambient_space(X))
    print(io, "  of ")
    Base.show(io, MIME("text/plain"), defining_ideal(X))
end


function Base.show(io::IO, X::AbsProjectiveAlgebraicSet)
    if get(io, :supercompact, false)
      print(io, "Projective algebraic set")
    else
      print(io, "Projective algebraic set in ")
      print(io, ambient_space(X))
    end
end

