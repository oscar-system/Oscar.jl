########################################################
# (1) Display
########################################################

# detailed printing
function Base.show(io::IO, ::MIME"text/plain", X::ProjectiveVariety{<:Field,<:MPolyQuoRing})
    println(io, "Projective variety")
    print(io, "  in ")
    println(io, ambient_space(X))
    println(io, "defined by")
    print(io, modulus(OO(X)))
end

# one line printing
function Base.show(io::IO, X::ProjectiveVariety)
    if get(io, :supercompact, false)
      print(io, "Projective variety")
    else
      print(io, "Projective variety in ")
      print(io, ambient_space(X))
    end
end

# For affine space
function Base.show(io::IO, ::MIME"text/plain", X::ProjectiveVariety{<:Field,<:MPolyRing})
    println(io, "Projective $(dim(X))-space")
    print(io, " over ")
    println(io, base_ring(X))
    println(io, "with coordinates")
    show(io, "text/plain", coordinates(X))
end

function Base.show(io::IO, X::ProjectiveVariety{<:Field,<:MPolyRing})
    if get(io, :supercompact, false)
      print(io, "Projective space")
    else
      print(io, "Projective $(dim(X))-space over ")
      print(IOContext(io, :supercompact => true), base_ring(X))
    end
end
