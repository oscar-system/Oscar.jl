########################################################
# (1) Display
########################################################

# detailed printing
function Base.show(io::IO, ::MIME"text/plain", X::AffineVariety{<:Field,<:MPolyQuoRing})
    println(io, "Affine variety")
    print(io, " in ")
    println(io, ambient_space(X))
    println(io, "defined by")
    print(io, modulus(OO(X)))
end

# one line printing
function Base.show(io::IO, X::AffineVariety)
    if get(io, :supercompact, false)
      print(io, "Affine variety")
    else
      print(io, "Affine variety in ")
      print(io, ambient_space(X))
    end
end

# For affine space
function Base.show(io::IO, ::MIME"text/plain", X::AffineVariety{<:Field,<:MPolyRing})
    println(io, "Affine $(dim(X))-space")
    print(io, " over ")
    println(io, base_ring(X))
    println(io, "with coordinates")
    show(io, "text/plain", coordinates(X))
end

function Base.show(io::IO, X::AffineVariety{<:Field,<:MPolyRing})
    if get(io, :supercompact, false)
      print(io, "Affine space")
    else
      print(io, "Affine $(dim(X))-space over ")
      print(IOContext(io, :supercompact => true), base_ring(X))
    end
end
