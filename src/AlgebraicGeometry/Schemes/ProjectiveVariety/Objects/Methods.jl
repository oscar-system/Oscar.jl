########################################################
# (1) Display
########################################################

# detailed printing
function Base.show(io::IO, ::MIME"text/plain", X::AbsProjectiveVariety{<:Field,<:MPolyQuoRing})
    println(io, "Projective variety")
    print(io, "  in ")
    println(io, ambient_space(X))
    print(io, "  defined by ")
    print(io, defining_ideal(X))
end

# one line printing
function Base.show(io::IO, X::AbsProjectiveVariety{<:Field,<:MPolyQuoRing})
    if get(io, :supercompact, false)
      print(io, "Projective variety")
    else
      print(io, "Projective variety defined by ")
      print(io, vanishing_ideal(X))
    end
end


# Projective space
Base.show(io::IO, ::MIME"text/plain", P::AbsProjectiveVariety{<:Field, <:MPolyDecRing}) = Base.show(io, MIME("text/plain"),underlying_scheme(P))
Base.show(io::IO, P::AbsProjectiveVariety{<:Field, <:MPolyDecRing}) = Base.show(io, underlying_scheme(P))
