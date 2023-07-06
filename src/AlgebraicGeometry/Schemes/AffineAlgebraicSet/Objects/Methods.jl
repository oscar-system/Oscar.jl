########################################################
# (1) Display
########################################################

function Base.show(io::IO, ::MIME"text/plain", X::AffineAlgebraicSet{<:Field,<:MPolyQuoRing})
  io = pretty(io)
  println(io, "Affine algebraic set")
  println(io, Indent(), "in ", Lowercase(), ambient_space(X))
  print(io, Dedent(), "defined by ", fat_ideal(X))
end

function Base.show(io::IO, X::AffineAlgebraicSet{<:Field,<:MPolyQuoRing})
  if get(io, :supercompact, false)
    print(io, "Scheme")
  elseif get_attribute(X, :is_empty, false)
    print(io, "Empty affine algebraic set")
  else
    print(io, "V(")
    join(io, gens(fat_ideal(X)), ", ")
    print(io,")")
  end
end

# special case for Zariski opens
function Base.show(io::IO, ::MIME"text/plain", X::AffineAlgebraicSet)
  io = pretty(io)
  println(io, "Reduced subscheme")
  print(io, Indent(),"of ", Lowercase(), fat_scheme(X))
  print(io, Dedent())
end

function Base.show(io::IO, X::AffineAlgebraicSet)
  io = pretty(io)
  if get(io, :supercompact, false)
    print(io, "Scheme")
  elseif get_attribute(X, :is_empty, false)
    print(io, "Empty affine algebraic set")
  else
    print(io, "Reduced subscheme of ", Lowercase(), fat_scheme(X))
  end
end

