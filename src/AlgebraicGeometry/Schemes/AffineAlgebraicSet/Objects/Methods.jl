########################################################
# (1) Display
########################################################

function Base.show(io::IO, ::MIME"text/plain", X::AffineAlgebraicSet{<:Field,<:MPolyQuoRing})
  io = pretty(io)
  println(io, "Affine algebraic set")  # at least one new line is needed
  println(io, Indent(), "in ", ambient_space(X))
  print(io, "defined by ", fat_ideal(X))
  # the last print statement must not add a new line
  print(io, Dedent()) # do not forget to Dedent() every Indent()
end

function Base.show(io::IO, X::AffineAlgebraicSet{<:Field,<:MPolyQuoRing})
  if get(io, :supercompact, false)
    # no nested printing
    print(io, "Affine algebraic set")
  else
    # nested printing allowed, preferably supercompact
    print(io, "V(")
    join(io, gens(fat_ideal(X)), ", ")
    print(io,")")
  end
end

# special case for Zariski opens
function Base.show(io::IO, ::MIME"text/plain", X::AffineAlgebraicSet)
  io = pretty(io)
  println(io, "Reduced subscheme")  # at least one new line is needed
  print(io, Indent(),"of ", Lowercase(), fat_scheme(X))
  print(io, Dedent())
  # the last print statement must not add a new line
end

function Base.show(io::IO, X::AffineAlgebraicSet)
  io = pretty(io)
  println(io, "Reduced subscheme of ", Lowercase(), fat_scheme(X))  # at least one new line is needed
  # the last print statement must not add a new line
end

