########################################################
# (1) Display
########################################################


function Base.show(io::IO, ::MIME"text/plain",
                   X::AbsProjectiveAlgebraicSet{<:Field,<:MPolyQuoRing})
  io = pretty(io)
  println(io, "Algebraic set ")
  println(io, Indent(), "in ", Lowercase(), ambient_space(X))
  print(io, "defined by ", Lowercase(), defining_ideal(X))
  print(io, Dedent())
end


function Base.show(io::IO, X::AbsProjectiveAlgebraicSet{<:Field, <:MPolyQuoRing})
  io = pretty(io)
  if get(io, :supercompact, false)
    print(io, "Projective algebraic set")
  else
    print(io, "V(")
    if isdefined(X, :Xred)
      I = vanishing_ideal(X)
    else
      I = fat_ideal(X)
    end
    join(io, gens(I), ",")
    print(io,")")
  end
end

