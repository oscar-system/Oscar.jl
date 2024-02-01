function Base.show(io::IO, ::MIME"text/plain", phi::AbsHyperComplexMorphism)
  io = pretty(io)
  println(io, "Morphism of complexes from")
  print(io, Indent())
  print(io, "$(domain(phi))")
  println(io, Dedent())
  println(io, "to")
  print(io, Indent())
  print(io, "$(codomain(phi))")
  print(io, Dedent())
end


