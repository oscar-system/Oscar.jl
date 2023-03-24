########################################################
# (1) Display
########################################################

function Base.show(io::IO, X::AffineVariety)
  d = dim(ambient_space(X))
  println(io, "Affine variety in $ambient_space(X)")
  println(io, "defined by")
  println(io, modulus(OO(X)))
end


