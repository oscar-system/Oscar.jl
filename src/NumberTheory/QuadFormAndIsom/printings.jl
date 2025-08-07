###############################################################################
#
#  Lattices with isometry
#
###############################################################################

function Base.show(io::IO,  ::MIME"text/plain", Lf::ZZLatWithIsom)
  io = pretty(io)
  println(io, lattice(Lf))
  n = order_of_isometry(Lf)
  print(io, Indent())
  if is_finite(n)
    println(io, "with isometry of finite order $n")
  else
    println(io, "with isometry of infinite order")
  end
  println(io, "given by")
  show(io, MIME"text/plain"(), isometry(Lf))
  print(io, Dedent())
end

function Base.show(io::IO, Lf::ZZLatWithIsom)
  if is_terse(io)
    print(io, "Integer lattice with isometry")
  else
    n = order_of_isometry(Lf)
    if is_finite(n)
      print(io, "Integer lattice with isometry of finite order $n")
    else
      print(io, "Integer lattice with isometry of infinite order")
    end
  end
end

###############################################################################
#
#  Quadratic space with isometry
#
###############################################################################

function Base.show(io::IO, ::MIME"text/plain", Vf::QuadSpaceWithIsom)
  io = pretty(io)
  println(io, space(Vf))
  n = order_of_isometry(Vf)
  print(io, Indent())
  if is_finite(n)
    println(io, "with isometry of finite order $n")
  else
    println(io, "with isometry of infinite order")
  end
  println(io, "given by")
  show(io, MIME"text/plain"(), isometry(Vf))
  print(io, Dedent())
end

function Base.show(io::IO, Vf::QuadSpaceWithIsom)
  if is_terse(io)
    print(io, "Quadratic space with isometry")
  else
    n = order_of_isometry(Vf)
    if is_finite(n)
      print(io, "Quadratic space with isometry of finite order $n")
    else
      print(io, "Quadratic space with isometry of infinite order")
    end
  end
end
