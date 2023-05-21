function Base.show(io::IO,  ::MIME"text/plain", Lf::LatWithIsom)
  io = AbstractAlgebra.pretty(io)
  println(io, lattice(Lf))
  n = order_of_isometry(Lf)
  print(io, AbstractAlgebra.Indent())
  if is_finite(n)
    println(io, "with isometry of finite order $n")
  else
    println(io, "with isometry of infinite order")
  end
  println(io, "given by")
  print(IOContext(io, :compact => true), isometry(Lf))
  print(io, AbstractAlgebra.Dedent())
end

function Base.show(io::IO, Lf::LatWithIsom)
  if get(io, :supercompact, false)
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
