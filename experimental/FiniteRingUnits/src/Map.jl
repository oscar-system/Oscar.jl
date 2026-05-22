domain(f::RingMultMap) = f.A

codomain(f::RingMultMap) = f.R

(f::RingMultMap)(x) = image(f, x)

function image(f::RingMultMap, x)
  @req parent(x) === domain(f) "Element must be contained in domain"
  return f.g(x)::elem_type(codomain(f))
end

function preimage(f::RingMultMap, y)
  @req parent(y) === codomain(f) "Element must be contained in codomain"
  @req Oscar.is_invertible(y)[1] "Element must be a unit"
  x = f.f(y)
  @assert parent(x) === domain(f)
  return x::elem_type(domain(f))
end

function Base.show(io::IO, M::RingMultMap)
   if Oscar.is_terse(io)
      # no nested printing
      print(io, "Map from multiplicative group of ring")
   else
      io = Oscar.pretty(io)
      io = Oscar.terse(io)
      print(io, "Map: ", Oscar.AbstractAlgebra.Lowercase(), domain(M))
      print(io, " -> ", Oscar.AbstractAlgebra.Lowercase(), codomain(M))
   end
end

function Oscar.AbstractAlgebra.show_map_data(io::IO, M::RingMultMap)
  #println(io)
  #println(io, "with", Oscar.AbstractAlgebra.Indent())
  #D = codomain(M)
  #for d in gens(D)
  #  println(io, preimage(M, d), " => ", d)
  #end
end
