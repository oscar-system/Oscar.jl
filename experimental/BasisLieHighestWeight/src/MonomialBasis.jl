@attributes mutable struct MonomialBasis{T<:ModuleData}
  V::T
  birational_seq::BirationalSequence
  monomial_ordering::MonomialOrdering
  monomials::Set{ZZMPolyRingElem}
  monomials_parent::ZZMPolyRing

  function MonomialBasis(
    V::T,
    birational_seq::BirationalSequence,
    monomial_ordering::MonomialOrdering,
    monomials::Set{ZZMPolyRingElem},
  ) where {T<:ModuleData}
    @req dim(V) == length(monomials) "dimesion mismatch"
    return new{T}(
      V,
      birational_seq,
      monomial_ordering,
      monomials,
      parent(first(monomials)),
    )
  end
end

base_lie_algebra(basis::MonomialBasis) = base_lie_algebra(basis.V)

highest_weight(basis::MonomialBasis{<:SimpleModuleData}) = highest_weight(basis.V)

dim(basis::MonomialBasis) = dim(basis.V)
length(basis::MonomialBasis) = dim(basis)

monomials(basis::MonomialBasis) = basis.monomials

monomial_ordering(basis::MonomialBasis) = basis.monomial_ordering

birational_sequence(basis::MonomialBasis) = basis.birational_seq

function Base.show(io::IO, ::MIME"text/plain", basis::MonomialBasis{<:SimpleModuleData})
  io = pretty(io)
  print(io, "Monomial basis of a highest weight module")
  print(
    io,
    Indent(),
    "\nof highest weight $(highest_weight(basis))",
    Dedent(),
  )
  print(io, Indent(), "\nof dimension $(dim(basis))", Dedent())
  print(io, Indent(), "\nwith monomial ordering $(monomial_ordering(basis))", Dedent())
  print(io, "\nover ", Lowercase(), base_lie_algebra(basis))
  if get_attribute(basis, :algorithm, nothing) in
    (basis_lie_highest_weight_compute, basis_coordinate_ring_kodaira_compute)
    println(
      io,
      "\nwhere the used birational sequence consists of the following roots:",
    )
    op_roots = operators_as_roots(birational_sequence(basis))
    print(IOContext(io, :typeinfo => typeof(op_roots)), Indent(), op_roots)
    print(io, Dedent())
  end
end

function Base.show(io::IO, basis::MonomialBasis{<:SimpleModuleData})
  if is_terse(io)
    print(io, "Monomial basis of a highest weight module")
  else
    io = pretty(io)
    print(
      io,
      "Monomial basis of a highest weight module with highest weight $(highest_weight(basis)) over ",
    )
    print(terse(io), Lowercase(), base_lie_algebra(basis))
  end
end

function Base.show(io::IO, ::MIME"text/plain", basis::MonomialBasis{<:DemazureModuleData})
  io = pretty(io)
  print(io, "Monomial basis of a Demazure module")
  print(
    io,
    Indent(),
    "\nof extremal weight ($(highest_weight(basis.V))) * $(weyl_group_elem(basis.V))",
    Dedent(),
  )
  print(io, Indent(), "\nof dimension $(dim(basis))", Dedent())
  print(io, Indent(), "\nwith monomial ordering $(monomial_ordering(basis))", Dedent())
  print(io, "\nover ", Lowercase(), base_lie_algebra(basis))
  if get_attribute(basis, :algorithm, nothing) in
    (basis_lie_highest_weight_compute, basis_coordinate_ring_kodaira_compute)
    println(
      io,
      "\nwhere the used birational sequence consists of the following roots:",
    )
    op_roots = operators_as_roots(birational_sequence(basis))
    print(IOContext(io, :typeinfo => typeof(op_roots)), Indent(), op_roots)
    print(io, Dedent())
  end
end

function Base.show(io::IO, basis::MonomialBasis{<:DemazureModuleData})
  if is_terse(io)
    print(io, "Monomial basis of a Demazure module")
  else
    io = pretty(io)
    print(
      io,
      "Monomial basis of a Demazure module with extremal weight ($(highest_weight(basis.V))) * $(weyl_group_elem(basis.V)) over ",
    )
    print(terse(io), Lowercase(), base_lie_algebra(basis))
  end
end
