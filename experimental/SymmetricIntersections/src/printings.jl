###############################################################################
#
#  Elevators
#
###############################################################################

function Base.show(io::IO, EC::ElevCtx{T, U}) where {T, U}
  d = degree_of_elevations(EC)
  if is_terse(io)
    print(io, "Elevator")
  else
    print(io, "Degree $d elevator of a list with entries of type $T")
  end
end

###############################################################################
#
#  Representations
#
###############################################################################

# Linear representations

function Base.show(io::IO, ::MIME"text/plain", LR::LinRep)
  println(io, "Linear representation")
  io = pretty(io)
  println(io, Indent(), "of ", underlying_group(LR))
  println(io, "over ", Lowercase(), base_field(representation_ring(LR)))
  print(io, Dedent(), "of dimension ", dimension_representation(LR))
end 

function Base.show(io::IO, LR::LinRep)
  if is_terse(io)
    print(io, "Linear representation")
  else
    print(io, "Linear representation of finite group of dimension $(dimension_representation(LR))")
  end
end

# Projective representations

function Base.show(io::IO, ::MIME"text/plain", PR::ProjRep)
  println(io, "Projective representation")
  io = pretty(io)
  println(io, Indent(), "of ", underlying_group(PR))
  println(io, "over ", Lowercase(), base_field(representation_ring_linear_lift(PR)))
  print(io, Dedent(), "of dimension ", dimension_representation(PR))
end

function Base.show(io::IO, PR::ProjRep)
  if is_terse(io)
    print(io, "Projective representation")
  else
    print(io, "Projective representation of finite group of dimension ", dimension_representation(PR))
  end
end

# Representation rings

function Base.show(io::IO, ::MIME"text/plain", RR::RepRing)
  println(io, "Representation ring")
  io = pretty(io)
  println(io, Indent(), "of ", underlying_group(RR))
  print(io, "over ", Lowercase(), base_field(RR))
  print(io, Dedent())
end

function Base.show(io::IO, RR::RepRing)
  if is_terse(io)
    print(io, "Representation ring")
  else
    print(io, "Representation ring of finite group over a field of characteristic 0")
  end
end

###############################################################################
#
#  Symmetric Grassmannians
#
###############################################################################

# Isotypical Grassmannians

function Base.show(io::IO, ::MIME"text/plain", M::IsotGrass)
  chi = submodule_character(M)
  println(io, "Symmetric Grassmannian of $(degree(chi))-dimensional submodules")
  io = pretty(io)
  println(io, Indent(), "of ", Lowercase(), module_representation(M))
  println(io, Dedent(), "with isotypical character")
  print(io, chi)
end

function Base.show(io::IO, M::IsotGrass)
  if is_terse(io)
    print(io, "Isotypical Grassmannian")
  else
    print(io, "Isotypical Grassmannian of dimension $(projective_dimension(M))")
  end
end

# Character Grassmannians

function Base.show(io::IO, ::MIME"text/plain", M::CharGrass)
  chi = submodule_character(M)
  println(io, "Symmetric Grassmannian of $(degree(chi))-dimensional submodules")
  io = pretty(io)
  print(io, Indent(), "of ", Lowercase(), module_representation(M))
  println(io)
  println(io, Dedent(), "with character")
  print(io, chi)
end

function Base.show(io::IO, M::CharGrass)
  if is_terse(io)
    print(io, "Character Grassmannian")
  else
    print(io, "Character Grassmannian of dimension $(projective_dimension(M))")
  end
end

# Determinant Grassmannians

function Base.show(io::IO, ::MIME"text/plain", M::DetGrass)
  chi = submodule_determinant_character(M)
  println(io, "Symmetric Grassmannian of $(M.d)-dimensional submodules")
  io = pretty(io)
  println(io, Indent(), "of ", Lowercase(), module_representation(M))
  println(io, Dedent(), "with determinant character")
  print(io, chi)
end

function Base.show(io::IO, M::DetGrass)
  if is_terse(io)
    print(io, "Determinant Grassmannian")
  else
    print(io, "Determinant Grassmannian of dimension $(projective_dimension(M))")
  end
end

# Invariant Grassmannians

function Base.show(io::IO, ::MIME"text/plain", M::InvGrass)
  println(io, "Symmetric Grassmannian of $(M.d)-dimensional submodules")
  io = pretty(io)
  print(io, Indent(), "of ", Lowercase(), module_representation(M))
  print(io, Dedent())
end

function Base.show(io::IO, M::InvGrass)
  if is_terse(io)
    print(io, "Invariant Grassmannian")
  else
    print(io, "Invariant Grassmannian of dimension $(projective_dimension(M))")
  end
end

###############################################################################
#
#  Symmetric intersections
#
###############################################################################

function Base.show(io::IO, ::MIME"text/plain", symci::SymInter)
  M = symci.para
  j = symci.j
  RR = representation_ring_linear_lift(symci.prep)
  F = base_field(RR)
  G = underlying_group(symci.prep)
  n = ngens(codomain(j))
  t = submodule_dimension(M)
  d = total_degree(j(gens(domain(j))[1]))
  if t == 1
    ty = "($d)"
  elseif t == 2
    ty = "($d, $d)"
  else
    ty = "($d, "
    for i in 2:t-1
      ty *= "$d, "
    end
    ty *= "$d)"
  end
  io = pretty(io)
  println(io, "Parameter space for intersections")
  println(io, Indent(), "of type $(ty)")
  println(io, "in projective $(n-1)-space")
  print(io, Indent(), "over ", Lowercase())
  Base.show(terse(io), F)
  println(io)
  print(io, Dedent(), Dedent(), "preserved under the action of ", Lowercase(), G)
end

function Base.show(io::IO, symci::SymInter)
  if is_terse(io)
    print(io, "Symmetric intersections")
  else
    print(io, "Parameter space for symmetric intersections")
  end
end

