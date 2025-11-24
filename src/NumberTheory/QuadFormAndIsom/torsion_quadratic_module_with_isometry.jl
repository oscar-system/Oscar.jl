###############################################################################
#
#  Accessors
#
###############################################################################

underlying_module(Tf::TorQuadModuleWithIsom) = Tf.T

isometry(Tf::TorQuadModuleWithIsom) = Tf.f

###############################################################################
#
#  Attributes
#
###############################################################################

@attr Int function order_of_isometry(Tf::TorQuadModuleWithIsom)
  T = underlying_module(Tf)
  f = isometry(Tf)
  if isone(order(T))
    return Int(1)
  end
  n = 1
  _f = f
  while !isone(matrix(_f))
    n += 1
    _f = compose(_f, f)
  end
  return n
end

###############################################################################
#
#  Constructors
#
###############################################################################

function torsion_quadratic_module_with_isometry(
  T::TorQuadModule,
  f::TorQuadModuleMap;
  check::Bool=true,
)
  if check
    @req domain(f) === codomain(f) === T "Wrong domain or codomain"
    @req Hecke.is_isometry(f) "Not an isometry"
  end

  if isone(order(T))
    return TorQuadModuleWithIsom(T, id_hom(T))
  end

  return TorQuadModuleWithIsom(T, f)
end

function torsion_quadratic_module_with_isometry(
  T::TorQuadModule,
  f::AutomorphismGroupElem{TorQuadModule};
)
  @req domain(f) === T "Wrong domain"
  return torsion_quadratic_module_with_isometry(T, hom(f); check=false)
end

function torsion_quadratic_module_with_isometry(
  T::TorQuadModule,
  m::ZZMatrix;
  check::Bool=true,
)
  f = hom(T, T, m; check)
  return torsion_quadratic_module_with_isometry(T, f; check)
end

function torsion_quadratic_module_with_isometry(
  T::TorQuadModule,
  g::FinGenAbGroupHom;
  check::Bool=true,
)
  if check
    @req domain(g) === codomain(g) === abelian_group(T) "Wrong domain or codomain"
  end
  f = TorQuadModuleMap(T, T, g)
  return torsion_quadratic_module_with_isometry(T, f; check)
end

function torsion_quadratic_module_with_isometry(
  T::TorQuadModule,
  g::MatrixGroupElem{QQFieldElem, QQMatrix};
  check::Bool=true,
)
  if check
    B = basis_matrix(relations(T))
    @req can_solve(B, B*matrix(g); side=:left) "Incompatible inputs"
  end
  f = hom(T, T, elem_type(T)[T(lift(t)*matrix(g)) for t in gens(T)])
  return torsion_quadratic_module_with_isometry(T, f; check)
end

function torsion_quadratic_module_with_isometry(
  q::QQMatrix,
  f::Union{TorQuadModuleMap, ZZMatrix, FinGenAbGroupHom, MatrixGroupElem{QQFieldElem, QQMatrix}};
  check::Bool=true,
)
  T = TorQuadModule(q)
  return torsion_quadratic_module_with_isometry(T, f; check)
end

###############################################################################
#
#  Submodules
#
###############################################################################

function _sub(Tf::TorQuadModuleWithIsom, i::TorQuadModuleMap)
  T = underlying_module(Tf)
  f = isometry(Tf)
  @req codomain(i) === T "Wrong codomain"
  S = domain(i)
  L = cover(S)
  imgs = elem_type(S)[]
  for a in gens(S)
    b = lift(f(i(a)))
    @req b in L "Submodule not preserved by the isometry"
    push!(imgs, S(b))
  end
  g = hom(S, S, gens(S), imgs)
  return torsion_quadratic_module_with_isometry(S, g; check=false)
end

function _is_submodule(Tf::TorQuadModuleWithIsom, S::TorQuadModule)
  is_sublattice(cover(underlying_module(Tf)), cover(S)) || return false
  f = isometry(Tf)
  for a in gens(S)
    lift(a) in cover(S) || return false
  end
  return true
end

function sub(
  Tf::TorQuadModuleWithIsom,
  gene::Vector{TorQuadModuleElem}
)
  T = underlying_module(Tf)
  S, i = sub(T, gene)
  Sg = _sub(Tf, i)
  return Sg, i
end

function primary_part(Tf::TorQuadModuleWithIsom, m::ZZRingElem)
  S, i = primary_part(underlying_module(Tf), m)
  Sg = _sub(Tf, i)
  return Sg, i
end

primary_part(Tf::TorQuadModuleWithIsom, m::Int) = primary_part(Tf, ZZ(m))

function orthogonal_submodule(
  Tf::TorQuadModuleWithIsom,
  S::TorQuadModule;
  check::Bool=true,
)
  if check
    @req _is_submodule(Tf, S) "Not a stable submodule"
  end
  K, i = orthogonal_submodule(underlying_module(Tf), S; check=false)
  Kg = _sub(Tf, i)
  return Kg, i
end

function submodules(
  Tf::TorQuadModuleWithIsom;
  kwargs...,
)
  T = underlying_module(Tf)
  f = isometry(Tf)
  _res = stable_submodules(T, [f]; kwargs...)
  return ((_sub(Tf, j), j) for (_, j) in _res)
end

###############################################################################
#
#  Isometry
#
###############################################################################

function automorphism_group(Tf::TorQuadModuleWithIsom)
  T = underlying_module(Tf)
  f = isometry(Tf)
  O, _ = centralizer_in_orthogonal_group(T, f; check=false)
  return O
end

function is_isomorphic_with_map(
  Tf::TorQuadModuleWithIsom,
  Sg::TorQuadModuleWithIsom,
)
  T = underlying_module(Tf)
  S = underlying_module(Sg)
  
  f = isometry(Tf)
  ok, phi = is_isometric_with_isometry(T, S)
  !ok && return false, phi

  fphi = inv(phi) * f * phi
  OS = orthogonal_group(S)
  gS = OS(isometry(Sg); check=false)
  fS = OS(fphi; check=false)
  ok, h = is_conjugate_with_data(OS, gS, fS)
  !ok && return false, hom(T, S, zero_matrix(ZZ, ngens(T), ngens(S)))

  phi = phi * hom(h)

  @hassert :ZZLatWithIsom 2 inv(phi) * f * phi == isometry(Sg)
  return true, phi
end

###############################################################################
#
#  Equality and hash
#
###############################################################################

function Base.:(==)(T1::TorQuadModuleWithIsom, T2::TorQuadModuleWithIsom)
  underlying_module(T1) == underlying_module(T2) || return false
  return isometry(T1) == isometry(T2)
end

function Base.hash(T::TorQuadModuleWithIsom, u::UInt)
  u = Base.hash(underlying_module(T), u)
  return Base.hash(isometry(T), u)
end
