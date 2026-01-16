###############################################################################
#
#  Accessors
#
###############################################################################

@doc raw"""
    underlying_module(Tf::TorQuadModuleWithIsom) -> TorQuadModule

Given a torsion quadratic module with isometry ``(T, f)``, return ``T``.

# Examples

```jldoctest
julia> Tf = torsion_quadratic_module_with_isometry(QQ[-1//60;], ZZ[11;])
Finite quadratic module of order 60
  with 1 generator
  with isometry given by
  [11]

julia> underlying_module(Tf)
Finite quadratic module
  over integer ring
Abelian group: Z/60
Bilinear value module: Q/Z
Quadratic value module: Q/2Z
Gram matrix quadratic form:
[119//60]
```
"""
underlying_module(Tf::TorQuadModuleWithIsom) = Tf.T

@doc raw"""
    torsion_quadratic_module(Tf::TorQuadModuleWithIsom) -> TorQuadModule

Alias for [`underlying_module(::TorQuadModuleWithIsom)`](@ref).
"""
torsion_quadratic_module(Tf::TorQuadModuleWithIsom) = underlying_module(Tf)

@doc raw"""
    isometry(Tf::TorQuadModuleWithIsom) -> TorQuadModuleMap

Given a torsion quadratic module with isometry ``(T, f)``, return ``f``.

# Examples

```jldoctest
julia> Tf = torsion_quadratic_module_with_isometry(QQ[-1//60;], ZZ[11;])
Finite quadratic module of order 60
  with 1 generator
  with isometry given by
  [11]

julia> isometry(Tf)
Map
  from finite quadratic module: Z/60 -> Q/2Z
  to finite quadratic module: Z/60 -> Q/2Z
```
"""
isometry(Tf::TorQuadModuleWithIsom) = Tf.f

###############################################################################
#
#  Attributes
#
###############################################################################

@doc raw"""
    order_of_isometry(Tf::TorQuadModuleWithIsom) -> Int

Given a torsion quadratic module with isometry ``(T, f)``, return the order
of the isometry ``f``.

# Examples

```jldoctest
julia> Tf = torsion_quadratic_module_with_isometry(QQ[-1//60;], ZZ[11;])
Finite quadratic module of order 60
  with 1 generator
  with isometry given by
  [11]

julia> order_of_isometry(Tf)
2
```
"""
@attr ZZRingElem function order_of_isometry(Tf::TorQuadModuleWithIsom)
  T = underlying_module(Tf)
  f = isometry(Tf)
  C = Oscar._orthogonal_group(T, [f]; check=false)
  return order(C)
end

###############################################################################
#
#  Constructors
#
###############################################################################

@doc raw"""
    torsion_quadratic_module_with_isometry(
      T::TorQuadModule,
      [f::U];
      check::Bool=true,
    ) -> TorQuadModuleWithIsom

Given a torsion quadratic module ``T`` and given ``f`` being either:
  * an isometry of ``T``
    (`U <: AutomorphismGroupElem{TorQuadModule}`),
  * an abelian group endomorphism of ``T``
    (`U <: Union{TorQuadModuleMap, FinGenAbGroupHom, ZZMatrix}`),
  * or an isometry of the cover ``L`` of ``T``
    (`U <: MatGroupElem{QQFieldElem, QQMatrix}`)
return the associated torsion quadratic module with isometry ``(T, f')``,
where ``f' <: TorQuadModuleMap`` is the isometry of ``T`` defined by ``f``.

If ``f`` is not supplied, it is set to be the identity map of ``T`` by
default.

If `check` is set to `true`, the function checks whether ``f`` indeed defines
an isometry of ``T``.
"""
torsion_quadratic_module_with_isometry(::TorQuadModule, ::TorQuadModuleMap)

function torsion_quadratic_module_with_isometry(
  T::TorQuadModule,
  f::TorQuadModuleMap;
  check::Bool=true,
)
  @req is_finite(abelian_group(T)) "Only available for finite modules"
  if check
    @req domain(f) === codomain(f) === T "Wrong domain or codomain"
    @req Hecke.is_isometry(f) "Not an isometry"
  end

  return TorQuadModuleWithIsom(T, f)
end

function torsion_quadratic_module_with_isometry(T::TorQuadModule)
  @req is_finite(abelian_group(T)) "Only available for finite modules"
  f = id_hom(T)
  return TorQuadModuleWithIsom(T, f)
end

function torsion_quadratic_module_with_isometry(
  T::TorQuadModule,
  f::AutomorphismGroupElem{TorQuadModule};
)
  @req is_finite(abelian_group(T)) "Only available for finite modules"
  @req domain(f) === T "Wrong domain"
  return torsion_quadratic_module_with_isometry(T, hom(f); check=false)
end

function torsion_quadratic_module_with_isometry(
  T::TorQuadModule,
  m::ZZMatrix;
  check::Bool=true,
)
  @req is_finite(abelian_group(T)) "Only available for finite modules"
  f = hom(T, T, m; check)
  return torsion_quadratic_module_with_isometry(T, f; check)
end

function torsion_quadratic_module_with_isometry(
  T::TorQuadModule,
  g::FinGenAbGroupHom;
  check::Bool=true,
)
  @req is_finite(abelian_group(T)) "Only available for finite modules"
  @req domain(g) === codomain(g) === abelian_group(T) "Wrong domain or codomain"
  f = TorQuadModuleMap(T, T, g)
  return torsion_quadratic_module_with_isometry(T, f; check)
end

function torsion_quadratic_module_with_isometry(
  T::TorQuadModule,
  g::MatGroupElem{QQFieldElem, QQMatrix};
  check::Bool=true,
)
  @req is_finite(abelian_group(T)) "Only available for finite modules"
  if check
    B = basis_matrix(relations(T))
    @req can_solve(B, B*matrix(g); side=:left) "Incompatible inputs"
  end
  f = hom(T, T, elem_type(T)[T(lift(t)*matrix(g)) for t in gens(T)])
  return torsion_quadratic_module_with_isometry(T, f; check)
end

@doc raw"""
    torsion_quadratic_module_with_isometry(
      q::QQMatrix,
      [f::ZZMatrix];
      check::Bool=true,
    ) -> TorQuadModuleWithIsom

Given a symmetric matrix ``q`` with rational entries and a matrix ``f``
with integer entries, return the pair ``(T, f')`` where:
  - ``T`` is a torsion quadratic module with Gram matrix ``T``, for the
    associated finite quadratic form;
  - ``f'`` is the isometry of ``T`` defined by ``f``.

If ``f`` is not supplied, then it is set to be the identity matrix by
default.

If `check` is set to `true`, the function tests whether ``f`` defines
an isometry of ``T``.
"""
function torsion_quadratic_module_with_isometry(
  q::QQMatrix,
  f::ZZMatrix;
  check::Bool=true,
)
  T = TorQuadModule(q)
  return torsion_quadratic_module_with_isometry(T, f; check)
end

function torsion_quadratic_module_with_isometry(q::QQMatrix)
  f = identity_matrix(ZZ, ncols(q))
  T = TorQuadModule(q)
  return torsion_quadratic_module_with_isometry(T, f; check=false)
end

###############################################################################
#
#  Submodules
#
###############################################################################

@doc raw"""
    _sub(Tf::TorQuadModuleWithIsom, i::TorQuadModuleMap) -> TorQuadModuleWithIsom

Given a torsion quadratic module with isometry ``(T, f)`` and an abelian group
homomorphism ``i``, of type ``TorQuadModuleMap``, whose image is preserved by
``f``, return the pair ``(S, g)`` where ``S`` is the image of ``i`` and ``g``
is the restriction of ``f`` to ``S``.
"""
function _sub(Tf::TorQuadModuleWithIsom, i::TorQuadModuleMap)
  T = underlying_module(Tf)
  f = isometry(Tf)
  @req codomain(i) === T "Wrong codomain"
  S = domain(i)
  L = cover(S)
  imgs = elem_type(S)[]
  for a in gens(S)
    ok, b = has_preimage_with_preimage(i, f(i(a)))
    @req ok "Submodule not preserved by the isometry"
    push!(imgs, b)
  end
  g = hom(S, S, gens(S), imgs)
  return torsion_quadratic_module_with_isometry(S, g; check=false)
end

@doc raw"""
    _is_submodule(Tf::TorQuadModuleWithIsom, S::TorQuadModule) -> Bool

Given a torsion quadratic module with isometry ``(T, f)``, return whether
``S`` is an ``f``-stable submodule of ``T``.
"""
function _is_submodule(Tf::TorQuadModuleWithIsom, S::TorQuadModule)
  T = underlying_module(Tf)
  is_sublattice(cover(T), cover(S)) || return false
  f = isometry(Tf)
  for a in gens(S)
    lift(f(T(lift(a)))) in cover(S) || return false
  end
  return true
end

@doc raw"""
    sub(
      Tf::TorQuadModuleWithIsom,
      gene::Vector{TorQuadModuleElem},
    ) -> TorQuadModuleWithIsom, TorQuadModuleMap

Given a torsion quadratic module with isometry ``(T, f)`` and a list `gene`
of elements of ``T``, return the pair ``(S, g)`` where ``S`` is the submodule
of ``T`` generated by `gene` and ``g`` is the restriction of ``f`` to ``S``.
The second output is the embedding of ``S`` in ``T``.

An error is thrown is ``S`` is not preserved by ``f``.
"""
function sub(
  Tf::TorQuadModuleWithIsom,
  gene::Vector{TorQuadModuleElem}
)
  T = underlying_module(Tf)
  S, i = sub(T, gene)
  Sg = _sub(Tf, i)
  return Sg, i
end

@doc raw"""
    primary_part(
      Tf::TorQuadModuleWithIsom,
      m::IntegerUnion,
    ) -> TorQuadModuleWithIsom, TorQuadModuleMap

Given a torsion quadratic module with isometry ``(T, f)``, return the
``m``-primary part ``S`` of ``T`` together with the restriction of ``f`` to it.
The second output is the embedding of ``S`` in ``T``.
"""
primary_part(::TorQuadModuleWithIsom, ::IntegerUnion)

function primary_part(Tf::TorQuadModuleWithIsom, m::ZZRingElem)
  S, i = primary_part(underlying_module(Tf), m)
  Sg = _sub(Tf, i)
  return Sg, i
end

primary_part(Tf::TorQuadModuleWithIsom, m::Int) = primary_part(Tf, ZZ(m))

@doc raw"""
    orthogonal_submodule(
      Tf::TorQuadModuleWithIsom,
      S::TorQuadModule;
      check::Bool=true,
    ) -> TorQuadModuleWithIsom, TorQuadModuleMap

Given a torsion quadratic module with isometry ``(T, f)`` and a submodule
``S`` of ``T`` which is preserved by ``f``, return the pair ``(K, g)`` where
``K`` is the orthogonal complement of ``S`` in ``T`` and ``g`` is the
restriction of ``f`` to ``K``. The second output is the embedding of ``K`` in
``T``.

If `check` is set to `true`, the function checks whether ``S`` is indeed an
``f``-stable submodule of ``T``.
"""
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

@doc raw"""
    submodules(
      Tf::TorQuadModuleWithIsom;
      kwargs...,
    )

Given a torsion quadratic module with isometry ``(T, f)``, return an iterator
on the ``f``-stable submodules of ``T``. The possible keyword arguments to
restrict the submodules (currently) are:
- `quotype::Vector{Int}`: only submodules whose quotient are isomorphic as an
  abelian group to `abelian_group(quotype)`.

See also [`stable_submodules(::TorQuadModule, ::Vector{TorQuadModuleMap})`](@ref).
"""
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

@doc raw"""
    automorphism_group_with_inclusion(
      Tf::TorQuadModuleWithIsom
    ) -> AutomorphismGroup{TorQuadModule}, GAPGroupHomomorphism

Given a torsion quadratic module with isometry ``(T, f)``, return the
automorphism group ``O(T, f)`` of ``(T, f)``. This is defined as the
centralizer of the isometry ``f`` in the orthogonal group ``O(T)`` of ``T``.
The second output is the embedding of ``O(T, f)`` in ``O(T)``.

# Examples

```jldoctest
julia> Tf = torsion_quadratic_module_with_isometry(QQ[-1//60;], ZZ[11;])
Finite quadratic module of order 60
  with 1 generator
  with isometry given by
  [11]

julia> O, _ = automorphism_group_with_inclusion(Tf)
(Group of isometries of finite quadratic module: Z/60 -> Q/2Z, Hom: group of isometries -> group of isometries)

julia> order(O)
8
```
"""
@attr Tuple{AutomorphismGroup{TorQuadModule}, GAPGroupHomomorphism{AutomorphismGroup{TorQuadModule}, AutomorphismGroup{TorQuadModule}}
} function automorphism_group_with_inclusion(Tf::TorQuadModuleWithIsom)
  T = underlying_module(Tf)
  f = isometry(Tf)
  OT = orthogonal_group(T)
  fT = OT(f; check=false)
  return centralizer(OT, fT)
end

@doc raw"""
    automorphism_group(
      Tf::TorQuadModuleWithMap,
    ) -> AutomorphismGroup{TorQuadModule}

Return the automorphism group of `Tf`, see also
[`automorphism_group_with_inclusion(::TorQuadModuleWithIsom)`](@ref) for more
details.
"""
automorphism_group(Tf::TorQuadModuleWithIsom) = first(automorphism_group_with_inclusion(Tf))

@doc raw"""
    is_isomorphic_with_map(
      Tf::TorQuadModuleWithIsom,
      Sg::TorQuadModuleWithIsom,
    ) -> Bool, TorQuadModuleMap

Given two torsion quadratic modules with isometry ``(T, f)`` and ``(S, g)``,
return whether there exists an isometry ``ψ`` between ``T`` and ``S`` that
conjugates ``f`` and ``g``, i.e. such that ``g ∘ ψ = ψ ∘ f``.

If such a ``ψ`` exists, return ``true, ψ``, otherwise return ``false, 0``
where here ``0`` denotes the zero map from ``T`` to ``S``.
"""
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

  @hassert :ZZLatWithIsom 2 matrix(inv(phi) * f * phi) == matrix(isometry(Sg))
  return true, phi
end

@doc raw"""
    is_anti_isomorphic_with_map(
      Tf::TorQuadModuleWithIsom,
      Sg::TorQuadModuleWithIsom,
    ) -> Bool, TorQuadModuleMap

Given two torsion quadratic modules with isometry ``(T, f)`` and ``(S, g)``,
return whether there exists an anti-isometry ``ψ`` between ``T`` and ``S``
that conjugates ``f`` and ``g``, i.e. such that ``g ∘ ψ = ψ ∘ f``.

If such a ``ψ`` exists, return ``true, ψ``, otherwise return ``false, 0``
where here ``0`` denotes the zero map from ``T`` to ``S``.

# Examples
```jldoctest
julia> E8  = root_lattice(:E, 8);

julia> m = QQ[ 1  0  0  0  0  0  0  0;
               0  1  0  0  0  0  0  0;
              -2 -4 -5 -4 -3 -2 -1 -2;
               1  2  3  3  2  1  0  1;
              -1 -2 -3 -3 -2 -1  0 -2;
               1  2  3  2  1  0  0  2;
               0  0  0  0  0  1  0  0;
               2  4  6  5  4  3  2  3];

julia> Lf = integer_lattice_with_isometry(E8, m);

julia> F = kernel_lattice(Lf, 3);

julia> C = orthogonal_submodule(Lf, basis_matrix(F));

julia> Tf = discriminant_group(TorQuadModuleWithIsom, F)
Finite quadratic module of order 12
  with 2 generators
  with isometry given by
  [0   3]
  [1   1]

julia> Sg = discriminant_group(TorQuadModuleWithIsom, C)
Finite quadratic module of order 12
  with 2 generators
  with isometry given by
  [0   3]
  [1   1]

julia> is_anti_isomorphic_with_map(Tf, Sg)
(true, Map: finite quadratic module -> finite quadratic module)
```
"""
function is_anti_isomorphic_with_map(
  Tf::TorQuadModuleWithIsom,
  Sg::TorQuadModuleWithIsom,
)
  T = underlying_module(Tf)
  S = underlying_module(Sg)

  f = isometry(Tf)
  ok, phi = is_anti_isometric_with_anti_isometry(T, S)
  !ok && return false, phi

  fphi = inv(phi) * f * phi
  OS = orthogonal_group(S)
  gS = OS(isometry(Sg); check=false)
  fS = OS(fphi; check=false)
  ok, h = is_conjugate_with_data(OS, gS, fS)
  !ok && return false, hom(T, S, zero_matrix(ZZ, ngens(T), ngens(S)))

  phi = phi * hom(h)

  @hassert :ZZLatWithIsom 2 matrix(inv(phi) * f * phi) == matrix(isometry(Sg))
  return true, phi
end

###############################################################################
#
#  Equality and hash
#
###############################################################################

function Base.:(==)(T1::TorQuadModuleWithIsom, T2::TorQuadModuleWithIsom)
  underlying_module(T1) === underlying_module(T2) || return false
  return matrix(isometry(T1)) == matrix(isometry(T2))
end

function Base.hash(T::TorQuadModuleWithIsom, u::UInt)
  u = Base.hash(underlying_module(T), u)
  return Base.hash(matrix(isometry(T)), u)
end
