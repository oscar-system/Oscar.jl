@doc raw"""
    QuadSpaceWithIsom

A container type for pairs `(V, f)` consisting on an rational quadratic space
`V` of type `QuadSpace` and an isometry `f` given as a `QQMatrix` representing
the action on the standard basis of `V`.

We store the order of `f` too, which can finite or of infinite order.

To construct an object of type `QuadSpaceWithIsom`, see the set of functions
called [`quadratic_space_with_isometry`](@ref)
"""
@attributes mutable struct QuadSpaceWithIsom
  V::Hecke.QuadSpace
  f::QQMatrix
  n::IntExt

  function QuadSpaceWithIsom(V::Hecke.QuadSpace, f::QQMatrix, n::IntExt)
    z = new()
    z.V = V
    z.f = f
    z.n = n
    return z
  end
end

@doc raw"""
    ZZLatWithIsom

A container type for pairs `(L, f)` consisting on an integer lattice `L` of
type `ZZLat` and an isometry `f` given as a `QQMatrix` representing the action
on a given basis of `L`.

We store the ambient space `V` of `L` together with an isometry `f_ambient`
inducing `f` on `L` seen as a pair $(V, f_ambient)$ of type `QuadSpaceWithIsom`.
We moreover store the order `n` of `f`, which can be finite or infinite.

To construct an object of type `ZZLatWithIsom`, see the set of functions called
[`integer_lattice_with_isometry`](@ref)
"""
@attributes mutable struct ZZLatWithIsom
  Vf::QuadSpaceWithIsom
  Lb::ZZLat
  f::QQMatrix
  n::IntExt

  function ZZLatWithIsom(Vf::QuadSpaceWithIsom, Lb::ZZLat, f::QQMatrix, n::IntExt)
    z = new()
    z.Lb = Lb
    z.f = f
    z.Vf = Vf
    z.n = n
    return z
  end
end
