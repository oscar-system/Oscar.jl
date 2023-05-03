@doc raw"""
    LatWithIsom

A container type for pairs `(L, f)` consisting on an integer lattice `L` of
type `ZLat` and an isometry `f` given as a `QQMatrix` representing the action
on a given basis of `L`.

The associated action `f_ambient` on the ambient space of `L` as well as the order
`n` of `f` are also stored.

To construct an object of type `LatWithIsom`, see the set of functions called
[`lattice_with_isometry`](@ref)
"""
@attributes mutable struct LatWithIsom
  Lb::ZLat
  f::QQMatrix
  f_ambient::QQMatrix
  n::IntExt

  function LatWithIsom(Lb::ZLat, f::QQMatrix, f_ambient::QQMatrix, n::IntExt)
    z = new()
    z.Lb = Lb
    z.f = f
    z.f_ambient = f_ambient
    z.n = n
    return z
  end
end

