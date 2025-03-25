@doc raw"""
    QuadSpaceWithIsom

A container type for pairs $(V, f)$ consisting of a rational quadratic space
$V$ of type `QuadSpace` and an isometry $f$ given as a `QQMatrix`
representing the action on the standard basis of $V$.

We store the order of $f$ too, which can finite or infinite.

To construct an object of type *QuadSpaceWithIsom*, see the set of functions
called [`quadratic_space_with_isometry`](@ref)

# Examples
```jldoctest
julia> V = quadratic_space(QQ, 4);

julia> quadratic_space_with_isometry(V; neg=true)
Quadratic space of dimension 4
  with isometry of finite order 2
  given by
  [-1    0    0    0]
  [ 0   -1    0    0]
  [ 0    0   -1    0]
  [ 0    0    0   -1]

julia> L = root_lattice(:E, 6);

julia> V = ambient_space(L);

julia> f = matrix(QQ, 6, 6, [ 1  2  3  2  1  1;
                             -1 -2 -2 -2 -1 -1;
                              0  1  0  0  0  0;
                              1  0  0  0  0  0;
                             -1 -1 -1  0  0 -1;
                              0  0  1  1  0  1]);

julia> Vf = quadratic_space_with_isometry(V, f)
Quadratic space of dimension 6
  with isometry of finite order 8
  given by
  [ 1    2    3    2    1    1]
  [-1   -2   -2   -2   -1   -1]
  [ 0    1    0    0    0    0]
  [ 1    0    0    0    0    0]
  [-1   -1   -1    0    0   -1]
  [ 0    0    1    1    0    1]
```
"""
@attributes mutable struct QuadSpaceWithIsom
  V::Hecke.QuadSpace
  f::QQMatrix
  n::IntExt

  function QuadSpaceWithIsom(
      V::Hecke.QuadSpace,
      f::QQMatrix,
      n::IntExt
    )
    return new(V, f, n)
  end
end

@doc raw"""
    ZZLatWithIsom

A container type for pairs $(L, f)$ consisting of an integer lattice $L$ of
type `ZZLat` and an isometry $f$ given as a `QQMatrix` representing
the action on the basis matrix of $L$.

We store the ambient space $V$ of $L$ together with an isometry $f_a$
inducing $f$ on $L$ seen as a pair $(V, f_a)$ of type
[`QuadSpaceWithIsom`](@ref). We moreover store the order $n$ of $f$, which can
be finite or infinite.

To construct an object of type *ZZLatWithIsom*, see the following examples:

# Examples

One first way to construct such object, is by entering directly the lattice
with an isometry. The isometry can be a honnest isometry of the lattice, or
it can be an isometry of the ambient space preserving the lattice. Depending on
this choice, one should enter the appropriate boolean value
*ambient_representation*. This direct construction is done through the
constructors [`integer_lattice_with_isometry`](@ref).

```jldoctest
julia> L = root_lattice(:E, 6);

julia> f = matrix(QQ, 6, 6, [ 1  2  3  2  1  1;
                             -1 -2 -2 -2 -1 -1;                                                
                              0  1  0  0  0  0;      
                              1  0  0  0  0  0;
                             -1 -1 -1  0  0 -1;
                              0  0  1  1  0  1]);

julia> Lf = integer_lattice_with_isometry(L, f; ambient_representation=false)
Integer lattice of rank 6 and degree 6
  with isometry of finite order 8
  given by
  [ 1    2    3    2    1    1]
  [-1   -2   -2   -2   -1   -1]
  [ 0    1    0    0    0    0]
  [ 1    0    0    0    0    0]
  [-1   -1   -1    0    0   -1]
  [ 0    0    1    1    0    1]

julia> B = matrix(QQ,1,6, [1   2   3   1   -1   3]);

julia> I = lattice_in_same_ambient_space(L, B); # This is the invariant sublattice L^f

julia> If = integer_lattice_with_isometry(I, ambient_isometry(Lf))
Integer lattice of rank 1 and degree 6
  with isometry of finite order 1
  given by
  [1]

julia> integer_lattice_with_isometry(I; neg=true)
Integer lattice of rank 1 and degree 6
  with isometry of finite order 2
  given by
  [-1]
```

Another way to construct such objects is to see them as sub-objects of their
ambient space, of type [`QuadSpaceWithIsom`](@ref). Through the constructors
[`lattice(::QuadSpaceWithIsom)`](@ref) and
[`lattice_in_same_ambient_space(::ZZLatWithIsom, ::MatElem)`](@ref), one can
then construct lattices with isometry for free, in a given space, as long as
the module they define is preserved by the fixed isometry of the ambient space.

# Examples
```jldoctest
julia> G = matrix(QQ, 6, 6 , [ 3 1 -1 1 0 0;
                               1 3  1 1 1 1;
                              -1 1  3 0 0 1;
                               1 1  0 4 2 2;
                               0 1  0 2 4 2;
                               0 1  1 2 2 4]);

julia> V = quadratic_space(QQ, G);

julia> f = matrix(QQ, 6, 6, [ 1 0  0 0 0  0
                              0 0 -1 0 0  0
                             -1 1 -1 0 0  0
                              0 0  0 1 0 -1
                              0 0  0 0 0 -1
                              0 0  0 0 1 -1]);

julia> Vf = quadratic_space_with_isometry(V, f);

julia> Lf = lattice(Vf)
Integer lattice of rank 6 and degree 6
  with isometry of finite order 3
  given by
  [ 1   0    0   0   0    0]
  [ 0   0   -1   0   0    0]
  [-1   1   -1   0   0    0]
  [ 0   0    0   1   0   -1]
  [ 0   0    0   0   0   -1]
  [ 0   0    0   0   1   -1]

julia> B = matrix(QQ, 4, 6, [1 0 3 0 0 0;
                             0 1 1 0 0 0;
                             0 0 0 0 1 0;
                             0 0 0 0 0 1]);

julia> Cf = lattice(Vf, B)  # coinvariant sublattice L_f
Integer lattice of rank 4 and degree 6
  with isometry of finite order 3
  given by
  [-2   3   0    0]
  [-1   1   0    0]
  [ 0   0   0   -1]
  [ 0   0   1   -1]

julia> Cf2 = lattice_in_same_ambient_space(Lf, B)
Integer lattice of rank 4 and degree 6
  with isometry of finite order 3
  given by
  [-2   3   0    0]
  [-1   1   0    0]
  [ 0   0   0   -1]
  [ 0   0   1   -1]

julia> Cf == Cf2
true
```

The last equality of the last example shows why we care about
*ambient context*: the two pairs of lattice with isometry `Cf` and `Cf2` are
basically the same mathematical objects. Indeed, they lie in the same space,
defines the same module and their respective isometries are induced by the same
isometry of the ambient space. As for regular `ZZLat`, as soon as the lattices
are in the same ambient space, we can compare them as $\mathbb Z$-modules,
endowed with an isometry.
"""
@attributes mutable struct ZZLatWithIsom
  Vf::QuadSpaceWithIsom
  Lb::ZZLat
  f::QQMatrix
  n::IntExt

  function ZZLatWithIsom(
      Vf::QuadSpaceWithIsom,
      Lb::ZZLat,
      f::QQMatrix,
      n::IntExt
    )
    return new(Vf, Lb, f, n)
  end
end
