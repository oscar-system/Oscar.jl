#
# Some basic algorithmic invariant theory
# References:
# [S] B. Sturmfels, "Algorithms in invariant theory", 2nd ed., Springer, 2008.
# [DK] H,. Derksen and G. Kemper, "Computational Invariant Theory", 2nd ed., Springer, 2015.

module Invariants

using Oscar
import Base: ^, +, -, *


# TODO: for now, everything works with matrix groups, but eventually it
# might be better to work with representations resp. with G-modules
# instead


"""

action of a permutation on an MPolyElement by permuting the ring variables

```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> f = x^2 + 2y + 3
x^2 + 2*y + 3

julia> g = cperm(1:2)
(1,2)

julia> f^g
2*x + y^2 + 3
```
"""
function ^(f::MPolyElem, s::Oscar.PermGroupElem)
  G = parent(s)
  @assert ngens(parent(f)) == degree(G)

  g = Generic.MPolyBuildCtx(parent(f))
  for (c, e) = Base.Iterators.zip(Generic.MPolyCoeffs(f), Generic.MPolyExponentVectors(f))
    s_e = zeros(Int, degree(G))
    for i=1:degree(G)
      s_e[s(i)] = e[i]
    end
    push_term!(g, c, s_e)
  end
  return finish(g)
end


# action on an spoly
"""

action of a permutation on an spoly by permuting the ring variables

```jldoctest
julia> R, (x, y) = Singular.PolynomialRing(QQ, ["x", "y"])
(Singular Polynomial Ring (Coeffs(19)),(x,y),(dp(2),C), Singular.spoly{Singular.n_unknown{fmpq}}[x, y])

julia> f = x^2 + 2y + 3
x^2+2*y+3

julia> g = cperm(1:2)
(1,2)

julia> f^g
y^2+2*x+3
```
"""
function ^(f::Singular.spoly, g::Oscar.PermGroupElem)
  G = parent(g)
  @assert ngens(parent(f)) == degree(G)
  return Singular.permute_variables(f, Vector(g), parent(f))
end

# TODO: implement the above actions now for matrix groups over Q (later: over number fields)

# Reynolds operator ([S, Section 2.1])
# Naive implementation straight from the definition
reynolds(G::Oscar.GAPGroup, f::MPolyElem) = sum(g -> f^g, G) / order(G)


function minpoly(g::GAP.GapObj)
    # TODO: check whether g is a matrix / matrixobj
    # TODO: support ring argument?
    return GAP.Globals.MinimalPolynomial(g)
end

function charpoly(g::GAP.GapObj)
    # TODO: check whether g is a matrix / matrixobj
    # TODO: support ring argument
    return GAP.Globals.CharacteristicPolynomial(g)
end

function charpoly(g::Oscar.MatrixGroupElem)
    return charpoly(g.X)
end

#function charpoly(g::Oscar.PermGroupElem)
#  return charpoly(GAP.Globals.PermutationMat(g))
#end

#
# Theorem 2.2.1
# TODO: restrict this to matrix groups in char 0
function molien(G::Oscar.GAPGroup)
    cc = conjugacy_classes(G)
    # TODO: need to convert GAP matrix to Julia matrix; ensure the charpoly can be
    # inverted to get a rational function
    # TODO: for perm groups, compute the charpoly of the permutation
    s = sum(c -> length(c) / charpoly(representative(c)), cc)
    return s / order(G)
end


# TODO: Algorithm 2.2.5



end # module
