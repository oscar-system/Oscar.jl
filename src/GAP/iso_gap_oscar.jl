# Basically the same as the usual preimage function but without a type check
# since we don't have elem_type(D) in this case
function preimage(M::MapFromFunc{D, C}, a; check::Bool = true) where {D <: GapObj, C}
  parent(a) === codomain(M) || error("the element is not in the map's codomain")
  if isdefined(M.header, :preimage)
    p = M.header.preimage(a)
    return p
  end
  error("No preimage function known")
end

# needed in order to do a generic argument check on the GAP side
function image(M::MapFromFunc{D, C}, a; check::Bool = true) where {D <: GapObj, C}
  check && (a in domain(M) || error("the element is not in the map's domain"))
  if isdefined(M, :header)
    if isdefined(M.header, :image)
      return M.header.image(a)::elem_type(C)
    else
      error("No image function known")
    end
  else
    return M(a)
  end
end

################################################################################
#
#  Ring isomorphism
#
################################################################################

function _iso_gap_oscar_field_of_cyclotomics(F::GapObj)
   GAPWrap.IsCyclotomicField(F) && return _iso_gap_oscar_field_cyclotomic(F)
   F === GAP.Globals.Cyclotomics && return _iso_gap_oscar_abelian_closure(F)
   GAPWrap.DegreeOverPrimeField(F) == 2 && return _iso_gap_oscar_field_quadratic(F)
   GAPWrap.IsNumberField(F) && return _iso_gap_oscar_number_field(F)
   error("no method found")
end

function _iso_gap_oscar_residue_ring(RG::GapObj)
   n = GAPWrap.Size(RG)
   if n isa GapObj
     n = ZZRingElem(n)
   end
   RO = residue_ring(ZZ, n)[1]

   finv, f = _iso_oscar_gap_residue_ring_functions(RO, RG)

   return MapFromFunc(RG, RO, f, finv)
end

function _iso_gap_oscar_field_finite(FG::GapObj)
   FO = GF(characteristic(FG), GAPWrap.DegreeOverPrimeField(FG))

   finv, f = _iso_oscar_gap_field_finite_functions(FO, FG)

   return MapFromFunc(FG, FO, f, finv)
end

function _iso_gap_oscar_field_rationals(FG::GapObj)
   FO = QQ
   finv, f = _iso_oscar_gap_field_rationals_functions(FO, FG)

   return MapFromFunc(FG, FO, f, finv)
end

function _iso_gap_oscar_ring_integers(FG::GapObj)
   FO = ZZ
   finv, f = _iso_oscar_gap_ring_integers_functions(FO, FG)

   return MapFromFunc(FG, FO, f, finv)
end

function _iso_gap_oscar_field_cyclotomic(FG::GapObj)
   FO = cyclotomic_field(GAPWrap.Conductor(FG))[1]
   finv, f = _iso_oscar_gap_field_cyclotomic_functions(FO, FG)

   return MapFromFunc(FG, FO, f, finv)
end

function _iso_gap_oscar_field_quadratic(FG::GapObj)
   if GAPWrap.IsCyclotomicCollection(FG)
     # Determine the root that is contained in `FG`.
     N = GAPWrap.Conductor(FG)
     if mod(N, 4) == 3
       N = -N
     elseif mod(N, 8) == 4
       N = div(N, 4)
       if mod(N, 4 ) == 1
         N = -N
       end
     elseif mod(N, 8) == 0
       N = div(N, 4)
       FGgens = GAPWrap.GeneratorsOfField(FG)
       if any(x -> GAPWrap.GaloisCyc(x,-1) != x, FGgens)
         N = -N
       end
     end
     FO = quadratic_field(N)[1]
     finv, f = _iso_oscar_gap_field_quadratic_functions(FO, FG)
   elseif GAPWrap.IsAlgebraicExtension(FG)
     # For the moment, do not try to interpret an algebraic extension.
     return _iso_gap_oscar_number_field(FG)
   else
     error("do not know how to handle FG")
   end

   return MapFromFunc(FG, FO, f, finv)
end

# If `FG` is a number field that is not cyclotomic then
# it can be in `IsAlgebraicExtension` or in `IsCyclotomicCollection`
# or in `IsNumberFieldByMatrices`.
function _iso_gap_oscar_number_field(FG::GapObj)
   if GAPWrap.IsCyclotomicCollection(FG)
     # The abelian number fields of GAP store one generator.
     gensFG = GAPWrap.GeneratorsOfField(FG)
     length(gensFG) == 1 || error("do not know how to handle FG")
     z = gensFG[1]
     pol = GAPWrap.MinimalPolynomial(GAP.Globals.Rationals::GapObj, z)
     N = GAPWrap.DegreeOfLaurentPolynomial(pol)
     cfs = GAPWrap.CoefficientsOfUnivariatePolynomial(pol)
     R = Hecke.Globals.Qx
     polFO = R(Vector{QQFieldElem}(cfs))
     FO, _ = number_field(polFO, "z")
     powers = GapObj([z^i for i in 0:(N-1)])::GapObj
     B = GAPWrap.Basis(FG, powers)

     finv = function(x::Nemo.AbsSimpleNumFieldElem)
        coeffs = GapObj(coefficients(x); recursive = true)::GapObj
        return (coeffs * powers)::GAP.Obj
     end

     f = function(x::GAP.Obj)
        coeffs = Vector{QQFieldElem}(GAPWrap.Coefficients(B, x))
        return FO(coeffs)
     end
   elseif GAPWrap.IsAlgebraicExtension(FG)
     # This is the analogon of `_iso_oscar_gap(FO::AbsSimpleNumField)`.
     pol = GAPWrap.DefiningPolynomial(FG)
     cfs = GAPWrap.CoefficientsOfUnivariatePolynomial(pol)
     R = Hecke.Globals.Qx
     polFO = R(Vector{QQFieldElem}(cfs))
     FO, _ = number_field(polFO, "z")
     fam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(FG))

     finv = function(x::Nemo.AbsSimpleNumFieldElem)
        coeffs = GapObj(coefficients(x); recursive = true)::GapObj
        return GAPWrap.ObjByExtRep(fam, coeffs)
     end

     f = function(x::GAP.Obj)
        coeffs = Vector{QQFieldElem}(GAPWrap.ExtRepOfObj(x))
        return FO(coeffs)
     end
   elseif GAPWrap.IsNumberFieldByMatrices(FG)
     # `pol` is the minimal polynomial of `M` because the GAP code
     # caches the objects.
     M = GAPWrap.PrimitiveElement(FG)
     pol = GAPWrap.DefiningPolynomial(FG)

     # Construct the isomorphic number field in Oscar.
     cfs = GAPWrap.CoefficientsOfUnivariatePolynomial(pol)
     R = Hecke.Globals.Qx
     polFO = R(Vector{QQFieldElem}(cfs))
     FO, _ = number_field(polFO, "z", cached = false)

     # The canonical basis of `FG` does in general not consist of the
     # powers of `M`.
     # Thus we cache the powers of `M` in order to compute preimages
     # of elements from `FO`,
     # and we compute the decomposition of the elements of the canonical
     # basis of `FG` in terms of powers of `M` in order to compute images
     # in `FO`.
     pow = GAPWrap.One(M)
     Mpowers = [pow]
     for i in 2:degree(FO)
       pow = pow * M
       push!(Mpowers, pow)
     end

     C = GAPWrap.CanonicalBasis(FG)
     A = inv(GapObj([GAPWrap.Coefficients(C, m) for m in Mpowers]))
     Mpowers = GapObj(Mpowers)

     finv = function(x::Nemo.AbsSimpleNumFieldElem)
        coeffs = GapObj(coefficients(x); recursive = true)::GapObj
        return GAPWrap.LinearCombination(coeffs, Mpowers)
     end

     f = function(x::GAP.Obj)
        coeffs = GAPWrap.Coefficients(C, x)
        coeffs = Vector{QQFieldElem}(coeffs * A)
        return FO(coeffs)
     end
   else
     error("do not know how to handle FG")
   end

   return MapFromFunc(FG, FO, f, finv)
end

function _iso_gap_oscar_abelian_closure(FG::GapObj)
   FO, _ = abelian_closure(QQ)
   finv, f = _iso_oscar_gap_abelian_closure_functions(FO, FG)

   return MapFromFunc(FG, FO, f, finv)
end

function _iso_gap_oscar_univariate_polynomial_ring(RG::GapObj)
   coeffs_iso = iso_gap_oscar(GAPWrap.LeftActingDomain(RG))
   RO, x = polynomial_ring(codomain(coeffs_iso), :x, cached = false)
   finv, f = _iso_oscar_gap_polynomial_ring_functions(RO, RG, inv(coeffs_iso))

   return MapFromFunc(RG, RO, f, finv)
end

function _iso_gap_oscar_multivariate_polynomial_ring(RG::GapObj)
   coeffs_iso = iso_gap_oscar(GAPWrap.LeftActingDomain(RG))
   nams = [string(x) for x in GAPWrap.IndeterminatesOfPolynomialRing(RG)]
   RO, x = polynomial_ring(codomain(coeffs_iso), nams, cached = false)
   finv, f = _iso_oscar_gap_polynomial_ring_functions(RO, RG, inv(coeffs_iso))

   return MapFromFunc(RG, RO, f, finv)
end


"""
    Oscar.iso_gap_oscar(R) -> Map{GapObj, T}

Return an isomorphism `f` with `domain` the GAP object `R`
and `codomain` an Oscar object `S`.

Elements `x` of `R` are mapped to `S` via `f(x)`,
and elements `y` of `S` are mapped to `R` via `preimage(f, y)`.

Matrices `m` over `R` are mapped to matrices over `S` via
`map_entries(f, m)`,
and matrices `n` over `S` are mapped to matrices over `R` via
`Oscar.preimage_matrix(f, n)`.

Admissible values of `R` and the corresponding `S` are currently as follows.

| `S` (in `GAP.Globals`)             | `R`                        |
|:---------------------------------- |:-------------------------- |
| `Integers`                         | `ZZ`                       |
| `Rationals`                        | `QQ`                       |
| `mod(Integers, n)`                 | `residue_ring(ZZ, n)[1]`   |
| `GF(p, d)`                         | `finite_field(p, d)[1]`    |
| `CF(n)`                            | `cyclotomic_field(n)[1]`   |
| `AlgebraicExtension(Rationals, f)` | `number_field(g)[1]`       |
| `Cyclotomics`                      | `abelian_closure(QQ)[1]`   |
| `PolynomialRing(F)`                | `polynomial_ring(G)[1]`    |
| `PolynomialRing(F, n)`             | `polynomial_ring(G, n)[1]` |

(Here `g` is the polynomial over `QQ` that corresponds to the polynomial `f`,
and `G` is equal to `Oscar.iso_gap_oscar(F)`.)

# Examples
```jldoctest
julia> f = Oscar.iso_gap_oscar(GAP.Globals.Integers);

julia> x = ZZ(2)^100;  y = preimage(f, x)
GAP: 1267650600228229401496703205376

julia> f(y) == x
true

julia> m = matrix(ZZ, 2, 3, [1, 2, 3, 4, 5, 6]);

julia> n = Oscar.preimage_matrix(f, m)
GAP: [ [ 1, 2, 3 ], [ 4, 5, 6 ] ]

julia> map_entries(f, n) == m
true

julia> R = GAP.Globals.PolynomialRing(GAP.Globals.Rationals);

julia> f = Oscar.iso_gap_oscar(R);

julia> x = gen(codomain(f));

julia> pol = x^2 + x + 1;

julia> y = preimage(f, pol)
GAP: x_1^2+x_1+1

julia> f(y) == pol
true
```

!!! warning
    The functions `Oscar.iso_gap_oscar` and [`Oscar.iso_oscar_gap`](@ref)
    are not injective.
    Due to caching, it may happen that `S` stores an attribute value
    of `Oscar.iso_oscar_gap(S)`,
    but that the codomain of this map is not identical with
    or even not equal to the given `R`.
"""
iso_gap_oscar(F::GapObj) = GAP.Globals.IsoGapOscar(F)


# Compute the isomorphism between a GAP domain
# and a corresponding Oscar object.
# In order to admit adding new types of GAP objects,
# we let GAP's method selection decide which Julia function shall do the work.
# The methods must be installed in GAP at runtime,
# this will be done in the OscarInterface package.
const _iso_gap_oscar_methods = []

push!(_iso_gap_oscar_methods, "IsField and IsFinite" => _iso_gap_oscar_field_finite)
push!(_iso_gap_oscar_methods, "IsRationals" => _iso_gap_oscar_field_rationals)
push!(_iso_gap_oscar_methods, "IsField and IsCyclotomicCollection" => _iso_gap_oscar_field_of_cyclotomics)
push!(_iso_gap_oscar_methods, "IsField and IsAlgebraicExtension" => _iso_gap_oscar_number_field)
push!(_iso_gap_oscar_methods, "IsIntegers" => _iso_gap_oscar_ring_integers)
push!(_iso_gap_oscar_methods, "IsRing and IsZmodnZObjNonprimeCollection" => _iso_gap_oscar_residue_ring)
push!(_iso_gap_oscar_methods, "IsUnivariatePolynomialRing" => _iso_gap_oscar_univariate_polynomial_ring)
push!(_iso_gap_oscar_methods, "IsPolynomialRing" => _iso_gap_oscar_multivariate_polynomial_ring)

# Note that `IsNumberFieldByMatrices` is a GAP property,
# but the alnuth package uses it like a GAP category,
# i.e., one can rely on the fact that
# this filter is set in a field generated by matrices .
push!(_iso_gap_oscar_methods, "IsField and IsNumberFieldByMatrices" => _iso_gap_oscar_number_field)


################################################################################
#
#  Matrix space isomorphism
#
#  Using the known ring isomorphism from a GAP ring to an Oscar ring,
#  we can map matrices from GAP to Oscar using `map_entries`.
#  (The generic `map_entries` method cannot be used because the concepts of
#  `parent`and `_change_base_ring` do not fit to the situation in GAP.)
#  For the direction from Oscar to GAP, we introduce a generic function
#  `preimage_matrix` that takes the `ring_iso` and an Oscar matrix.
#
################################################################################

function preimage_matrix(f::Map{GapObj, T}, a::MatElem) where T
   isdefined(f.header, :preimage) || error("No preimage function known")
   @assert base_ring(a) === codomain(f)
   rows = Vector{GapObj}(undef, nrows(a))
   for i in 1:nrows(a)
      rows[i] = GapObj([preimage(f, a[i, j]) for j in 1:ncols(a)])
   end
   return GAPWrap.ImmutableMatrix(domain(f), GapObj(rows), true)
end

function AbstractAlgebra.map_entries(f::Map{GapObj, T}, a::GapObj) where T
   m = GAPWrap.NrRows(a)
   n = GAPWrap.NrCols(a)
   L = [f(a[i, j]) for i in 1:m for j in 1:n]
   return matrix(codomain(f), m, n, L)
end
