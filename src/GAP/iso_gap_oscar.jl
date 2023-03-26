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

# Compute the isomorphism between the GAP domain `F`
# and a corresponding Oscar object.
function _iso_gap_oscar(F::GAP.GapObj)
   if GAPWrap.IsField(F)
     if GAPWrap.IsFinite(F)
       return _iso_gap_oscar_field_finite(F)
     else
       if GAPWrap.IsRationals(F)
         return _iso_gap_oscar_field_rationals(F)
       elseif GAPWrap.IsCyclotomicCollection(F)
         if GAPWrap.IsCyclotomicField(F)
           return _iso_gap_oscar_field_cyclotomic(F)
         elseif F === GAP.Globals.Cyclotomics
           return _iso_gap_oscar_abelian_closure(F)
         elseif GAPWrap.DegreeOverPrimeField(F) == 2
           return _iso_gap_oscar_field_quadratic(F)
         elseif GAPWrap.IsNumberField(F)
           return _iso_gap_oscar_number_field(F)
         end
       elseif GAPWrap.IsAlgebraicExtension(F)
         return _iso_gap_oscar_number_field(F)
       end
     end
   elseif GAPWrap.IsZmodnZObjNonprimeCollection(F)
     return _iso_gap_oscar_residue_ring(F)
   elseif GAPWrap.IsIntegers(F)
     return _iso_gap_oscar_ring_integers(F)
   elseif GAPWrap.IsUnivariatePolynomialRing(F)
     return _iso_gap_oscar_univariate_polynomial_ring(F)
   elseif GAPWrap.IsPolynomialRing(F)
     return _iso_gap_oscar_multivariate_polynomial_ring(F)
   end

   error("no method found")
end

function _iso_gap_oscar_residue_ring(RG::GAP.GapObj)
   n = GAPWrap.Size(RG)
   if n isa GAP.GapObj
     n = ZZRingElem(n)
   end
   RO = residue_ring(ZZ, n)

   finv, f = _iso_oscar_gap_residue_ring_functions(RO, RG)

   return MapFromFunc(f, finv, RG, RO)
end

function _iso_gap_oscar_field_finite(FG::GAP.GapObj)
   FO = Nemo._GF(characteristic(FG), GAPWrap.DegreeOverPrimeField(FG))

   finv, f = _iso_oscar_gap_field_finite_functions(FO, FG)

   return MapFromFunc(f, finv, FG, FO)
end

function _iso_gap_oscar_field_rationals(FG::GAP.GapObj)
   FO = QQ
   finv, f = _iso_oscar_gap_field_rationals_functions(FO, FG)

   return MapFromFunc(f, finv, FG, FO)
end

function _iso_gap_oscar_ring_integers(FG::GAP.GapObj)
   FO = ZZ
   finv, f = _iso_oscar_gap_ring_integers_functions(FO, FG)

   return MapFromFunc(f, finv, FG, FO)
end

function _iso_gap_oscar_field_cyclotomic(FG::GAP.GapObj)
   FO = cyclotomic_field(GAPWrap.Conductor(FG))[1]
   finv, f = _iso_oscar_gap_field_cyclotomic_functions(FO, FG)

   return MapFromFunc(f, finv, FG, FO)
end

function _iso_gap_oscar_field_quadratic(FG::GAP.GapObj)
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
       FGgens = GAP.Globals.GeneratorsOfField(FG)
       if any(x -> GAP.Globals.GaloisCyc(x,-1) != x, FGgens)
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

   return MapFromFunc(f, finv, FG, FO)
end

# If `FG` is a number field that is not cyclotomic then
# it can be in `IsAlgebraicExtension` or in `IsCyclotomicCollection`.
function _iso_gap_oscar_number_field(FG::GapObj)
   if GAPWrap.IsCyclotomicCollection(FG)
     # The abelian number fields of GAP store one generator.
     gensFG = GAPWrap.GeneratorsOfField(FG)
     length(gensFG) == 1 || error("do not know how to handle FG")
     z = gensFG[1]
     pol = GAPWrap.MinimalPolynomial(GAP.Globals.Rationals::GapObj, z)
     N = GAPWrap.DegreeOfLaurentPolynomial(pol)
     cfs = GAPWrap.CoefficientsOfUnivariatePolynomial(pol)
     R, x = polynomial_ring(QQ, "x")
     polFO = R(Vector{QQFieldElem}(cfs))
     FO, _ = number_field(polFO, "z")
     powers = GapObj([z^i for i in 0:(N-1)])::GapObj
     B = GAPWrap.Basis(FG, powers)

     finv = function(x::Nemo.nf_elem)
        coeffs = GAP.GapObj(coefficients(x), recursive = true)::GapObj
        return (coeffs * powers)::GAP.Obj
     end

     f = function(x::GAP.Obj)
        coeffs = Vector{QQFieldElem}(GAPWrap.Coefficients(B, x))
        return FO(coeffs)
     end
   elseif GAPWrap.IsAlgebraicExtension(FG)
     # This is the analogon of `_iso_oscar_gap(FO::AnticNumberField)`.
     pol = GAPWrap.DefiningPolynomial(FG)
     cfs = GAPWrap.CoefficientsOfUnivariatePolynomial(pol)
     R, x = polynomial_ring(QQ, "x")
     polFO = R(Vector{QQFieldElem}(cfs))
     FO, _ = number_field(polFO, "z")
     fam = GAPWrap.ElementsFamily(GAP.Globals.FamilyObj(FG))

     finv = function(x::Nemo.nf_elem)
        coeffs = GAP.GapObj(coefficients(x), recursive = true)::GapObj
        return GAPWrap.ObjByExtRep(fam, coeffs)
     end

     f = function(x::GAP.Obj)
        coeffs = Vector{QQFieldElem}(GAPWrap.ExtRepOfObj(x))
        return FO(coeffs)
     end
   else
     error("do not know how to handle FG")
   end

   return MapFromFunc(f, finv, FG, FO)
end

function _iso_gap_oscar_abelian_closure(FG::GAP.GapObj)
   FO, _ = abelian_closure(QQ)
   finv, f = _iso_oscar_gap_abelian_closure_functions(FO, FG)

   return MapFromFunc(f, finv, FG, FO)
end

function _iso_gap_oscar_univariate_polynomial_ring(RG::GAP.GapObj)
   coeffs_iso = iso_gap_oscar(GAPWrap.LeftActingDomain(RG))
   RO, x = polynomial_ring(codomain(coeffs_iso), "x")
   finv, f = _iso_oscar_gap_polynomial_ring_functions(RO, RG, inv(coeffs_iso))

   return MapFromFunc(f, finv, RG, RO)
end

function _iso_gap_oscar_multivariate_polynomial_ring(RG::GAP.GapObj)
   coeffs_iso = iso_gap_oscar(GAPWrap.LeftActingDomain(RG))
   nams = [string(x) for x in GAPWrap.IndeterminatesOfPolynomialRing(RG)]
   RO, x = polynomial_ring(codomain(coeffs_iso), nams)
   finv, f = _iso_oscar_gap_polynomial_ring_functions(RO, RG, inv(coeffs_iso))

   return MapFromFunc(f, finv, RG, RO)
end


"""
    Oscar.iso_gap_oscar(R) -> Map{GAP.GapObj, T}

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
| `mod(Integers, n)`                 | `residue_ring(ZZ, n)`      |
| `GF(p, d)`                         | `FiniteField(p, d)[1]`     |
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
"""
iso_gap_oscar(F::GAP.GapObj) = GAP.Globals.IsoGapOscar(F)


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
   return GAP.Globals.ImmutableMatrix(domain(f), GapObj(rows), true)
end

function AbstractAlgebra.map_entries(f::Map{GapObj, T}, a::GapObj) where T
   m = GAPWrap.NrRows(a)
   n = GAPWrap.NrCols(a)
   L = [f(a[i, j]) for i in 1:m for j in 1:n]
   return matrix(codomain(f), m, n, L)
end
