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
   if GAP.Globals.IsField(F)
     if GAPWrap.IsFinite(F)
       return _iso_gap_oscar_field_finite(F)
     else
       if GAP.Globals.IsRationals(F)
         return _iso_gap_oscar_field_rationals(F)
       elseif GAP.Globals.IsCyclotomicCollection(F)
         if GAP.Globals.IsCyclotomicField(F)
           return _iso_gap_oscar_field_cyclotomic(F)
         end
       end
     end
   elseif GAP.Globals.IsIntegers(F)
     return _iso_gap_oscar_ring_integers(F)
   elseif GAP.Globals.IsUnivariatePolynomialRing(F)
     return _iso_gap_oscar_univariate_polynomial_ring(F)
   end

   error("no method found")
end

function _iso_gap_oscar_field_finite(FG::GAP.GapObj)
   p = characteristic(FG)  # of type `fmpz`
   d = GAP.Globals.DegreeOverPrimeField(FG)
   if d == 1
     if p < fmpz(2)^64
       p = UInt64(p)
     end
     FO = GF(p)
   else
     FO = GF(p, d)
   end

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
   FO = CyclotomicField(GAPWrap.Conductor(FG))[1]
   finv, f = _iso_oscar_gap_field_cyclotomic_functions(FO, FG)

   return MapFromFunc(f, finv, FG, FO)
end

function _iso_gap_oscar_univariate_polynomial_ring(RG::GAP.GapObj)
   coeffs_iso = iso_gap_oscar(GAP.Globals.LeftActingDomain(RG))
   RO, x = PolynomialRing(codomain(coeffs_iso), "x")
   finv, f = _iso_oscar_gap_polynomial_ring_functions(RO, RG, inv(coeffs_iso))

   return MapFromFunc(f, finv, RG, RO)
end


# Use a GAP attribute for caching the mapping.
# The following must be executed at runtime,
# the function gets called in Oscar's `__init__`.
function __init_IsoGapOscar()
    if ! hasproperty(GAP.Globals, :IsoGapOscar)
      GAP.Globals.DeclareAttribute(GAP.Obj("IsoGapOscar"), GAP.Globals.IsDomain);
      GAP.Globals.InstallMethod(GAP.Globals.IsoGapOscar,
        GAP.Obj([GAP.Globals.IsDomain]), GAP.GapObj(_iso_gap_oscar));
    end
end

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
