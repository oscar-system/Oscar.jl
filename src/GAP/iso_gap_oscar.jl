# Basically the same as the usual preimage function but without a type check
# since we don't have elem_type(C) in this case
function preimage(M::Map{D, C}, a) where {D <: GapObj, C}
  if isdefined(M.header, :preimage)
    p = M.header.preimage(a)
    return p
  end
  error("No preimage function known")
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
       if GAP.Globals.IsPrimeField(F)
         return _iso_gap_oscar_field_finite_prime(F)
       else
         return _iso_gap_oscar_field_finite_ext(F)
       end
     else
       if GAP.Globals.IsRationals(F)
         return _iso_gap_oscar_field_rationals(F)
       elseif GAP.Globals.IsCyclotomicCollection(F)
         if GAP.Globals.IsCyclotomicField(F)
           return _iso_gap_oscar_field_cyclotomic(F)
         end
       end
     end
   end

   error("no method found")
end

# For finite prime fields, choose `GF(p)`
# (of type `Nemo.GaloisField` if `p` is small (fits into `UInt64`),
# `Nemo.GaloisFmpzField` otherwise)).
function _iso_gap_oscar_field_finite_prime(F::GAP.GapObj)
   p = characteristic(F)  # of type `fmpz`
   if p < fmpz(2)^64
     p = UInt64(p)
   end
   G = GF(p)
   e = GAPWrap.One(F)

   f(x) = G(fmpz(GAPWrap.IntFFE(x)))
   finv(x::Union{Nemo.gfp_elem, Nemo.gfp_fmpz_elem}) = GAP.Obj(lift(x))*e

   return MapFromFunc(f, finv, F, G)
end

# For finite non-prime fields, choose `FiniteField(pol, "z")`
# (of type `FqNmodFiniteField` if `p` is small, `FqFiniteField` otherwise)
# because we can choose the same polynomial as is used in GAP.
function _iso_gap_oscar_field_finite_ext(F::GAP.GapObj)
   p = characteristic(F)  # of type `fmpz`
   if p < fmpz(2)^64
     p = UInt64(p)
   end
   d = GAP.Globals.DegreeOverPrimeField(F)
   e = GAPWrap.One(F)

   R, x = PolynomialRing(GF(p), "x")
   if GAP.Globals.IsPrimeField(GAP.Globals.LeftActingDomain(F))
     pol_gap = GAP.Globals.DefiningPolynomial(F)
     if GAPWrap.IsAlgebraicExtension(F)
#FIXME:
# As soon as the problem from https://github.com/gap-system/gap/issues/4694
# is fixed, go back to `Basis_F = GAP.Globals.Basis(F)` also in this case.
# Note that we want a basis in the filter `IsCanonicalBasisAlgebraicExtension`.
        Basis_F = GAP.Globals.Objectify(GAP.Globals.NewType(
                  GAP.Globals.FamilyObj(F),
                  GAP.Globals.IsCanonicalBasisAlgebraicExtension),
                  GAP.NewPrecord(0))
        GAP.Globals.SetUnderlyingLeftModule(Basis_F, F)
     else
        Basis_F = GAP.Globals.Basis(F)
     end
   else
     # The GAP field is not an extension of the prime field.
     # What is a reasonable way to compute (on the GAP side) a polynomial
     # w.r.t. the prime field, and to decompose field elements w.r.t.
     # the corresponding basis?
     error("extensions of extension fields are not supported")
   end

   coeffs = GAP.Globals.CoefficientsOfUnivariatePolynomial(pol_gap)
   coeffs = [fmpz(Oscar.GAPWrap.IntFFE(x)) for x in coeffs]
   pol = R(coeffs)
   G = FiniteField(pol, "z")[1]
   Basis_G = Vector{elem_type(G)}(undef, d)
   Basis_G[1] = G(1)
   for i = 2:d
      Basis_G[i] = Basis_G[i - 1]*gen(G)
   end

   function f(x)
      v = GAPWrap.Coefficients(Basis_F, x)
      v_int = [fmpz(GAPWrap.IntFFE(v[i])) for i = 1:length(v)]
      return sum([v_int[i]*Basis_G[i] for i = 1:d])
   end

   function finv(x::Union{fq_nmod, fq})
      v = [GAP.Obj(coeff(x, i)) for i in 0:(d - 1)]
      return sum([v[i]*Basis_F[i] for i in 1:d])
   end

   return MapFromFunc(f, finv, F, G)
end

function _iso_gap_oscar_field_rationals(F::GAP.GapObj)
   return MapFromFunc(x -> fmpq(x), x -> GAP.Obj(x), F, QQ)
end

function _iso_gap_oscar_field_cyclotomic(F::GAP.GapObj)
   N = GAPWrap.Conductor(F)
   cycpol = GAP.Globals.CyclotomicPol(N)
   dim = length(cycpol)-1
   G = CyclotomicField(N)[1]

   function f(x)
      GAPWrap.IsCyc(x) || error("$x is not a GAP cyclotomic")
      denom = GAPWrap.DenominatorCyc(x)
      n = GAPWrap.Conductor(x)
      mod(N, n) == 0 || error("$x does not lie in the $N-th cyclotomic field")
      coeffs = GAP.Globals.CoeffsCyc(x * denom, N)
      GAP.Globals.ReduceCoeffs(coeffs, cycpol)
      coeffs = Vector{fmpz}(coeffs)
      coeffs = coeffs[1:dim]
      return G(coeffs) // fmpz(denom)
   end

   function finv(x::Nemo.nf_elem)
      coeffs = [Nemo.coeff(x, i) for i in 0:(N-1)]
      return GAPWrap.CycList(GAP.GapObj(coeffs; recursive=true))
   end

   return MapFromFunc(f, finv, F, G)
end

# Use a GAP attribute for caching the mapping.
# The following must be executed at runtime,
# the function gets called in Oscar's `__init__`.
function __init_IsoGapOscar()
    if ! hasproperty(GAP.Globals, :IsoGapOscar)
      GAP.evalstr("""
DeclareAttribute( "IsoGapOscar", IsDomain );
InstallMethod( IsoGapOscar,
[ IsDomain ],
D -> Julia.Oscar._iso_gap_oscar( D ) );
""")
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
