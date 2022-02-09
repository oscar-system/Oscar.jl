# Basically the same as the usual image function but without a type check since
# we don't have elem_type(C) in this case
function image(M::Map{D, C}, a) where {D, C <: GapObj}
  if isdefined(M, :header)
    if isdefined(M.header, :image)
      return M.header.image(a)
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

# Assume that `FO` and `FG` are finite fields of the same order
# in Oscar and GAP, respectively.
function _iso_oscar_gap_field_finite_functions(FO::Union{Nemo.GaloisField, Nemo.GaloisFmpzField}, FG::GAP.GapObj)
   p = characteristic(FO)
   e = GAPWrap.One(FG)

   f(x) = GAP.Obj(lift(x))*e

   finv = function(x::GAP.Obj)
     y = GAPWrap.IntFFE(x)
     return y isa Int ? FO(y) : FO(fmpz(y))
   end

   return (f, finv)
end

function _iso_oscar_gap_field_finite_functions(FO::Union{FqFiniteField, FqDefaultFiniteField, FqNmodFiniteField}, FG::GAP.GapObj)
   p = characteristic(FO)
   d = degree(FO)

   # Compute the canonical basis of `FG`.
   if ! GAP.Globals.IsPrimeField(GAP.Globals.LeftActingDomain(FG))
     # The GAP field is not an extension of the prime field.
     # What is a reasonable way to compute (on the GAP side) a polynomial
     # w.r.t. the prime field, and to decompose field elements w.r.t.
     # the corresponding basis?
     error("extensions of extension fields are not supported")
   elseif GAPWrap.IsAlgebraicExtension(FG)
#FIXME:
# As soon as the problem from https://github.com/gap-system/gap/issues/4694
# is fixed, go back to `basis_FG = GAP.Globals.Basis(FG)` also in this case.
# Note that we want a basis in the filter `IsCanonicalBasisAlgebraicExtension`.
     basis_FG = GAP.Globals.Objectify(GAP.Globals.NewType(
                GAP.Globals.FamilyObj(FG),
                GAP.Globals.IsCanonicalBasisAlgebraicExtension),
                GAP.NewPrecord(0))
     GAP.Globals.SetUnderlyingLeftModule(basis_FG, FG)
   else
     basis_FG = GAP.Globals.Basis(FG)
   end

   # Check whether the two fields have compatible polynomials.
   polFO = defining_polynomial(FO)
   coeffsFO = collect(coefficients(polFO))

   polFG = GAP.Globals.DefiningPolynomial(FG)
   coeffsFG = [fmpz(Oscar.GAPWrap.IntFFE(x)) for x in
               GAP.Globals.CoefficientsOfUnivariatePolynomial(polFG)]

   if coeffsFO == coeffsFG
     # The two fields are compatible.
     F = FO

     f = function(x)
       v = [GAP.Obj(coeff(x, i)) for i in 0:(d - 1)]
       return sum([v[i]*basis_FG[i] for i in 1:d])
     end

     finv = function(x::GAP.Obj)
       v = GAPWrap.Coefficients(basis_FG, x)
       v_int = [fmpz(GAPWrap.IntFFE(v[i])) for i = 1:d]
       return sum([v_int[i]*basis_F[i] for i = 1:d])
     end
   else
     # Create an Oscar field `FO2` that is compatible with `FG`
     # and has the same type as `FO` ...
     if FO isa FqNmodFiniteField
       p = Int(p)
     end
     R, x = PolynomialRing(GF(p), "x")
     FO2 = FiniteField(R(coeffsFG), "z")[1]

     # ... and an isomorphism between the two Oscar fields.
     emb = embed(FO2, FO)
     F = FO2

     f = function(x)
       v = [GAP.Obj(coeff(preimage(emb, x), i)) for i in 0:(d - 1)]
       return sum([v[i]*basis_FG[i] for i in 1:d])
     end

     finv = function(x::GAP.Obj)
       v = GAPWrap.Coefficients(basis_FG, x)
       v_int = [fmpz(GAPWrap.IntFFE(v[i])) for i = 1:d]
       return emb(sum([v_int[i]*basis_F[i] for i = 1:d]))
     end
   end

   # Compute the canonical basis of `FO` or `FO2`.
   basis_F = Vector{elem_type(F)}(undef, d)
   basis_F[1] = F(1)
   for i = 2:d
     basis_F[i] = basis_F[i - 1]*gen(F)
   end

   return (f, finv)
end

# Compute the isomorphism between the Oscar field `FO` and a corresponding
# GAP field.
# Try to avoid finite fields of the kind `IsAlgebraicExtension` on the GAP side,
# because they do not permit a lot of interesting computations.
# (Matrices over these fields are not suitable for MeatAxe functions,
# `FieldOfMatrixList` does not work, etc.)
# This means that we do not attempt to create a field on the GAP side
# whose defining polynomial fits to the one on the Oscar side;
# instead, we adjust this on the Oscar side if necessary.
# However, if an `IsAlgebraicExtension` field is the (currently) only supported
# field on the GAP side then we choose such a field,
# with the defining polynomial of the Oscar field.
function _iso_oscar_gap(FO::FinField)
   p = GAP.Obj(characteristic(FO))
   d = degree(FO)
   if GAP.Globals.IsCheapConwayPolynomial(p, d)
     FG = GAPWrap.GF(p, d)
   else
     # Calling `GAPWrap.GF(p, d)` would throw a GAP error.
     polFO = defining_polynomial(FO)
     coeffsFO = collect(coefficients(polFO))

     e = one(GAP.Globals.Z(p))
     fam = GAP.Globals.FamilyObj(e)
     coeffsFG = GAP.GapObj([GAP.Obj(lift(x))*e for x in coeffsFO])
     polFG = GAP.Globals.UnivariatePolynomialByCoefficients(fam, coeffsFG, 1)
     FG = GAPWrap.GF(p, polFG)
   end
   f, finv = _iso_oscar_gap_field_finite_functions(FO, FG)

   return MapFromFunc(f, finv, FO, FG)
end


@attributes FlintRationalField # TODO: port this to Nemo
function _iso_oscar_gap(F::FlintRationalField)
   return MapFromFunc(x -> GAP.Obj(x), x -> fmpq(x), F, GAP.Globals.Rationals)
end

@attributes FlintIntegerRing # TODO: port this to Nemo
function _iso_oscar_gap(F::FlintIntegerRing)
   return MapFromFunc(x -> GAP.Obj(x), x -> fmpz(x), F, GAP.Globals.Integers)
end

# For the moment, support only cyclotomic fields.
# General number fields will require caching,
# in order to achieve that Oscar matrix groups over the same number field
# are mapped to GAP matrix groups over the same `AlgebraicExtension` field.

# Assume that `FO` and `FG` are cyclotomic fields with the same conductor
# in Oscar and GAP, respectively.
function _iso_oscar_gap_field_cyclotomic_functions(FO::AnticNumberField, FG::GAP.GapObj)
   N = GAPWrap.Conductor(FG)
   cycpol = GAP.Globals.CyclotomicPol(N)
   dim = length(cycpol)-1

   f = function(x::Nemo.nf_elem)
      coeffs = [Nemo.coeff(x, i) for i in 0:(N-1)]
      return GAPWrap.CycList(GAP.GapObj(coeffs; recursive=true))
   end

   finv = function(x)
      GAPWrap.IsCyc(x) || error("$x is not a GAP cyclotomic")
      denom = GAPWrap.DenominatorCyc(x)
      n = GAPWrap.Conductor(x)
      mod(N, n) == 0 || error("$x does not lie in the $N-th cyclotomic field")
      coeffs = GAP.Globals.CoeffsCyc(x * denom, N)
      GAP.Globals.ReduceCoeffs(coeffs, cycpol)
      coeffs = Vector{fmpz}(coeffs)
      coeffs = coeffs[1:dim]
      return FO(coeffs) // fmpz(denom)
   end

   return (f, finv)
end

function _iso_oscar_gap(FO::AnticNumberField)
   flag, N = Hecke.iscyclotomic_type(FO)
   flag || error("$FO is not a cyclotomic field")

   FG = GAPWrap.CF(GAP.Obj(N))
   f, finv = _iso_oscar_gap_field_cyclotomic_functions(FO, FG)

   return MapFromFunc(f, finv, FO, FG)
end

function iso_oscar_gap(F)
   return get_attribute!(F, :iso_oscar_gap) do
      return _iso_oscar_gap(F)
   end
end

#TODO function iso_oscar_gap(F::T) where T <: QabField


################################################################################
#
#  Matrix space isomorphism
#
#  Using the known ring isomorphism from an Oscar ring to a GAP ring,
#  we can map matrices from Oscar to GAP using `map_entries`.
#  (The generic `map_entries` method cannot be used because the concepts of
#  `parent`and `_change_base_ring` do not fit to the situation in GAP.)
#  For the direction from GAP to Oscar, we introduce a generic function
#  `preimage_matrix` that takes the `ring_iso` and a GAP matrix.
#
################################################################################

function AbstractAlgebra.map_entries(f::Map{T, GapObj}, a::MatElem) where T
   isempty(a) && error("empty matrices are not supported by GAP")
   @assert base_ring(a) === domain(f)
   rows = Vector{GapObj}(undef, nrows(a))
   for i in 1:nrows(a)
      rows[i] = GapObj([f(a[i, j]) for j in 1:ncols(a)])
   end
   return GAP.Globals.ImmutableMatrix(codomain(f), GapObj(rows), true)
end

function preimage_matrix(f::Map{T, GapObj}, a::GapObj) where T
   isdefined(f.header, :preimage) || error("No preimage function known")
   m = GAPWrap.NrRows(a)
   n = GAPWrap.NrCols(a)
   L = [f.header.preimage(a[i, j]) for i in 1:m for j in 1:n]
   return matrix(domain(f), m, n, L)
end
