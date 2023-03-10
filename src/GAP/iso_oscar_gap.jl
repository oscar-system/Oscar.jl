# Basically the same as the usual image function but without a type check since
# we don't have elem_type(C) in this case
function image(M::MapFromFunc{D, C}, a; check::Bool = true) where {D, C <: GapObj}
  parent(a) === domain(M) || error("the element is not in the map's domain")
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

# needed in order to do a generic argument check on the GAP side
function preimage(M::MapFromFunc{D, C}, a; check::Bool = true) where {D, C <: GapObj}
  if isdefined(M.header, :preimage)
    check && (a in codomain(M) || error("the element is not in the map's codomain"))
    p = M.header.preimage(a)::elem_type(D)
    @assert parent(p) === domain(M)
    return p
  end
  error("No preimage function known")
end

################################################################################
#
#  Ring isomorphism
#
################################################################################

# Assume that `RO` and `RG` are residue rings of the same size
# in Oscar and GAP, respectively.
function _iso_oscar_gap_residue_ring_functions(RO::Union{Nemo.zzModRing, Nemo.ZZModRing}, RG::GAP.GapObj)
   e = GAPWrap.One(RG)
   f(x) = GAP.Obj(lift(x))*e

   finv = function(x::GAP.Obj)
     @assert GAPWrap.IsFFE(x) || GAPWrap.IsZmodnZObj(x)
     y = GAP.Globals.Int(x)
     return y isa Int ? RO(y) : RO(ZZRingElem(y))
   end

   return (f, finv)
end

# Compute the isomorphism between the Oscar residue ring `RO`
# and a corresponding GAP residue ring.
function _iso_oscar_gap(RO::Union{Nemo.zzModRing, Nemo.ZZModRing})
   n = ZZRingElem(modulus(RO))
   RG = GAPWrap.mod(GAP.Globals.Integers::GapObj, GAP.Obj(n))
   f, finv = _iso_oscar_gap_residue_ring_functions(RO, RG)

   return MapFromFunc(f, finv, RO, RG)
end

# Assume that `FO` and `FG` are finite fields of the same order
# in Oscar and GAP, respectively.
function _iso_oscar_gap_field_finite_functions(FO::Union{Nemo.fpField, Nemo.FpField}, FG::GAP.GapObj)
   e = GAPWrap.One(FG)

   f(x) = GAP.Obj(lift(x))*e

   finv = function(x::GAP.Obj)
     y = GAPWrap.IntFFE(x)
     return y isa Int ? FO(y) : FO(ZZRingElem(y))
   end

   return (f, finv)
end

function _iso_oscar_gap_field_finite_functions(FO::Union{FqPolyRepField, FqField, fqPolyRepField}, FG::GAP.GapObj)
   p = characteristic(FO)
   d = degree(FO)

   # Compute the canonical basis of `FG`.
   if ! GAPWrap.IsPrimeField(GAPWrap.LeftActingDomain(FG))
     # The GAP field is not an extension of the prime field.
     # What is a reasonable way to compute (on the GAP side) a polynomial
     # w.r.t. the prime field, and to decompose field elements w.r.t.
     # the corresponding basis?
     error("extensions of extension fields are not supported")
   end
   basis_FG = GAPWrap.Basis(FG)
   # Test that we do not run into the problem from
   # https://github.com/gap-system/gap/issues/4694.
   @assert (!GAPWrap.IsAlgebraicExtension(FG)) ||
           GAPWrap.IsCanonicalBasisAlgebraicExtension(basis_FG)

   # Check whether the two fields have compatible polynomials.
   polFO = modulus(FO)
   coeffsFO = collect(coefficients(polFO))

   polFG = GAPWrap.DefiningPolynomial(FG)
   coeffsFG = [ZZRingElem(GAPWrap.IntFFE(x)) for x in
               GAPWrap.CoefficientsOfUnivariatePolynomial(polFG)]

   if coeffsFO == coeffsFG
     # The two fields are compatible.
     F = FO

     f = function(x)
       v = [GAP.Obj(coeff(x, i)) for i in 0:(d - 1)]
       return sum([v[i]*basis_FG[i] for i in 1:d])
     end

     finv = function(x::GAP.Obj)
       v = GAPWrap.Coefficients(basis_FG, x)
       v_int = [ZZRingElem(GAPWrap.IntFFE(v[i])) for i = 1:d]
       return sum([v_int[i]*basis_F[i] for i = 1:d])
     end
   else
     # Create an Oscar field `FO2` that is compatible with `FG`
     # and has the same type as `FO` ...
     R = parent(modulus(FO))
     FO2 = typeof(FO)(R(coeffsFG), :z, true; check = true)

     # ... and an isomorphism between the two Oscar fields.
     emb = embed(FO2, FO)
     F = FO2

     f = function(x)
       v = [GAP.Obj(coeff(preimage(emb, x), i)) for i in 0:(d - 1)]
       return sum([v[i]*basis_FG[i] for i in 1:d])
     end

     finv = function(x::GAP.Obj)
       v = GAPWrap.Coefficients(basis_FG, x)
       v_int = [ZZRingElem(GAPWrap.IntFFE(v[i])) for i = 1:d]
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
   p = GAP.Obj(characteristic(FO))::GAP.Obj
   d = degree(FO)
   if GAPWrap.IsCheapConwayPolynomial(p, d)
     FG = GAPWrap.GF(p, d)
   else
     # Calling `GAPWrap.GF(p, d)` would throw a GAP error.
     polFO = modulus(FO)
     coeffsFO = collect(coefficients(polFO))

     e = one(GAPWrap.Z(p))
     fam = GAPWrap.FamilyObj(e)
     coeffsFG = GAP.GapObj([GAP.Obj(lift(x))*e for x in coeffsFO])
     polFG = GAPWrap.UnivariatePolynomialByCoefficients(fam, coeffsFG, 1)
     FG = GAPWrap.GF(p, polFG)
   end
   f, finv = _iso_oscar_gap_field_finite_functions(FO, FG)

   return MapFromFunc(f, finv, FO, FG)
end


function _iso_oscar_gap_field_rationals_functions(FO::QQField, FG::GapObj)
#TODO   return (GAP.Obj, QQFieldElem)
   return (x -> GAP.Obj(x), x -> QQFieldElem(x))
end

function _iso_oscar_gap(FO::QQField)
   FG = GAP.Globals.Rationals::GapObj

   f, finv = _iso_oscar_gap_field_rationals_functions(FO, FG)

   return MapFromFunc(f, finv, FO, FG)
end

function _iso_oscar_gap_ring_integers_functions(FO::ZZRing, FG::GapObj)
#TODO  return (GAP.Obj, ZZRingElem)
   return (x -> GAP.Obj(x), x -> ZZRingElem(x))
end

function _iso_oscar_gap(FO::ZZRing)
   FG = GAP.Globals.Integers::GapObj

   f, finv = _iso_oscar_gap_ring_integers_functions(FO, FG)

   return MapFromFunc(f, finv, FO, FG)
end

# Assume that `FO` and `FG` are cyclotomic fields with the same conductor
# in Oscar and GAP, respectively.
# (Cyclotomic fields are easier to handle than general number fields.)
function _iso_oscar_gap_field_cyclotomic_functions(FO::AnticNumberField, FG::GAP.GapObj)
   N = conductor(FO)
   cycpol = GAPWrap.CyclotomicPol(N)
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
      coeffs = GAPWrap.CoeffsCyc(x * denom, N)
      GAPWrap.ReduceCoeffs(coeffs, cycpol)
      coeffs = Vector{ZZRingElem}(coeffs)
      coeffs = coeffs[1:dim]
      return FO(coeffs) // ZZRingElem(denom)
   end

   return (f, finv)
end

function _iso_oscar_gap(FO::AnticNumberField)
   flag, N = Hecke.is_cyclotomic_type(FO)
   if flag
     FG = GAPWrap.CF(GAP.Obj(N))
     f, finv = _iso_oscar_gap_field_cyclotomic_functions(FO, FG)
   else
     polFO = FO.pol
     N = degree(polFO)
     coeffs_polFO = collect(coefficients(polFO))
     fam = GAP.Globals.CyclotomicsFamily::GapObj
     cfs = GAP.GapObj(coeffs_polFO, recursive = true)::GapObj
     polFG = GAPWrap.UnivariatePolynomialByCoefficients(fam, cfs, 1)
     FG = GAPWrap.AlgebraicExtension(GAP.Globals.Rationals::GapObj, polFG)
     fam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(FG))

     f = function(x::Nemo.nf_elem)
        coeffs = GAP.GapObj(coefficients(x), recursive = true)::GapObj
        return GAPWrap.ObjByExtRep(fam, coeffs)
     end

     finv = function(x::GapObj)
        coeffs = Vector{QQFieldElem}(GAPWrap.ExtRepOfObj(x))
        return FO(coeffs)
     end
   end

   return MapFromFunc(f, finv, FO, FG)
end

# Assume that `FO` is a `QQAbField` and `FG` is `GAP.Globals.Cyclotomics`.
function _iso_oscar_gap_abelian_closure_functions(FO::QQAbField, FG::GAP.GapObj)
   return (GAP.julia_to_gap, QQAbElem)
end

function _iso_oscar_gap(FO::QQAbField)
   FG = GAP.Globals.Cyclotomics::GapObj
   f, finv = _iso_oscar_gap_abelian_closure_functions(FO, FG)

   return MapFromFunc(f, finv, FO, FG)
end

"""
    Oscar.iso_oscar_gap(R) -> Map{T, GAP.GapObj}

Return an isomorphism `f` with domain `R`
and `codomain` a GAP object `S`.

Elements `x` of `R` are mapped to `S` via `f(x)`,
and elements `y` of `S` are mapped to `R` via `preimage(f, y)`.

Matrices `m` over `R` are mapped to matrices over `S` via
`map_entries(f, m)`,
and matrices `n` over `S` are mapped to matrices over `R` via
`Oscar.preimage_matrix(f, n)`.

Admissible values of `R` and the corresponding `S` are currently as follows.

| `R`                                  | `S` (in `GAP.Globals`)             |
|:------------------------------------ |:---------------------------------- |
| `ZZ`                                 | `Integers`                         |
| `QQ`                                 | `Rationals`                        |
| `residue_ring(ZZ, n)`                | `mod(Integers, n)`                 |
| `FiniteField(p, d)[1]`               | `GF(p, d)`                         |
| `cyclotomic_field(n)[1]`             | `CF(n)`                            |
| `number_field(f::QQPolyRingElem)[1]` | `AlgebraicExtension(Rationals, g)` |
| `abelian_closure(QQ)[1]`             | `Cyclotomics`                      |
| `polynomial_ring(F)[1]`              | `PolynomialRing(G)`                |
| `polynomial_ring(F, n)[1]`           | `PolynomialRing(G, n)`             |

(Here `g` is the polynomial over `GAP.Globals.Rationals` that corresponds
to `f`,
and `G` is equal to `Oscar.iso_oscar_gap(F)`.)

# Examples
```jldoctest
julia> f = Oscar.iso_oscar_gap(ZZ);

julia> x = ZZ(2)^100;  y = f(x)
GAP: 1267650600228229401496703205376

julia> preimage(f, y) == x
true

julia> m = matrix(ZZ, 2, 3, [1, 2, 3, 4, 5, 6]);

julia> n = map_entries(f, m)
GAP: [ [ 1, 2, 3 ], [ 4, 5, 6 ] ]

julia> Oscar.preimage_matrix(f, n) == m
true

julia> R, x = polynomial_ring(QQ);

julia> f = Oscar.iso_oscar_gap(R);

julia> pol = x^2 + x - 1;

julia> y = f(pol)
GAP: x_1^2+x_1-1

julia> preimage(f, y) == pol
true
```
"""
@attr Map function iso_oscar_gap(F)
   return _iso_oscar_gap(F)
end


################################################################################
#
# Univariate polynomial rings
#
function _iso_oscar_gap_polynomial_ring_functions(RO::PolyRing{T}, RG::GAP.GapObj, coeffs_iso::MapFromFunc) where T
   fam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(codomain(coeffs_iso)))
   ind = GAPWrap.IndeterminateNumberOfUnivariateRationalFunction(
           GAPWrap.IndeterminatesOfPolynomialRing(RG)[1])

   f = function(x::PolyRingElem{T})
      cfs = GAP.GapObj([coeffs_iso(x) for x in coefficients(x)])
      return GAPWrap.UnivariatePolynomialByCoefficients(fam, cfs, ind)
   end

   finv = function(x)
      GAPWrap.IsPolynomial(x) || error("$x is not a GAP polynomial")
      cfs = Vector{GAP.Obj}(GAPWrap.CoefficientsOfUnivariatePolynomial(x))
      return RO([preimage(coeffs_iso, c) for c in cfs])
   end

   return (f, finv)
end

function _iso_oscar_gap(RO::PolyRing{T}) where T
   coeffs_iso = iso_oscar_gap(base_ring(RO))
   RG = GAPWrap.PolynomialRing(codomain(coeffs_iso))

   f, finv = _iso_oscar_gap_polynomial_ring_functions(RO, RG, coeffs_iso)

   return MapFromFunc(f, finv, RO, RG)
end


################################################################################
#
# Multivariate polynomial rings
#
@attributes AbstractAlgebra.Generic.MPolyRing # TODO: port this to AA

function _iso_oscar_gap_polynomial_ring_functions(RO::MPolyRing{T}, RG::GAP.GapObj, coeffs_iso::MapFromFunc) where T
   fam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(RG))
   n = nvars(RO)
   indets = GAPWrap.IndeterminatesOfPolynomialRing(RG)
   ind = [GAPWrap.IndeterminateNumberOfUnivariateRationalFunction(x)
          for x in indets]::Vector{Int}

   f = function(x::MPolyRingElem{T})
      extrep = []
      for (c, l) in zip(AbstractAlgebra.coefficients(x), AbstractAlgebra.exponent_vectors(x))
        v = []
        for i in 1:n
          if l[i] != 0
            append!(v, [i, l[i]])
          end
        end
        push!(extrep, GAP.GapObj(v))
        push!(extrep, coeffs_iso(c))
      end
      return GAPWrap.PolynomialByExtRep(fam, GAP.GapObj(extrep))
   end

   finv = function(x)
      GAPWrap.IsPolynomial(x) || error("$x is not a GAP polynomial")
      extrep = Vector{GAP.Obj}(GAPWrap.ExtRepPolynomialRatFun(x))
      M = Generic.MPolyBuildCtx(RO)
      for i in 1:2:length(extrep)
        v = fill(0, n)
        l = extrep[i]
        for j in 1:2:length(l)
          v[l[j]] = l[j+1]
        end
        push_term!(M, preimage(coeffs_iso, extrep[i+1]), v)
      end
      return finish(M)
   end

   return (f, finv)
end

function _iso_oscar_gap(RO::MPolyRing{T}) where T
   coeffs_iso = iso_oscar_gap(base_ring(RO))
   RG = GAPWrap.PolynomialRing(codomain(coeffs_iso), nvars(RO))

   f, finv = _iso_oscar_gap_polynomial_ring_functions(RO, RG, coeffs_iso)

   return MapFromFunc(f, finv, RO, RG)
end


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

function AbstractAlgebra.map_entries(f::Map{T, GapObj}, a::MatrixElem{S}) where {S <: RingElement, T}
   isempty(a) && error("empty matrices are not supported by GAP")
   @assert base_ring(a) === domain(f)
   rows = Vector{GapObj}(undef, nrows(a))
   for i in 1:nrows(a)
      rows[i] = GapObj([f(a[i, j]) for j in 1:ncols(a)])
   end
   return GAPWrap.ImmutableMatrix(codomain(f), GapObj(rows), true)
end

function preimage_matrix(f::Map{T, GapObj}, a::GapObj) where T
   isdefined(f.header, :preimage) || error("No preimage function known")
   m = GAPWrap.NrRows(a)
   n = GAPWrap.NrCols(a)
   L = [f.header.preimage(a[i, j]) for i in 1:m for j in 1:n]
   return matrix(domain(f), m, n, L)
end
