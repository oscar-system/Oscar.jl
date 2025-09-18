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
function _iso_oscar_gap_residue_ring_functions(RO::Union{zzModRing, ZZModRing}, RG::GapObj)
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
function _iso_oscar_gap(RO::Union{zzModRing, ZZModRing})
   n = ZZRingElem(modulus(RO))
   RG = GAPWrap.mod(GAP.Globals.Integers::GapObj, GAP.Obj(n))
   f, finv = _iso_oscar_gap_residue_ring_functions(RO, RG)

   return MapFromFunc(RO, RG, f, finv)
end

_ffe_to_int(a::FqFieldElem) = Nemo._coeff(a, 0)
_ffe_to_int(a::FqPolyRepFieldElem) = coeff(a, 0)
_ffe_to_int(a::fqPolyRepFieldElem) = coeff(a, 0)
_ffe_to_int(a::Union{fpFieldElem,FpFieldElem}) = lift(a)

function _make_prime_field_functions(FO, FG)
   e = GAPWrap.One(FG)

   f = function(x)
     y = GapObj(_ffe_to_int(x))::GapInt
     return y*e
   end

   finv = function(x::GAP.Obj)
     y = GAPWrap.IntFFE(x)
     return y isa Int ? FO(y) : FO(ZZRingElem(y))
   end

   return (f, finv)
end

# Assume that `FO` and `FG` are finite fields of the same order
# in Oscar and GAP, respectively.
function _iso_oscar_gap_field_finite_functions(FO::Union{fpField, FpField}, FG::GapObj)
   return _make_prime_field_functions(FO, FG)
end

function _iso_oscar_gap_field_finite_functions(FO::Union{FqPolyRepField, FqField, fqPolyRepField}, FG::GapObj)
   p = characteristic(FO)
   d = degree(FO)

   if degree(FO) != absolute_degree(FO) ||
      ! GAPWrap.IsPrimeField(GAPWrap.LeftActingDomain(FG))
     # The Oscar field or the GAP field is not an extension of the prime field.
     # What is a reasonable way to compute (on the GAP side) a polynomial
     # w.r.t. the prime field, and to decompose field elements w.r.t.
     # the corresponding basis?
     error("extensions of extension fields are not supported")
   end

   # handle prime fields first
   if d == 1
     return _make_prime_field_functions(FO, FG)
   end

   # Compute the canonical basis of `FG`.
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

     if FO isa FqField
       f = function(x)
         v = [GAP.Obj(Nemo._coeff(x, i)) for i in 0:(d - 1)]
         return sum([v[i]*basis_FG[i] for i in 1:d])
       end
     else
       f = function(x)
         v = [GAP.Obj(coeff(x, i)) for i in 0:(d - 1)]
         return sum([v[i]*basis_FG[i] for i in 1:d])
       end
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

     if FO isa FqField
       f = function(x)
         y = preimage(emb, x)
         v = [GAP.Obj(Nemo._coeff(y, i)) for i in 0:(d - 1)]
         return sum([v[i]*basis_FG[i] for i in 1:d])
       end
     else
       f = function(x)
         y = preimage(emb, x)
         v = [GAP.Obj(coeff(y, i)) for i in 0:(d - 1)]
         return sum([v[i]*basis_FG[i] for i in 1:d])
       end
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
   if d == 1 || GAPWrap.IsCheapConwayPolynomial(p, d)
     FG = GAPWrap.GF(p, d)
   else
     # Calling `GAPWrap.GF(p, d)` would throw a GAP error.
     polFO = modulus(FO)
     coeffsFO = collect(coefficients(polFO))

     e = one(GAPWrap.Z(p))
     fam = GAPWrap.FamilyObj(e)
     coeffsFG = GapObj([GAP.Obj(lift(x))*e for x in coeffsFO])
     polFG = GAPWrap.UnivariatePolynomialByCoefficients(fam, coeffsFG, 1)
     FG = GAPWrap.GF(p, polFG)
   end
   f, finv = _iso_oscar_gap_field_finite_functions(FO, FG)

   return MapFromFunc(FO, FG, f, finv)
end


function _iso_oscar_gap_field_rationals_functions(FO::QQField, FG::GapObj)
   return (GAP.Obj, QQFieldElem)
end

function _iso_oscar_gap(FO::QQField)
   FG = GAP.Globals.Rationals::GapObj

   f, finv = _iso_oscar_gap_field_rationals_functions(FO, FG)

   return MapFromFunc(FO, FG, f, finv)
end

function _iso_oscar_gap_ring_integers_functions(FO::ZZRing, FG::GapObj)
   return (GAP.Obj, ZZRingElem)
end

function _iso_oscar_gap(FO::ZZRing)
   FG = GAP.Globals.Integers::GapObj

   f, finv = _iso_oscar_gap_ring_integers_functions(FO, FG)

   return MapFromFunc(FO, FG, f, finv)
end

# Assume that `FO` and `FG` are cyclotomic fields with the same conductor
# in Oscar and GAP, respectively.
# (Cyclotomic fields are easier to handle than general number fields.)
function _iso_oscar_gap_field_cyclotomic_functions(FO::AbsSimpleNumField, FG::GapObj)
   flag, N = Hecke.is_cyclotomic_type(FO)
   @req flag "FO was not constructed as a cyclotomic field"
   cycpol = GAPWrap.CyclotomicPol(N)
   dim = length(cycpol)-1

   f = function(x::AbsSimpleNumFieldElem)
      coeffs = [Nemo.coeff(x, i) for i in 0:(N-1)]
      return GAPWrap.CycList(GapObj(coeffs; recursive = true))
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

# Assume that `FO` and `FG` are quadratic fields with the same square root
# in Oscar and GAP, respectively.
# (Quadratic fields are easier to handle than general number fields.)
function _iso_oscar_gap_field_quadratic_functions(FO::AbsSimpleNumField, FG::GapObj)
   flag, N = Hecke.is_quadratic_type(FO)
   @assert flag

   oO = one(FO)
   zO = gen(FO)
   oG = 1
   zG = GAPWrap.Sqrt(GAP.Obj(N))
   B = GAPWrap.BasisNC(FG, GapObj([oG, zG]))

   f = function(x::AbsSimpleNumFieldElem)
      return GAP.Obj(coeff(x,0)) * oG + GAP.Obj(coeff(x,1)) * zG
   end

   finv = function(x::GAP.Obj)
      GAPWrap.IsCyc(x) || error("$x is not a GAP cyclotomic")
      coeffs = GAPWrap.Coefficients(B, x)
      @req coeffs !== GAP.Globals.fail "$x is not an element oof $FG"
      return QQFieldElem(coeffs[1]) * oO + QQFieldElem(coeffs[2]) * zO
   end

   return (f, finv)
end

# Deal with simple extensions of Q.
function _iso_oscar_gap(FO::SimpleNumField{QQFieldElem})
   flag1, N1 = Hecke.is_cyclotomic_type(FO)
   flag2, N2 = Hecke.is_quadratic_type(FO)
   if flag1
     FG = GAPWrap.CF(GAP.Obj(N1))
     f, finv = _iso_oscar_gap_field_cyclotomic_functions(FO, FG)
   elseif flag2
     FG = GAPWrap.Field(GAPWrap.Sqrt(GAP.Obj(N2)))
     f, finv = _iso_oscar_gap_field_quadratic_functions(FO, FG)
   elseif degree(FO) == 1
     FG = GAP.Globals.Rationals::GapObj
     f = x -> GAP.Obj(coeff(x, 0))
     finv = x -> FO(QQ(x))
   else
     polFO = defining_polynomial(FO)
     coeffs_polFO = collect(coefficients(polFO))
     fam = GAP.Globals.CyclotomicsFamily::GapObj
     cfs = GapObj(coeffs_polFO; recursive = true)::GapObj
     polFG = GAPWrap.UnivariatePolynomialByCoefficients(fam, cfs, 1)
     FG = GAPWrap.AlgebraicExtension(GAP.Globals.Rationals::GapObj, polFG)
     fam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(FG))

     f = function(x::SimpleNumFieldElem{QQFieldElem})
        coeffs = GapObj(coefficients(x); recursive = true)::GapObj
        return GAPWrap.AlgExtElm(fam, coeffs)
     end

     finv = function(x::GapObj)
        coeffs = Vector{QQFieldElem}(GAPWrap.ExtRepOfObj(x))
        return FO(coeffs)
     end
   end

   return MapFromFunc(FO, FG, f, finv)
end

# Deal with simple extensions of proper extensions of Q.
function _iso_oscar_gap(FO::SimpleNumField{T}) where T <: FieldElem
   B = base_field(FO)
   isoB = iso_oscar_gap(B)
   BG = codomain(isoB)::GapObj

   polFO = defining_polynomial(FO)
   coeffs_polFO = collect(coefficients(polFO))
   fam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(BG))
   cfs = GapObj([isoB(x) for x in coeffs_polFO])::GapObj
   polFG = GAPWrap.UnivariatePolynomialByCoefficients(fam, cfs, 1)
   FG = GAPWrap.AlgebraicExtension(BG, polFG)
   fam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(FG))

   f = function(x::SimpleNumFieldElem{T})
      coeffs = GapObj([isoB(x) for x in coefficients(x)])::GapObj
      return GAPWrap.AlgExtElm(fam, coeffs)
   end

   finv = function(x::GapObj)
      coeffs = [preimage(isoB, x) for x in GapObj(GAPWrap.ExtRepOfObj(x))]
      return FO(coeffs)
   end

   return MapFromFunc(FO, FG, f, finv)
end

# Deal with non-simple extensions of Q or of extensions of Q.
function _iso_oscar_gap(FO::NumField)
   @assert ! is_simple(FO)
   if is_absolute(FO)
     F, emb = absolute_simple_field(FO)
   else
     F, emb = simple_extension(FO)
   end
   iso = iso_oscar_gap(F)
   FG = codomain(iso)
   fam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(FG))
   B = base_field(F)
   isoB = iso_oscar_gap(B)

   f = function(x::NumFieldElem)
      coeffs = GapObj([isoB(x) for x in coefficients(preimage(emb, x))])::GapObj
      return GAPWrap.AlgExtElm(fam, coeffs)
   end

   if is_absolute(FO)
     finv = function(x::GapObj)
        coeffs = Vector{QQFieldElem}(GAPWrap.ExtRepOfObj(x))
        return emb(F(coeffs))
     end
   else
     finv = function(x::GapObj)
        coeffs = [preimage(isoB, y) for y in GAPWrap.ExtRepOfObj(x)]
        return emb(F(coeffs))
     end
   end

   return MapFromFunc(FO, FG, f, finv)
end


# Assume that `FO` is a `QQAbField` and `FG` is `GAP.Globals.Cyclotomics`.
function _iso_oscar_gap_abelian_closure_functions(FO::QQAbField, FG::GapObj)
   return (GapObj, QQAbFieldElem)
end

function _iso_oscar_gap(FO::QQAbField)
   FG = GAP.Globals.Cyclotomics::GapObj
   f, finv = _iso_oscar_gap_abelian_closure_functions(FO, FG)

   return MapFromFunc(FO, FG, f, finv)
end

"""
    Oscar.iso_oscar_gap(R::T) -> Map{T, GapObj}

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
| `residue_ring(ZZ, n)[1]`             | `mod(Integers, n)`                 |
| `finite_field(p, d)[1]`              | `GF(p, d)`                        |
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

!!! warning
    The functions `Oscar.iso_oscar_gap` and [`Oscar.iso_gap_oscar`](@ref)
    are not injective.
    Due to caching, it may happen that `S` stores an attribute value
    of `Oscar.iso_gap_oscar(S)`,
    but that the codomain of this map is not identical with
    or even not equal to the given `R`.

    Note also that `R` and `S` may differ w.r.t. some structural properties
    because GAP does not support all kinds of constructions that are
    possible in Oscar.
    For example, if `R` is a non-simple number field then `S` will be a
    simple extension because GAP knows only simple field extensions.
    Thus using `Oscar.iso_oscar_gap(R)` for objects `R` whose recursive
    structure is not fully supported in GAP will likely cause overhead
    at runtime.
"""
@attr Map{T, GapObj} function iso_oscar_gap(F::T) where T
   return _iso_oscar_gap(F)
end


################################################################################
#
# Univariate polynomial rings
#
function _iso_oscar_gap_polynomial_ring_functions(RO::PolyRing{T}, RG::GapObj, coeffs_iso::MapFromFunc) where T
   fam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(codomain(coeffs_iso)))
   ind = GAPWrap.IndeterminateNumberOfUnivariateRationalFunction(
           GAPWrap.IndeterminatesOfPolynomialRing(RG)[1])

   f = function(x::PolyRingElem{T})
      cfs = GapObj([coeffs_iso(x) for x in coefficients(x)])
      return GAPWrap.UnivariatePolynomialByCoefficients(fam, cfs, ind)
   end

   finv = function(x)
      GAPWrap.IsPolynomial(x) || error("$x is not a GAP polynomial")
      cfs = Vector{GAP.Obj}(GAPWrap.CoefficientsOfUnivariatePolynomial(x))
      return RO([preimage(coeffs_iso, c) for c in cfs])
   end

   return (f, finv)
end

function _iso_oscar_gap(RO::PolyRing)
   coeffs_iso = iso_oscar_gap(base_ring(RO))
   RG = GAPWrap.PolynomialRing(codomain(coeffs_iso))

   f, finv = _iso_oscar_gap_polynomial_ring_functions(RO, RG, coeffs_iso)

   return MapFromFunc(RO, RG, f, finv)
end


################################################################################
#
# Multivariate polynomial rings
#
function _iso_oscar_gap_polynomial_ring_functions(RO::MPolyRing{T}, RG::GapObj, coeffs_iso::MapFromFunc) where T
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
        push!(extrep, GapObj(v))
        push!(extrep, coeffs_iso(c))
      end
      return GAPWrap.PolynomialByExtRep(fam, GapObj(extrep))
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

   return MapFromFunc(RO, RG, f, finv)
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
      rows[i] = GapObj([f(a[i, j])::GAP.Obj for j in 1:ncols(a)])
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
