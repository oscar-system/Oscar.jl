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

# computes the isomorphism between the Oscar field F and the corresponding GAP field
function _iso_oscar_gap(F::Union{Nemo.GaloisField, Nemo.GaloisFmpzField})
   p = characteristic(F)
   G = GAPWrap.GF(GAP.Obj(p))
   e = GAPWrap.One(G)

   f(x::Union{Nemo.gfp_elem, Nemo.gfp_fmpz_elem}) = GAP.Obj(lift(x))*e
   function finv(x::GAP.Obj)
      y = GAPWrap.IntFFE(x)
      return y isa Int ? F(y) : F(fmpz(y))
   end

   return MapFromFunc(f, finv, F, G)
end

function _iso_oscar_gap(F::Union{Nemo.FqNmodFiniteField, Nemo.FqFiniteField})
   p = characteristic(F)
   d = degree(F)
   p_gap = GAP.Obj(p)
   Fp_gap = GAPWrap.GF(p_gap) # the prime field in GAP
   e = GAPWrap.One(Fp_gap)

   # prime fields are easy and efficient to deal with, handle them separately
   if d == 1
      f1(x::Union{Nemo.fq_nmod, Nemo.fq}) = GAP.Obj(coeff(x, 0))*e
      function finv1(x::GAP.Obj)
         y = GAPWrap.IntFFE(x)
         return y isa Int ? F(y) : F(fmpz(y))
      end
      return MapFromFunc(f1, finv1, F, Fp_gap)
   end

   # non-prime fields: convert the defining polynomial to GAP...
   L = [ GAP.Obj(lift(coeff(defining_polynomial(F), i)))*e for i in 0:d - 1 ]
   push!(L, e)
   L_gap = GapObj(L)
   f_gap = GAP.Globals.UnivariatePolynomial(Fp_gap, L_gap)

   # ... and compute a GAP field G defined via this polynomial
   # (If the given polynomial is a Conway polynomial then we may call
   # GAP's `GF(p, d)`, which avoids GAP's `AlgebraicExtension`.)
   if GAPWrap.IsCheapConwayPolynomial(p_gap, d) &&
      f_gap == GAP.Globals.ConwayPolynomial(p_gap, d)
      G = GAPWrap.GF(p_gap, d)
   else
      G = GAPWrap.GF(Fp_gap, f_gap)
   end

   # compute matching bases of both fields
   if GAPWrap.IsAlgebraicExtension(G)
#FIXME:
# As soon as the problem from https://github.com/gap-system/gap/issues/4694
# is fixed, go back to `Basis_G = GAP.Globals.Basis(G)` also in this case.
# Note that the above `GF` call delegates to `FieldExtension`,
# and we want a basis in the filter `IsCanonicalBasisAlgebraicExtension`.
      Basis_G = GAP.Globals.Objectify(GAP.Globals.NewType(
                GAP.Globals.FamilyObj(G),
                GAP.Globals.IsCanonicalBasisAlgebraicExtension),
                GAP.NewPrecord(0))
      GAP.Globals.SetUnderlyingLeftModule(Basis_G, G)
   else
      Basis_G = GAP.Globals.Basis(G)
   end
   Basis_F = Vector{elem_type(F)}(undef, d)
   Basis_F[1] = F(1)
   for i = 2:d
      Basis_F[i] = Basis_F[i - 1]*gen(F)
   end

   function f(x::Union{fq_nmod, fq})
      v = [ GAP.Obj(coeff(x, i)) for i in 0:d - 1 ]
      return sum([ v[i]*Basis_G[i] for i in 1:d ])
   end

   # For "small" primes x should be an FFE, but for bigger ones it's GAP_jll.Mptr (?)
   function finv(x::GAP.Obj)
      v = GAPWrap.Coefficients(Basis_G, x)
      v_int = [ fmpz(GAPWrap.IntFFE(v[i])) for i = 1:length(v) ]
      return sum([ v_int[i]*Basis_F[i] for i = 1:d ])
   end

   return MapFromFunc(f, finv, F, G)
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
function _iso_oscar_gap(F::AnticNumberField)

   flag, N = Hecke.iscyclotomic_type(F)
   flag || error("$F is not a cyclotomic field")

   G = GAPWrap.CF(GAP.Obj(N))
   cycpol = GAP.Globals.CyclotomicPol(N)
   dim = length(cycpol)-1

   function f(x::Nemo.nf_elem)
      coeffs = [Nemo.coeff(x, i) for i in 0:(N-1)]
      return GAPWrap.CycList(GAP.GapObj(coeffs; recursive=true))
   end

   function finv(x::GAP.Obj)
      GAPWrap.IsCyc(x) || error("$x is not a GAP cyclotomic")
      denom = GAPWrap.DenominatorCyc(x)
      n = GAPWrap.Conductor(x)
      mod(N, n) == 0 || error("$x does not lie in the $N-th cyclotomic field")
      coeffs = GAP.Globals.CoeffsCyc(x * denom, N)
      GAP.Globals.ReduceCoeffs(coeffs, cycpol)
      coeffs = Vector{fmpz}(coeffs)
      coeffs = coeffs[1:dim]
      denom = fmpz(denom)
      return F(coeffs) // denom
   end

   return MapFromFunc(f, finv, F, G)
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
