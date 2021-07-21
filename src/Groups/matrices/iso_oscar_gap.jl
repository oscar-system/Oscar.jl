import AbstractAlgebra: MatElem
import GAP: FFE

# Basically the same as the usual image function but without a type check since
# we don't have elem_type(C) in this case
function image(M::Map{D, C}, a) where {D, C <: GAP.GAP_jll.MPtr}
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
function ring_iso_oscar_gap(F::T) where T <: Union{Nemo.GaloisField, Nemo.GaloisFmpzField}
   p = characteristic(F)
   G = GAP.Globals.GF(GAP.julia_to_gap(p))

   f(x::Union{Nemo.gfp_elem, Nemo.gfp_fmpz_elem}) = GAP.julia_to_gap(lift(x))*GAP.Globals.One(G)

   # For "small" primes x should be an FFE, but for bigger ones it's GAP_jll.Mptr (?)
   finv(x) = F(GAP.gap_to_julia(GAP.Globals.IntFFE(x)))

   return MapFromFunc(f, finv, F, G)
end

function ring_iso_oscar_gap(F::T) where T <: Union{Nemo.FqNmodFiniteField, Nemo.FqFiniteField}
   p = characteristic(F)
   d = degree(F)
   Fp_gap = GAP.Globals.GF(GAP.julia_to_gap(p))

   if d == 1
      f1(x::Union{Nemo.fq_nmod, Nemo.fq}) = GAP.julia_to_gap(coeff(x, 0))*GAP.Globals.One(Fp_gap)

      # For "small" primes x should be an FFE, but for bigger ones it's GAP_jll.Mptr (?)
      finv1(x) = F(GAP.gap_to_julia(GAP.Globals.IntFFE(x)))

      return MapFromFunc(f1, finv1, F, Fp_gap)
   end

   e = GAP.Globals.One(Fp_gap)
   L = [ GAP.julia_to_gap(lift(coeff(defining_polynomial(F), i)))*e for i in 0:d - 1 ]
   push!(L, e)
   L_gap = GAP.julia_to_gap(L)
   f_gap = GAP.Globals.UnivariatePolynomial(Fp_gap, L_gap)
   G = GAP.Globals.GF(Fp_gap, f_gap)

   # TODO: Somehow Basis does not work for some (?) field. Should it or is there
   # a better way (on the GAP side)?
   # For example
   # > F, a = FiniteField(next_prime(2^8), 2, "a")
   # does not work.
   Basis_G = GAP.Globals.Basis(G)
   Basis_F = Vector{elem_type(F)}(undef, d)
   Basis_F[1] = F(1)
   for i = 2:d
      Basis_F[i] = Basis_F[i - 1]*gen(F)
   end

   function f(x::Union{fq_nmod, fq})
      v = [ GAP.julia_to_gap(coeff(x, i)) for i in 0:d - 1 ]
      return sum([ v[i]*Basis_G[i] for i in 1:d ])
   end

   # For "small" primes x should be an FFE, but for bigger ones it's GAP_jll.Mptr (?)
   function finv(x)
      v = GAP.Globals.Coefficients(Basis_G, x)
      # TODO: IntFFE apparently does not work for "big" primes, e.g.
      # > F, a = FiniteField(next_prime(fmpz(2)^64), 2, "a")
      # does not work.
      v_int = [ GAP.gap_to_julia(GAP.Globals.IntFFE(v[i])) for i = 1:length(v) ]
      return sum([ v_int[i]*Basis_F[i] for i = 1:d ])
   end

   return MapFromFunc(f, finv, F, G)
end

################################################################################
#
#  Matrix space isomorphism
#
################################################################################

# computes the isomorphism between the Oscar matrix space of dimension deg over
# F and the corresponding GAP matrix space

mat_iso_oscar_gap(F::Union{Nemo.GaloisField, Nemo.GaloisFmpzField, Nemo.FqNmodFiniteField, Nemo.FqFiniteField}, deg) = mat_iso_oscar_gap(F, deg, ring_iso_oscar_gap(F))

function mat_iso_oscar_gap(F::T, deg::Int, FtoGAP::MapFromFunc{T, S}) where {T <: Union{Nemo.GaloisField, Nemo.GaloisFmpzField, Nemo.FqNmodFiniteField, Nemo.FqFiniteField}, S}
   function f(x::MatElem)
      @assert base_ring(x) === domain(FtoGAP)
      rows = Vector{GapObj}(undef, nrows(x))
      for i in 1:nrows(x)
         rows[i] = GAP.julia_to_gap([ FtoGAP(x[i, j]) for j in 1:ncols(x) ])
      end

      return GAP.Globals.ImmutableMatrix(codomain(FtoGAP), GAP.julia_to_gap(rows), true)
   end

   function finv(x::GapObj)
      m = GAP.Globals.Size(x)# == nrows
      n = GAP.Globals.Size(x[1])# == ncols
      L = [ preimage(FtoGAP, x[i, j]) for i in 1:m for j in 1:n]

      return matrix(domain(FtoGAP), m, n, L)
   end

   M = MatrixSpace(F, deg, deg)
   return MapFromFunc(f, finv, M, GAP.Globals.MatrixAlgebra(codomain(FtoGAP), deg))
end
