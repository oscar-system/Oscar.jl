import AbstractAlgebra: MatElem
import GAP: FFE

# isomorphism from the Oscar field to the GAP field
# T: type of the domain
# S: element type of the domain (S == elem_type(T))
struct GenRingIso{T, S} <: Map{T,GapObj,SetMap,GenRingIso}
   domain::T      # Oscar field
   codomain::GapObj          # GAP field
   f                      # function from the Oscar field to the GAP field
   f_inv               # its inverse
end

# isomorphism from the Oscar matrix space to the GAP matrix space
# T: type of the domain
# S: element type base ring (!) of the domain, that is, the entries of the
#    matrices are of type S (T <: MatSpace{S})
struct GenMatIso{T, S} <: Map{T,GapObj,SetMap,GenMatIso}
   domain::T        # Oscar matrix space
   codomain::GapObj             # GAP matrix space
   fr::GenRingIso             # see type above
   f                         # function from the Oscar matrix space to the GAP matrix space
   f_inv                  # its inverse
end

################################################################################
#
#  Ring isomorphism
#
################################################################################

# computes the isomorphism between the Oscar field F and the corresponding GAP field F_gap
# the output has type GenRingIso (see above)

function gen_ring_iso(F::T) where T <: Union{Nemo.GaloisField, Nemo.GaloisFmpzField}
   p = characteristic(F)
   G = GAP.Globals.GF(GAP.julia_to_gap(p))

   f(x::Union{Nemo.gfp_elem, Nemo.gfp_fmpz_elem}) = GAP.julia_to_gap(lift(x))*GAP.Globals.One(G)

   # For "small" primes x should be an FFE, but for bigger ones it's GAP_jll.Mptr (?)
   finv(x) = F(GAP.gap_to_julia(GAP.Globals.IntFFE(x)))

   return GenRingIso{typeof(F), elem_type(F)}(F, G, f, finv)
end

function gen_ring_iso(F::T) where T <: Union{Nemo.FqNmodFiniteField, Nemo.FqFiniteField}
   p = characteristic(F)
   d = degree(F)
   Fp_gap = GAP.Globals.GF(GAP.julia_to_gap(p))

   if d == 1
      f1(x::Union{Nemo.fq_nmod, Nemo.fq}) = GAP.julia_to_gap(coeff(x, 0))*GAP.Globals.One(Fp_gap)

      # For "small" primes x should be an FFE, but for bigger ones it's GAP_jll.Mptr (?)
      finv1(x) = F(GAP.gap_to_julia(GAP.Globals.IntFFE(x)))

      return GenRingIso{T, elem_type(F)}(F, Fp_gap, f1, finv1)
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

   return GenRingIso{T, elem_type(F)}(F, G, f, finv)
end

function (g::GenRingIso{T, S})(x::S) where { T, S }
   @assert parent(x) === g.domain "Not an element of the domain"
   return g.f(x)
end

# TODO: Do we really want it like this? I does not make any sense intuitively
# (to me) to call the map itself for the preimage/inverse.
(g::GenRingIso)(x) = g.f_inv(x)

Base.show(io::IO, f::GenRingIso) = print(io, "Ring isomorphism between ", f.domain, " and the corresponding GAP field")

################################################################################
#
#  Matrix space isomorphism
#
################################################################################

# return the GAP matrix corresponding to the Oscar matrix x
# assumes x is a square matrix
function mat_oscar_gap(x::MatElem{S}, riso::GenRingIso{T, S}) where { T, S }
   rows = Vector{GapObj}(undef, nrows(x))
   for i in 1:nrows(x)
      rows[i] = GAP.julia_to_gap([ riso(x[i,j]) for j in 1:ncols(x) ])
   end

   return GAP.Globals.ImmutableMatrix(riso.codomain, GAP.julia_to_gap(rows), true)
end

# return the Oscar matrix corresponding to the GAP matrix x
# assumes x is a square matrix
function mat_gap_oscar(x::GapObj, riso::GenRingIso)
   m = GAP.Globals.Size(x)# == nrows
   n = GAP.Globals.Size(x[1])# == ncols
   L = [ riso(x[i, j]) for i in 1:m for j in 1:n]

   return matrix(riso.domain, m, n, L)
end

# computes the isomorphism between the Oscar matrix space of dimension deg over F and the corresponding GAP matrix space
# the output has type GenMatIso (see above)

function gen_mat_iso(deg::Int, F::Union{Nemo.GaloisField, Nemo.GaloisFmpzField, Nemo.FqNmodFiniteField, Nemo.FqFiniteField})
   riso = gen_ring_iso(F)                                      # "riso" = Ring ISOmorphism
   homom(x::MatElem) = mat_oscar_gap(x, riso)
   homominv(x::GapObj) = mat_gap_oscar(x, riso)
   M = MatrixSpace(F, deg, deg)
   return GenMatIso{typeof(M), elem_type(F)}(M, GAP.Globals.MatrixAlgebra(riso.codomain, deg), riso, homom, homominv)
end

function (g::GenMatIso{T, S})(x::MatElem{S}) where { T, S }
   @assert parent(x) === g.domain "Not an element of the domain"
   return g.f(x)
end

(g::GenMatIso)(x::GapObj) = g.f_inv(x)

Base.show(io::IO, f::GenMatIso) = print(io, "Matrix algebra homomorphism from Oscar algebra to GAP algebra")
