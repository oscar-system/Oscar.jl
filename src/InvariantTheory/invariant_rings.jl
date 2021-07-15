
mutable struct InvRing{S, T, U, V, W, X}
   field::S
   poly_ring::T

   group::U
   action::Vector{V}          # Needs to be set (so far)
   action_singular::Vector{W} # Needs to be set (so far)

   modular::Bool

   primary::Vector{X}
   secondary::Vector{X}
   irreducible_secondary::Vector{X}
   fundamental::Vector{X}

   # Cache some stuff on the Singular side
   # (possibly removed at some point)
   reynolds_singular::Singular.smatrix
   molien_singular::Singular.smatrix
   primary_singular::Singular.smatrix

   function InvRing(K::S, R::T, G::U, action::Vector{V}) where {S, T, U, V}
     n = degree(G)
     R, = PolynomialRing(K, "x" => 1:n)
     R_sing = singular_ring(R)
     action_singular = identity.([change_base_ring(R_sing, g) for g in action])
     X = elem_type(R)
     W = eltype(action_singular)
     z = new{S, T, U, V, W, X}()
     z.field = K
     z.poly_ring = R
     z.group = G
     z.action = action
     z.action_singular = action_singular
     z.modular = true
     if iszero(characteristic(K))
       z.modular = false
     else
       if !iszero(mod(order(G), characteristic(K)))
         z.modular = false
       end
     end
     return z
   end
end

################################################################################
#
#  Field access
#
################################################################################

coefficient_ring(I::InvRing) = I.field

polynomial_ring(I::InvRing) = I.poly_ring

action(I::InvRing) = I.action

_action_singular(I::InvRing) = I.action_singular

group(I::InvRing) = I.group

ismodular(I::InvRing) = I.modular

function invariant_ring(M::Vector{<: MatrixElem})
  return invariant_ring(base_ring(M[1]), M)
end

function invariant_ring(K::Field, M::Vector{<: MatrixElem})
  return invariant_ring(matrix_group([change_base_ring(K, g) for g in M]))
end

function invariant_ring(G::MatrixGroup)
  n = degree(G)
  R, = PolynomialRing(base_ring(G), "x" => 1:n, cached = false)
  S = singular_ring(R)
  action = mat_elem_type(typeof(G))[g.elm for g in gens(G)]
  return InvRing(base_ring(G), R, G, action)
end

invariant_ring(matrices::MatrixElem{T}...) where {T} = invariant_ring(collect(matrices))

function Base.show(io::IO, IR::InvRing)
  print(io, "Invariant ring of\n")
  print(io, group(IR), "\n")
  print(io, "with generators\n")
  print(io, action(IR))
end

function reynolds_molien_via_singular(IR::InvRing)
   if !isdefined(IR, :reynolds_singular) || !isdefined(IR, :molien_singular)
      R = polynomial_ring(IR)
      singular_matrices = _action_singular(IR)
      rey, mol = Singular.LibFinvar.reynolds_molien(singular_matrices...)
      IR.reynolds_singular = rey
      IR.molien_singular = mol
   end
   return IR.reynolds_singular, IR.molien_singular
end

function primary_invariants_via_singular(IR::InvRing)
   if !isdefined(IR, :primary_singular)
      if ismodular(IR)
         rey, mol = reynolds_molien_via_singular(IR)
         P = Singular.LibFinvar.primary_char0(rey, mol)
         R = polynomial_ring(IR)
         p = Vector{elem_type(R)}()
         for i = 1:ncols(P)
            push!(p, R(P[1, i]))
         end
         IR.primary_singular = P
         IR.primary = p
      else
         P = Singular.LibFinvar.primary_invariants(_action_singular(IR)...)[1]
         R = polynomial_ring(IR)
         p = Vector{elem_type(R)}()
         for i = 1:ncols(P)
           push!(p, R(P[1, i]))
         end
         IR.primary_singular = P
         IR.primary = p
      end
   end
   return IR.primary
end

function primary_invariants(IR::InvRing)
   if !isdefined(IR, :primary)
      return primary_invariants_via_singular(IR)
   end
   return copy(IR.primary)
end

function secondary_invariants_via_singular(IR::InvRing)
   if !isdefined(IR, :secondary)
      rey, mol = reynolds_molien_via_singular(IR)
      primary_invariants_via_singular(IR)
      P = IR.primary_singular
      S, IS = Singular.LibFinvar.secondary_char0(P, rey, mol)
      R = polynomial_ring(IR)
      s = Vector{elem_type(R)}()
      for i = 1:ncols(S)
         push!(s, R(S[1, i]))
      end
      is = Vector{elem_type(R)}()
      for i = 1:ncols(IS)
         push!(is, R(IS[1, i]))
      end
      IR.secondary = s
      IR.irreducible_secondary = is
   end
   return IR.secondary
end

function secondary_invariants(IR::InvRing)
   if !isdefined(IR, :secondary)
      return secondary_invariants_via_singular(IR)
   end
   return IR.secondary
end

function irreducible_secondary_invariants(IR::InvRing)
   if !isdefined(IR, :irreducible_secondary)
      secondary_invariants_via_singular(IR)
   end
   return IR.irreducible_secondary
end

# Doesn't belong here...
# Matrices act from the left here!
function heisenberg_group(n::Int)
   K, a = CyclotomicField(n, "a")
   M1 = zero_matrix(K, n, n)
   M1[1, n] = one(K)
   for i = 2:n
      M1[i, i - 1] = one(K)
   end

   M2 = zero_matrix(K, n, n)
   M2[1, 1] = one(K)
   for i = 2:n
      M2[i, i] = M2[i - 1, i - 1]*a
   end
   return MatrixGroup(n, K, [ M1, M2 ])
end
