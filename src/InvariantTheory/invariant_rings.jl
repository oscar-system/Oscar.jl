export invariant_ring, primary_invariants, secondary_invariants

###############################################

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

   function InvRing(K::S, G::U, action::Vector{V}) where {S, U, V}
     n = degree(G)
     R, = PolynomialRing(K, "x" => 1:n, cached = false)
     R_sing = singular_ring(R)
     action_singular = identity.([change_base_ring(R_sing, g) for g in action])
     T = typeof(R)
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

#######################################################

@doc Markdown.doc"""
    invariant_ring(G::MatrixGroup)

Return the invariant ring of the finite matrix group `G`.
    
CAVEAT: The creation of invariant rings is lazy in the sense that no explicit computations are done until specifically invoked (for example by the `primary_invariants` function).

# Examples
```jldoctest
julia> K, a = CyclotomicField(3, "a")
(Cyclotomic field of order 3, a)

julia> M1 = matrix(K, [0 0 1; 1 0 0; 0 1 0])
[0   0   1]
[1   0   0]
[0   1   0]

julia> M2 = matrix(K, [1 0 0; 0 a 0; 0 0 -a-1])
[1   0        0]
[0   a        0]
[0   0   -a - 1]

julia> G = MatrixGroup(3, K, [ M1, M2 ])
Matrix group of degree 3 over Cyclotomic field of order 3

julia> IR = invariant_ring(G)
Invariant ring of
Matrix group of degree 3 over Cyclotomic field of order 3
with generators
AbstractAlgebra.Generic.MatSpaceElem{nf_elem}[[0 0 1; 1 0 0; 0 1 0], [1 0 0; 0 a 0; 0 0 -a-1]]
```
"""
function invariant_ring(G::MatrixGroup)
  n = degree(G)
  action = mat_elem_type(typeof(G))[g.elm for g in gens(G)]
  return InvRing(base_ring(G), G, action)
end

#######################################################

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
      if iszero(characteristic(coefficient_ring(IR)))
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

#######################################################

@doc Markdown.doc"""
    primary_invariants(IR::InvRing)

Return a system of primary invariants of `IR`.

If a system of primary invariants of `IR` is already cached, return the cached system. 
Otherwise, compute and cache such a system first.

NOTE: The primary invariants are sorted by increasing degree.

# Examples
```jldoctest
julia> K, a = CyclotomicField(3, "a")
(Cyclotomic field of order 3, a)

julia> M1 = matrix(K, [0 0 1; 1 0 0; 0 1 0])
[0   0   1]
[1   0   0]
[0   1   0]

julia> M2 = matrix(K, [1 0 0; 0 a 0; 0 0 -a-1])
[1   0        0]
[0   a        0]
[0   0   -a - 1]

julia> G = MatrixGroup(3, K, [ M1, M2 ])
Matrix group of degree 3 over Cyclotomic field of order 3

julia> IR = invariant_ring(G)
Invariant ring of
Matrix group of degree 3 over Cyclotomic field of order 3
with generators
AbstractAlgebra.Generic.MatSpaceElem{nf_elem}[[0 0 1; 1 0 0; 0 1 0], [1 0 0; 0 a 0; 0 0 -a-1]]

julia> primary_invariants(IR)
3-element Vector{AbstractAlgebra.Generic.MPoly{nf_elem}}:
 x[1]*x[2]*x[3]
 x[1]^3 + x[2]^3 + x[3]^3
 x[1]^3*x[2]^3 + x[1]^3*x[3]^3 + x[2]^3*x[3]^3
```
"""    
function primary_invariants(IR::InvRing)
   if !isdefined(IR, :primary)
      primary_invariants_via_singular(IR)
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

#######################################################

@doc Markdown.doc"""
    secondary_invariants(IR::InvRing)

Return a system of secondary invariants of `IR` with respect to the currently cached system of primary invariants of `IR`
(if no system of primary invariants of `IR` is cached, compute and cache such a system first).

If a system of secondary invariants is already cached, return the cached system. 
Otherwise, compute and cache such a system first.

NOTE: The secondary invariants are sorted by increasing degree.

# Examples
```jldoctest
julia> K, a = CyclotomicField(3, "a")
(Cyclotomic field of order 3, a)

julia> M1 = matrix(K, [0 0 1; 1 0 0; 0 1 0])
[0   0   1]
[1   0   0]
[0   1   0]

julia> M2 = matrix(K, [1 0 0; 0 a 0; 0 0 -a-1])
[1   0        0]
[0   a        0]
[0   0   -a - 1]

julia> G = MatrixGroup(3, K, [ M1, M2 ])
Matrix group of degree 3 over Cyclotomic field of order 3

julia> IR = invariant_ring(G)
Invariant ring of
Matrix group of degree 3 over Cyclotomic field of order 3
with generators
AbstractAlgebra.Generic.MatSpaceElem{nf_elem}[[0 0 1; 1 0 0; 0 1 0], [1 0 0; 0 a 0; 0 0 -a-1]]

julia> secondary_invariants(IR)
2-element Vector{AbstractAlgebra.Generic.MPoly{nf_elem}}:
 1
 x[1]^6*x[3]^3 + x[1]^3*x[2]^6 + x[2]^3*x[3]^6
```
"""    
function secondary_invariants(IR::InvRing)
   if !isdefined(IR, :secondary)
      secondary_invariants_via_singular(IR)
   end
   return copy(IR.secondary)
end

function irreducible_secondary_invariants(IR::InvRing)
   if !isdefined(IR, :irreducible_secondary)
      secondary_invariants_via_singular(IR)
   end
   return copy(IR.irreducible_secondary)
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
