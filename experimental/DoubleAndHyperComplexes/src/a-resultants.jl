export a_resultant, a_resultant_complex, discriminant_complex

# When computing resultants with the Weyman complex via 
# 
#   Δ = det Rπ_*(K* ⊗ 𝒪(-α))
# 
# on some toric variety `X` for the projection `π : X × Spec R → R` 
# there are a lot of objects which need to first be constructed 
# from the input (e.g. monomial support sets) and should then 
# probably be kept. For instance, one has the degree of freedom to
# chose an arbitrary twist `α`. This can be used to obtain simpler 
# complexes for which the determinant is easier to compute. Also, 
# one maybe wants to recycle context objects from toric geometry 
# throughout multiple computations of resultants, etc. 
#   To this end we provide yet another context object. Without it
# it would be very difficult to pass a twist `α` for the computation,
# because it needs to live in the very specific grading group of 
# the cox ring of the toric variety which is being used. With this 
# infrastructure, the user can first create the context object 
# from their input data and then derive all other information from 
# that. 
mutable struct ResultantCtx
  A::Vector{<:Union{Matrix, MatrixElem}} # the monomial support sets
  f::Vector{<:MPolyRingElem} # the tautological polynomials (which might already have specialized coefficients)
  R::Ring # the parameter ring (the coefficient ring of the `parent` of the `f`s)
  X::NormalToricVariety # the toric variety which comes out of the support sets
  F::Vector{<:MPolyRingElem} # the homogenizations of the `f` in the Cox ring
  inner_ctx::NewToricCtx # inner and outer context objects for the direct image
  outer_ctx::ToricCtxWithParams
  complexes::Dict{FinGenAbGroupElem, <:AbsHyperComplex} # results for different twists

  function ResultantCtx(A::Vector{<:Union{Matrix, MatrixElem}})
    @assert !is_empty(A) "resultants of empty support sets are not allowed"
    @assert length(A) == ncols(first(A)) + 1 "wrong number of support sets"
    @assert all(ncols(a) == ncols(first(A)) for a in A) "support sets need to have the same length"
    return new(A)
  end

  function ResultantCtx(f::Vector{T}) where {T<:MPolyRingElem}
    @assert !is_empty(f) "resultants of empty collections are not allowed"
    @assert all(parent(ff) === parent(first(f)) for ff in f) "all polynomials need to have the same parent"
    @assert ngens(parent(first(f))) + 1 == length(f) "exactly one more equation than number of variables must be provided"
    A = _support_sets(f)
    return new(A, f, coefficient_ring(first(f)))
  end
end

function ResultantCtx(f::MPolyRingElem)
  return ResultantCtx(pushfirst!([derivative(f, i) for i in 1:ngens(parent(f))], f))
end

support_sets(ctx::ResultantCtx) = ctx.A

function tautological_polynomials(ctx::ResultantCtx)
  if !isdefined(ctx, :f)
    error("tautological polynomials have not been defined")
  end
  return ctx.f
end

function homogeneous_tautological_polynomials(ctx::ResultantCtx)
  if !isdefined(ctx, :F)
    S_ext = extended_cox_ring(ctx)
    R = coefficient_ring(ctx)
    supps = support_sets(ctx)
    var_groups = Vector{elem_type(R)}[]

    if !isdefined(ctx, :f)
      # If these polynomials are not defined, then we need the `coefficient_ring` 
      # to be a polynomial ring with the coefficients as variables. 
      @assert R isa MPolyRing "incompatible ring"
      offset = 0
      for A in supps
        push!(var_groups, gens(R)[offset+1:offset+nrows(A)])
        offset += nrows(A)
      end
    else
      var_groups = [collect(coefficients(f)) for f in ctx.f]
    end
    
    X = toric_variety(ctx)

    # Coordinates of Cartier divisors in X of the system's Newton polytopes
    div_coords = [map(u -> -minimum( grad*primitive_generator(u) ), rays(X)) for grad in supps];
    U = primitive_generator.(rays(X))

    ctx.F = elem_type(S_ext)[
      sum(c[j]*prod(S_ext[i]^(dot(A[j,:], U[i]) + dc[i]) for i in 1:n_rays(X); 
                    init=one(S_ext)) 
          for j in 1:size(A, 1); init=zero(S_ext)) 
      for (A, c, dc) in zip(supps, var_groups, div_coords)];
  end
  return ctx.F
end

function toric_variety(ctx::ResultantCtx) 
  if !isdefined(ctx, :X)
    ctx.X = _get_toric_variety(support_sets(ctx))
  end
  return ctx.X::NormalToricVariety
end

function coefficient_ring(ctx::ResultantCtx)
  if !isdefined(ctx, :R)
    ctx.R = _parent_for_resultant(support_sets(ctx))
  end
  return ctx.R
end

function cox_ring(ctx::ResultantCtx)
  return cox_ring(toric_variety(ctx))
end

function extended_cox_ring(ctx::ResultantCtx)
  return graded_ring(outer_ctx(ctx))
end

function inner_ctx(ctx::ResultantCtx)
  if !isdefined(ctx, :inner_ctx)
    ctx.inner_ctx = NewToricCtx(toric_variety(ctx))
  end
  return ctx.inner_ctx
end

function outer_ctx(ctx::ResultantCtx)
  if !isdefined(ctx, :outer_ctx)
    ctx.outer_ctx = _get_outer_toric_ctx(inner_ctx(ctx), coefficient_ring(ctx))
  end
  return ctx.outer_ctx
end

grading_group(ctx::ResultantCtx) = grading_group(cox_ring(ctx))


@doc raw"""
    a_resultant_complex(support_sets::Vector{T};
        ctx::ResultantCtx, twist
      ) where {T <: Union{Matrix, MatrixElem}}

Given a collection of ``n+1`` support sets ``Aᵢ ∈ ℤʳˣⁿ``, i.e. a matrix whose ``r = r(i)`` rows are exponent vectors 
of Laurent polynomials, compute a complex ``C`` of free modules over the ring ``R = ℤ[aᵢᵥ : i = 1,…,n+1, ν ∈ 1,…,r(i)]`` 
of coefficients ``aᵢᵥ`` of those polynomials such that ``\det(C)`` is a (possibly non-reduced) equation for the 
resultant ``Δ``. See [GZ26](@cite) for more details.
"""
function a_resultant_complex(support_sets::Vector{T}; 
    ctx::ResultantCtx=ResultantCtx(support_sets),
    twist::FinGenAbGroupElem=zero(grading_group(ctx))
  ) where {T <: Union{Matrix, MatrixElem}}

  n = ncols(first(support_sets))
  @assert all(ncols(A) == n for A in support_sets) "the number of columns (variables) must coincide"
  @assert length(support_sets) == n+1 "wrong number of support sets"

  S_ext = extended_cox_ring(ctx)
  K = HomogKoszulComplex(S_ext, homogeneous_tautological_polynomials(ctx))
  t = ZeroDimensionalComplex(graded_free_module(S_ext, [-twist]))
  Kt = tensor_product(K, t)
  return DirectImageComplex(outer_ctx(ctx), Kt)
end

### helper functions for filling the `ResultantCtx`
function _parent_for_resultant(support_sets::Vector{T}) where {T}
  # create a polynomial ring for the generic coefficients of the support matrices
  R, a = polynomial_ring(ZZ, vcat([[Symbol("a_$(k)_$(i)") for i in 1:nrows(A)] for (k, A) in enumerate(support_sets)]...))
  return R
end

function _get_toric_variety(support_sets::Vector{T}) where {T}
  # Newton polytopes of support sets
  Q = map(convex_hull, support_sets);
  # Toric variety specified by the normal fan of the Minkowski sum of all Newton polytopes
  return normal_toric_variety(reduce(minkowski_sum, Q));
end

function _get_toric_ctx(support_sets::Vector{T}, X::NormalToricVariety) where {T}
  return 
end

function _get_outer_toric_ctx(inner_ctx::NewToricCtx, R::MPolyRing)
  X = toric_variety(inner_ctx)
  S = cox_ring(X)
  kk = coefficient_ring(R)
  coef_map = MapFromFunc(QQ, R, x->R(kk(x)))
  S_ext, phi = change_base_ring(coef_map, S)
  return ToricCtxWithParams(inner_ctx, phi)
end

function _get_outer_toric_ctx(inner_ctx::NewToricCtx, R::Ring)
  X = toric_variety(inner_ctx)
  S = cox_ring(X)
  S_ext, phi = change_base_ring(R, S)
  return ToricCtxWithParams(inner_ctx, phi)
end

@doc raw"""
    a_resultant(support_sets::Vector{T};
        ctx::ResultantCtx, twist
      ) where {T <: Union{Matrix, MatrixElem}}

Compute the ``A``-resultant using the Weyman complex as outlined 
in [GKZ08](@cite). The input `support_sets` 
is a list of integer matrices, whose rows stand for the exponent 
vectors of the monomial support sets. The output is a polynomial
in the coefficients `a_i_j` associated to the `j`-th row of the 
`i`-th support set. 

One can use different twists to compute the resultant and the 
determinant of a specific complex. This twist can be specified 
with the respective keyword argument. The `ResultantCtx` is a 
context object which maintains data required for the computation 
of the resultant. In particular it can be used to obtain the 
correct parent for the twist and other information on intermediate 
data.

# Examples
```jldoctest
julia> supps = [[0 0; 0 1; 1 0], [0 0; 2 0; 1 1; 0 2], [0 0; 1 0; 0 1]]
3-element Vector{Matrix{Int64}}:
 [0 0; 0 1; 1 0]
 [0 0; 2 0; 1 1; 0 2]
 [0 0; 1 0; 0 1]

julia> delta = a_resultant(supps)
a_1_1^2*a_2_2*a_3_3^2 - a_1_1^2*a_2_3*a_3_2*a_3_3 + a_1_1^2*a_2_4*a_3_2^2 - 2*a_1_1*a_1_2*a_2_2*a_3_1*a_3_3 + a_1_1*a_1_2*a_2_3*a_3_1*a_3_2 + a_1_1*a_1_3*a_2_3*a_3_1*a_3_3 - 2*a_1_1*a_1_3*a_2_4*a_3_1*a_3_2 + a_1_2^2*a_2_1*a_3_2^2 + a_1_2^2*a_2_2*a_3_1^2 - 2*a_1_2*a_1_3*a_2_1*a_3_2*a_3_3 - a_1_2*a_1_3*a_2_3*a_3_1^2 + a_1_3^2*a_2_1*a_3_3^2 + a_1_3^2*a_2_4*a_3_1^2

a_1_1^2*a_2_2*a_3_3^2 - a_1_1^2*a_2_3*a_3_2*a_3_3 + a_1_1^2*a_2_4*a_3_2^2 - 2*a_1_1*a_1_2*a_2_2*a_3_1*a_3_3 + a_1_1*a_1_2*a_2_3*a_3_1*a_3_2 + a_1_1*a_1_3*a_2_3*a_3_1*a_3_3 - 2*a_1_1*a_1_3*a_2_4*a_3_1*a_3_2 + a_1_2^2*a_2_1*a_3_2^2 + a_1_2^2*a_2_2*a_3_1^2 - 2*a_1_2*a_1_3*a_2_1*a_3_2*a_3_3 - a_1_2*a_1_3*a_2_3*a_3_1^2 + a_1_3^2*a_2_1*a_3_3^2 + a_1_3^2*a_2_4*a_3_1^2

julia> ctx = Oscar.ResultantCtx(supps);

julia> G = grading_group(ctx)
Z

julia> is_associated(delta, a_resultant(supps; ctx, twist=G[1]))
true


```
"""
function a_resultant(support_sets::Vector{T};
    ctx::ResultantCtx=ResultantCtx(support_sets),
    twist::FinGenAbGroupElem=zero(grading_group(ctx))
  ) where {T <: Union{Matrix, MatrixElem}}
  return det(a_resultant_complex(support_sets; ctx, twist); upper_bound=length(support_sets))
end

@doc raw"""
    a_resultant_complex(F::Vector{T};
        ctx::ResultantCtx, twist
      ) where {T<:MPolyRingElem}

Given a system of ``n+1`` polynomials ``F = (f₀,…,fₙ)`` in ``n`` variables over a ring ``R``, 
compute a complex of ``R``-modules ``C``, such that ``\det(C) = 0`` is the resultant, 
i.e. the locus ``Δ ⊂ \mathrm{Spec} R`` over which a solution to the system ``F = 0`` exists. 
See [GKZ08](@cite) for more details. 

# Examples
```jldoctest
julia> supps = [[0 0; 0 1; 1 0], [0 0; 2 0; 1 1; 0 2], [0 0; 1 0; 0 1]]
3-element Vector{Matrix{Int64}}:
 [0 0; 0 1; 1 0]
 [0 0; 2 0; 1 1; 0 2]
 [0 0; 1 0; 0 1]

julia> cplx = a_resultant_complex(supps);

julia> matrix(map(cplx, 1))
[                                        -a_1_2*a_2_1*a_3_2 + a_1_3*a_2_1*a_3_3   -a_3_1   -a_1_1]
[                                         a_1_1*a_2_4*a_3_2 - a_1_3*a_2_4*a_3_1   -a_3_3   -a_1_2]
[-a_1_1*a_2_2*a_3_3 + a_1_1*a_2_3*a_3_2 + a_1_2*a_2_2*a_3_1 - a_1_3*a_2_3*a_3_1   -a_3_2   -a_1_3]

```
"""
function a_resultant_complex(F::Vector{T};
    ctx::ResultantCtx=ResultantCtx(F),
    twist::FinGenAbGroupElem=zero(grading_group(ctx))
  ) where {T<:MPolyRingElem}
  S_ext = extended_cox_ring(ctx)
  K = HomogKoszulComplex(S_ext, homogeneous_tautological_polynomials(ctx))
  t = ZeroDimensionalComplex(graded_free_module(S_ext, [-twist]))
  Kt = tensor_product(K, t)
  return DirectImageComplex(outer_toric_ctx_object, Kt)
end

@doc raw"""
    a_resultant(F::Vector{T};
        ctx::ResultantCtx, twist
      ) where {T<:MPolyRingElem}

Given a system of ``n+1`` polynomials ``F = (f₀,…,fₙ)`` in ``n`` variables over a ring ``R``, 
compute the resultant for the polynomial system as the determinant of the corresponding 
`a_resultant_complex`. See [GKZ08](@cite) for more details. 

# Examples
```jldoctest
julia> supps = [[0 0; 0 1; 1 0], [0 0; 2 0; 1 1; 0 2], [0 0; 1 0; 0 1]]
3-element Vector{Matrix{Int64}}:
 [0 0; 0 1; 1 0]
 [0 0; 2 0; 1 1; 0 2]
 [0 0; 1 0; 0 1]

julia> delta = a_resultant(supps)
a_1_1^2*a_2_2*a_3_3^2 - a_1_1^2*a_2_3*a_3_2*a_3_3 + a_1_1^2*a_2_4*a_3_2^2 - 2*a_1_1*a_1_2*a_2_2*a_3_1*a_3_3 + a_1_1*a_1_2*a_2_3*a_3_1*a_3_2 + a_1_1*a_1_3*a_2_3*a_3_1*a_3_3 - 2*a_1_1*a_1_3*a_2_4*a_3_1*a_3_2 + a_1_2^2*a_2_1*a_3_2^2 + a_1_2^2*a_2_2*a_3_1^2 - 2*a_1_2*a_1_3*a_2_1*a_3_2*a_3_3 - a_1_2*a_1_3*a_2_3*a_3_1^2 + a_1_3^2*a_2_1*a_3_3^2 + a_1_3^2*a_2_4*a_3_1^2

```
"""
function a_resultant(F::Vector{T};
    ctx::ResultantCtx=ResultantCtx(F),
    twist::FinGenAbGroupElem=zero(grading_group(ctx))
  ) where {T<:MPolyRingElem}
  return det(a_resultant_complex(F; ctx, twist); upper_bound=dim(toric_variety(ctx)))
end

function _support_sets(F::Vector{T}) where {T<:MPolyRingElem}
  pre = [transpose(reduce(hcat, AbstractAlgebra.exponent_vectors(f))) for f in F]
  return [matrix_space(ZZ, size(A)...)(A) for A in pre]
end

function _tautological_polynomial(A::Union{Matrix, MatrixElem})
  P, x = polynomial_ring(ZZ, ncols(A))
  R, a = polynomial_ring(ZZ, [Symbol("a_$i") for i in 1:nrows(A)])
  PR, transf = change_base_ring(R, P)
  f = sum(a*prod(v^A[i, k] for (k, v) in enumerate(gens(PR)); init=one(PR)) for (i, a) in enumerate(a); init=zero(PR))
  return f
end

### Discriminants of a single polynomial `f`
# In this case the additional equations are given by the partial derivatives 
# and we have the corresponding specialization to the coefficients of `f`. 
@doc raw"""
    discriminant_complex(A::Union{Matrix, MatrixElem};
        tautological_polynomial, ctx::ResultantCtx, twist::FinGenAbGroupElem
      )

Compute the discriminant for the monomial support set, whose exponent vectors 
are given by the rows of `A`; cf. [GKZ08](@cite).
"""
function discriminant_complex(A::Union{Matrix, MatrixElem};
    tautological_polynomial::MPolyRingElem=_tautological_polynomial(A), 
    ctx::ResultantCtx=ResultantCtx(f), 
    twist::FinGenAbGroupElem=zero(grading_group(ctx))
  )
  return discriminant_complex(tautological_polynomial; ctx, twist)
end

function discriminant(A::Union{Matrix, MatrixElem};
    tautological_polynomial::MPolyRingElem=_tautological_polynomial(A), 
    ctx::ResultantCtx=ResultantCtx(pushfirst!([derivative(tautological_polynomial, i) for i in 1:ncols(A)], tautological_polynomial)),
    twist::FinGenAbGroupElem=zero(grading_group(ctx))
  )
  return det(discriminant_complex(tautological_polynomial; ctx, twist); upper_bound=ncols(A))
end

@doc raw"""
    discriminant_complex(f::MPolyRingElem; 
        ctx::ResultantCtx, twist::FinGenAbGroupElem
      )

Compute a complex of modules over the `coefficient_ring` of `f`, whose determinant 
is the discriminant of `f` as in [GKZ08](@cite). Note that the coefficients of `f` 
need not be in a field. 

# Examples
```jldoctest
julia> supp = [0; 1; 2;;]
3×1 Matrix{Int64}:
 0
 1
 2

julia> f = Oscar._tautological_polynomial(supp)
a_3*x1^2 + a_2*x1 + a_1

julia> cplx = discriminant_complex(f);

julia> matrix(map(cplx, 1))
[         -a_2*a_3   2*a_3]
[2*a_1*a_3 - a_2^2     a_2]

julia> det(ans)
-4*a_1*a_3^2 + a_2^2*a_3

julia> ctx = Oscar.ResultantCtx(f); # try with another twist

julia> G = grading_group(ctx)
Z

julia> cplx = discriminant_complex(f; ctx, twist=2*G[1]);

julia> matrix(map(cplx, 1))
[-a_1     -a_2     -a_3]
[-a_2   -2*a_3        0]
[   0     -a_2   -2*a_3]

julia> det(ans)
-4*a_1*a_3^2 + a_2^2*a_3


```
"""
function discriminant_complex(f::MPolyRingElem;
    ctx::ResultantCtx=ResultantCtx(f),
    twist::FinGenAbGroupElem=zero(grading_group(ctx))
  )
  return a_resultant_complex(support_sets(ctx); ctx, twist)
end

@doc raw"""
    discriminant(f::MPolyRingElem; 
        ctx::ResultantCtx, twist::FinGenAbGroupElem
      )

Compute the discriminant of `f` as in [GKZ08](@cite).
Note that the coefficients of `f` need not be in a field. 

# Examples
```jldoctest
julia> supp = [0; 1; 2;;]
3×1 Matrix{Int64}:
 0
 1
 2

julia> f = Oscar._tautological_polynomial(supp)
a_3*x1^2 + a_2*x1 + a_1

julia> discriminant(f)
-4*a_1*a_3^2 + a_2^2*a_3
```
"""
function discriminant(f::MPolyRingElem; 
    ctx::ResultantCtx=ResultantCtx(f), 
    twist::FinGenAbGroupElem=zero(grading_group(ctx))
  )
  return det(discriminant_complex(f; ctx, twist); upper_bound=ngens(parent(f)))
end

