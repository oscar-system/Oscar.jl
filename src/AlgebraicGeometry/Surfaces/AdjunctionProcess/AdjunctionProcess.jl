
##############################################################################
#
# adjunction process for surfaces
#
##############################################################################

#function _random_vector(m::Int, I::MPolyIdeal, bound::Int = 30000)
#   @assert m>0
#   R = base_ring(I)  
#   SM = Singular.LibRandom.randommat(m, 1, Oscar.singular_generators(I), bound)
#   return [R(SM[i,1]) for i=1:nrows(SM)]
#end

#function _random_matrix(m::Int, n::Int, I::MPolyIdeal, bound::Int = 30000)
#   @assert n>0
#   @assert m>0
#   R = base_ring(I)  
#   SM = Singular.LibRandom.randommat(m, n, Oscar.singular_generators(I), bound)
#   return Matrix(matrix(R, nrows(SM), ncols(SM), [R(SM[i,j]) for i=1:nrows(SM) for j=1:ncols(SM)]))
#end

function _random_vector(R::MPolyDecRing, m::Int, d::Int)
   @assert m>0
   rd =  _random_matrix(R, m, 1, d)
   return [R(rd[i,1]) for i=1:nrows(rd)]
end

function _random_matrix(R::MPolyDecRing, m::Int, n::Int, d::Int)
   @assert n>0
   @assert m>0
   F = graded_free_module(R, m)
   G = graded_free_module(R, n)
   G = twist(G, d)
   H, h = hom(F, G)
   phi = rand_homogeneous(H, 0)
   return matrix(h(phi))
end

function _arithmetic_genus(I::MPolyIdeal)
   S = base_ring(I)
    #req is_standard_graded(S) "The base ring must be standard graded."
   A, _ = quo(S, I)
   H = hilbert_polynomial(A)
   return (-1)^(dim(A)-1)*(ZZ(coeff(H, 0)) - 1)
end

@doc raw"""
    sectional_genus(X::AbsProjectiveVariety)

Given a subvariety `X` of some $\mathbb P^n$, return the arithmetic genus of the intersection of `X`
with a general linear subspace of $\mathbb P^n$ of dimension $c+1$.

# Examples
```jldoctest
julia> X = bordiga()
Projective variety
  in projective 4-space over GF(31991) with coordinates [x, y, z, u, v]
defined by ideal with 4 generators

julia> dim(X)
2

julia> codim(X)
2

julia> degree(X)
6

julia> sectional_genus(X)
3
```
"""
function sectional_genus(X::AbsProjectiveVariety)
   I = defining_ideal(X)
   S = base_ring(I)
   nl = length(gens(S))-codim(I)-2
   rd = _random_vector(S, nl, 1)
   J = ideal(S, rd)+I
   dd = dim(J)
   while !(dd == 2)
      rd = _random_vector(S, nl, 1)
      J = ideal(S, rd)+I
      dd = dim(J)
   end
   return _arithmetic_genus(J)
end

@doc raw"""
    is_linearly_normal(X::AbsProjectiveVariety)

Return `true` if `X` is linearly normal, and `false` otherwise.

# Examples
```jldoctest
julia> X = bordiga()
Projective variety
  in projective 4-space over GF(31991) with coordinates [x, y, z, u, v]
defined by ideal with 4 generators

julia> dim(X)
2

julia> codim(X)
2

julia> is_linearly_normal(X)
true
```
"""
function is_linearly_normal(X::AbsProjectiveVariety)
   I = defining_ideal(X)
   M = ideal_as_module(I)
   n = length(gens(base_ring(I)))
   tbl = sheaf_cohomology(M, -2, 1)
   return tbl[1, 1] == 0
end

function _adjoint_matrix(D::AbstractAlgebra.Generic.MatSpaceElem)
     #@req !any(x->is_zero(x) ? true : degree(Int, x) > 1, Oscar._vec(D)) "The code assumes that the input matrix has linear entries."
     R = base_ring(D)
     K = coefficient_ring(R)
     r = ncols(D)
     P, _ = graded_polynomial_ring(K, :z => 1:r; cached=false)
     RP, _ = graded_polynomial_ring(K, vcat(symbols(R), symbols(P)); cached=false)
     embRRP = hom(R, RP, gens(RP)[1:ngens(R)])    
     projRPP = hom(RP, P, vcat([zero(P) for i = 1:ngens(R)], gens(P)))
     M = map(embRRP, D)*matrix(gens(RP)[(ngens(R)+1):ngens(RP)])
     DM = [map_entries(x->derivative(x,i), M) for i = 1:ngens(R)]
     MM = transpose(hcat(DM...));
     return transpose(map(projRPP, MM))
end

@doc raw"""
    canonical_bundle(X::AbsProjectiveVariety)

Given a smooth projective variety `X`, return a module whose sheafification is the canonical bundle of `X`.

!!! note
    The function does not check smoothness. If you are uncertain, enter `is_smooth(X)` first.

# Examples

```jldoctest
julia> R, x = graded_polynomial_ring(QQ, :x => (1:6));

julia> I = ideal(R, [x[1]*x[6] - x[2]*x[5] + x[3]*x[4]]);

julia> GRASSMANNIAN = variety(I);

julia> Omega = canonical_bundle(GRASSMANNIAN)
Graded submodule of R^1 with 1 generator
  1: e[1]
represented as subquotient with no relations

julia> degrees_of_generators(Omega)
1-element Vector{FinGenAbGroupElem}:
 [4]
```

```
julia> R, (x, y, z) = graded_polynomial_ring(QQ,[:x, :y, :z]);

julia> I = ideal(R, [y^2*z + x*y*z - x^3 - x*z^2 - z^3]);

julia> ELLCurve = variety(I);

julia> Omega = canonical_bundle(ELLCurve)
Graded submodule of R^1 with 1 generator
  1: e[1]
represented as subquotient with no relations

julia> degrees_of_generators(Omega)
1-element Vector{FinGenAbGroupElem}:
 [0]
```

```jldoctest
julia> X = bordiga()
Projective variety
  in projective 4-space over GF(31991) with coordinates [x, y, z, u, v]
defined by ideal with 4 generators

julia> dim(X)
2

julia> codim(X)
2

julia> Omega = canonical_bundle(X);

julia> typeof(Omega)
SubquoModule{MPolyDecRingElem{fpFieldElem, fpMPolyRingElem}}
```
"""
function canonical_bundle(X::AbsProjectiveVariety)
  Pn = ambient_coordinate_ring(X)
  A = homogeneous_coordinate_ring(X)
  n = ngens(Pn)-1
  c = codim(X)
  FA = free_resolution(A, algorithm = :fres)
  C_simp = simplify(FA)
  C_shift = shift(C_simp, c)
  OmegaPn = graded_free_module(Pn, [n+1])
  D = hom(C_shift, OmegaPn)
  D_simp = simplify(D)
  Z, inc = kernel(D_simp, 0)
  B, inc_B = boundary(D_simp, 0)
  return prune_with_map(SubquoModule(D_simp[0], ambient_representatives_generators(Z), ambient_representatives_generators(B)))[1] 
end

@doc raw"""
    adjunction_process(X::AbsProjectiveVariety, steps::Int=0)

Given a smooth surface `X` and a non-negative integer `steps`, return data which describes the adjunction process for `X`: 
If `steps == 0`, carry out the complete process. Otherwise, carry out the indicated number of steps only.

More precisely, if $X^{(0)} = X \rightarrow X^{(1)}\rightarrow \dots \rightarrow X^{(r)}$ is the sequence of successive adjunction maps and 
adjoint surfaces in the completed adjunction process, return a quadruple `L`, say, where: 

`L[1]` is a vector of tuples of numerical data: For each step $X^{(i)}\rightarrow X^{(i+1)}$, return the tuple $(n^{(i)}, d^{(i)}, \pi^{(i)}, s^{(i)}),$ 
where $n^{(i)}$ is the dimension of the ambient projective space of $X^{(i)}$, $d^{(i)}$ is the degree of $X^{(i)}$, $\pi^{(i)}$ is the sectional genus of $X^{(i)}$, 
and $s^{(i)}$ is the number of exceptional $(-1)$-lines on $X^{(i)}$ which are blown down to points in $ X^{(i+1)}$.

`L[2]` is a vector of adjoint matrices: For each step $X^{(i)}\rightarrow X^{(i+1)}$, return a presentation matrix of $S_X^{(i)}(1)$ considered 
as a module over $S_X^{(i+1)}$, where the $S_X^{(i)}$ are the homogeneous coordinate rings of the $X^{(i)}$. If `X` is rational, these matrices
can be used to compute a rational parametrization of `X`.

`L[3]` is a vector of zero-dimensional projective algebraic sets: For each step $X^{(i)}\rightarrow X^{(i+1)}$, return the union of points in $ X^{(i+1)}$
which are obtained by blowing down the exceptional $(-1)$-lines on $X^{(i)}$.
    
`L[4]` is a projective variety: Return the last adjoint surface $X^{(r)}$. 


!!! note
    The function does not check whether `X` is smooth. If you are uncertain, enter `is_smooth(X)` first.

!!! warning
    At current state, the adjunction process is only implemented for rational and Enriques surfaces which are linearly normal in the given embedding. The function does not check whether `X` is rational or an Enriques surface. In fact, at current state, `OSCAR` does not offer direct checks for this. Note, however, that the adjunction process will give an answer to this question a posteriori in cases where it terminates with a surface which is known to be rational or an Enriques surface.

# Examples
```jldoctest
julia> X = bordiga()
Projective variety
  in projective 4-space over GF(31991) with coordinates [x, y, z, u, v]
defined by ideal with 4 generators


julia> dim(X)
2

julia> codim(X)
2

julia> L = adjunction_process(X);

julia> L[1]
2-element Vector{NTuple{4, ZZRingElem}}:
 (4, 6, 3, 0)
 (2, 1, 0, 10)

julia> L[4]
Projective variety
  in projective 2-space over GF(31991) with coordinates [z[1], z[2], z[3]]
defined by ideal (0)

julia> L[3][1]
Projective algebraic set
  in projective 2-space over GF(31991) with coordinates [z[1], z[2], z[3]]
defined by ideal with 5 generators

julia> dim(L[3][1])
0

julia> degree(L[3][1])
10
```

```
julia> X = rational_d9_pi7();

julia> L = adjunction_process(X);

julia> L[1]
3-element Vector{NTuple{4, ZZRingElem}}:
 (4, 9, 7, 0)
 (6, 9, 4, 6)
 (3, 3, 1, 3)
```

!!! note
    Inspecting the  returned numerical data in the first example above, we see that the Bordiga surface is the blow-up of the projective plane in 10 points, embedded into projective 4-space by the linear system $H = 4L -\sum_{i=1}^{10} E_i$. Here, $L$ is the preimage of a line and the $E_i$ are the exceptional divisors. In the second example, we see from the output that the terminal object of the adjunction process is a Del Pezzo surface in projective 3-space, that is, the blow-up of the projective plane in 6 points. In sum, we see that `X` is the blow-up of the projective plane in 15 points, embedded into projective 4-space by the linear system $H = 9L - \sum_{i=1}^{6} 3E_i - \sum_{i=7}^{9} 2E_i - \sum_{i=10}^{15} E_i$. 
"""
function adjunction_process(X::AbsProjectiveVariety, steps::Int = 0)
   @assert steps >= 0
   Pn = ambient_coordinate_ring(X)
   I = defining_ideal(X)
   Y = X
   @req dim(I) == 3 "The given variety is not a surface."
   @req is_linearly_normal(X) "The given variety is not linearly normal."
   # TODO If X is not linear normal, embed it as a linearly normal variety first.
   numlist = [(ZZ(ngens(Pn)-1), degree(X), sectional_genus(X), zero(ZZ))]
   ptslist = ProjectiveAlgebraicSet[]
   adjlist = AbstractAlgebra.Generic.MatSpaceElem[] 
   Omega = canonical_bundle(X)
   FOmega = free_resolution(Omega, length = 1, algorithm = :mres)
   D = matrix(map(FOmega,1))
   count = 1
   while nrows(D) > 2 && (steps == 0 || count <= steps)
      if !any(x -> total_degree(x) > 1, Oscar._vec(D)) 
         adj = _adjoint_matrix(D)
      else
        return (numlist, adjlist, ptslist, variety(I, check = false, is_radical = false))
      end
      rd = _random_matrix(Pn, 3, ncols(adj), 0)
      dd = dim(ideal(Pn, rd*gens(Pn))+I) 
      while !(dd == 0)
         rd = _random_matrix(Pn, 3, ncols(adj), 0)
         dd = dim(ideal(Pn, rd*gens(Pn))+I)
      end

      I = annihilator(cokernel(adj))
      Pn = base_ring(I)
    
      Ipts = ideal(Pn, [zero(Pn)])    
      RD =  map_entries(constant_coefficient, rd)   
      for i = 1:3
          KKi = map_entries(Pn, kernel(matrix(RD[i,:])))
          AAi = annihilator(cokernel(adj*transpose(KKi)))
	  Ipts = Ipts+AAi
      end
      Ipts = saturation(Ipts, ideal(Pn, gens(Pn)))
      pts = algebraic_set(Ipts,  is_radical = false, check = false) 
      if dim(pts) == 0
         l = degree(pts)
      else
         l = zero(ZZ)
      end
      Y = variety(I, check = false, is_radical = false)
      dY = degree(Y)
      piY = sectional_genus(Y)
      dummy = (ZZ(ngens(Pn)-1), dY, piY)
      if l==0 && dummy == numlist[count][1:3]   # Enriques surface
	  return (numlist, adjlist, ptslist, Y)
      else
        push!(numlist, (ZZ(ngens(Pn)-1), dY, piY, l))
        push!(adjlist, adj)
        push!(ptslist, pts)
        Omega = canonical_bundle(Y)
        FOmega = free_resolution(Omega, length = 1, algorithm = :mres)
        D = matrix(map(FOmega,1))
      end
      count = count+1
   end
   return (numlist, adjlist, ptslist, Y)
end



