# This file is an old version of the algorithm that can compute (not all cases) of 
# basis_lie_highest_weight and is used only in runtests.jl to check that the newer algorithm matches
# There is code doubling, but I am not sure how the src part is going to change when its integrated with the other
# lie algebra work and would prefer to change this file after doing the integration to make sure that everything stays
# correct.

module MBOld

using Oscar
using ..Oscar: GAPWrap

struct SparseVectorSpaceBasis
  A::Vector{SRow{ZZRingElem}}
  pivot::Vector{Int}
end

"""
divides vector by gcd of nonzero entries, returns vector and first nonzero index
used: addAndReduce!
"""
function normalize(v::SRow{ZZRingElem})
  if isempty(v)
    return (0, 0)
  end

  pivot = first(v)[1]

  return divexact(v, gcd(map(y -> y[2], union(v)))), pivot
end

reduceCol(a, b, i::Int) = b[i] * a - a[i] * b

"""
for each pivot of sp.A we make entry of v zero and return the result
0 => linear dependent
* => linear independent, new column element of sp.A since it increases basis
invariants: the row of a pivot element in any column in A is 0 (except the pivot element)
           elements of A are integers, gcd of each column is 1
"""
function addAndReduce!(sp::SparseVectorSpaceBasis, v::SRow{ZZRingElem})
  A = sp.A
  pivot = sp.pivot
  v, newPivot = normalize(v)
  if newPivot == 0
    return 0
  end

  for j in 1:length(A)
    i = pivot[j]
    if i != newPivot
      continue
    end
    v = reduceCol(v, A[j], i)
    v, newPivot = normalize(v)
    if newPivot == 0
      return 0
    end
  end

  pos = findfirst(pivot .> newPivot)
  if (pos === nothing)
    pos = length(pivot) + 1
  end

  insert!(A, pos, v)
  insert!(pivot, pos, newPivot)

  return v
end

#### Lie algebras

function lieAlgebra(t::String, n::Int)
  L = codomain(Oscar.iso_oscar_gap(Oscar.lie_algebra(QQ, Symbol(t), n)))
  return L, NTuple{3,Vector{GAP.Obj}}(GAP.Globals.ChevalleyBasis(L))
end

"""
used to create tensorMatricesForOperators
"""
function matricesForOperators(L, hw, ops)
  M = GAP.Globals.HighestWeightModule(L, GAP.Obj(hw))
  mats = map(
    o ->
      sparse_matrix(transpose(matrix(QQ, GAP.Globals.MatrixOfAction(GAPWrap.Basis(M), o)))),
    ops,
  )
  denominators = map(y -> denominator(y[2]), union(union(mats...)...))
  d = lcm(denominators)# // 1
  mats = (A -> change_base_ring(ZZ, d * A)).(mats)
  return mats
end

function weightsForOperators(L, cartan, ops)
  asVec(v) = Vector{Int}(GAPWrap.Coefficients(GAPWrap.Basis(L), v))
  if any(iszero.(asVec.(ops)))
    error("ops should be non-zero")
  end
  nzi(v) = findfirst(asVec(v) .!= 0)
  return [[asVec(h * v)[nzi(v)] / asVec(v)[nzi(v)] for h in cartan] for v in ops]
end

#### tensor model

function kron(A, B)
  res = sparse_matrix(ZZ, nrows(A) * nrows(B), ncols(A) * ncols(B))
  for i in 1:nrows(B)
    for j in 1:nrows(A)
      new_row_tuples = Vector{Tuple{Int,ZZRingElem}}([(1, ZZ(0))])
      for (index_A, element_A) in union(getindex(A, j))
        for (index_B, element_B) in union(getindex(B, i))
          push!(new_row_tuples, ((index_A - 1) * ncols(B) + index_B, element_A * element_B))
        end
      end
      new_row = sparse_row(ZZ, new_row_tuples)
      setindex!(res, new_row, (j - 1) * nrows(B) + i)
    end
  end
  return res
end

# temporary fix sparse in Oscar does not work
function tensorProduct(A, B)
  temp_mat = kron(A, spid(nrows(B))) + kron(spid(nrows(A)), B)
  res = sparse_matrix(ZZ, nrows(A) * nrows(B), ncols(A) * ncols(B))
  for i in 1:nrows(temp_mat)
    setindex!(res, getindex(temp_mat, i), i)
  end
  return res
end

spid(n) = identity_matrix(SMat, ZZ, n)
tensorProducts(As, Bs) = (AB -> tensorProduct(AB[1], AB[2])).(zip(As, Bs))
tensorPower(A, n) = (n == 1) ? A : tensorProduct(tensorPower(A, n - 1), A)
tensorPowers(As, n) = (A -> tensorPower(A, n)).(As)

"""
Calculates the matrices g_i corresponding to the operator ops[i].
"""
function tensorMatricesForOperators(L, hw, ops)
  mats = []

  for i in 1:length(hw)
    if hw[i] <= 0
      continue
    end
    wi = Int.(1:length(hw) .== i) # i-th fundamental weight
    _mats = matricesForOperators(L, wi, ops)
    _mats = tensorPowers(_mats, hw[i])
    if size(mats)[1] > 0
      A = _mats[1]
      B = mats[1]
    end
    mats = mats == [] ? _mats : tensorProducts(mats, _mats)
  end
  return mats
end

"""
    basisLieHighestWeight(t::String, n::Int, hw::Vector{Int}; parallel::Bool = true) :: Tuple{Vector{Vector{UInt8}},Vector{SRow{ZZRingElem}}}

Compute a monomial basis for the highest weight module with highest weight ``hw`` (in terms of the fundamental weights), for a simple Lie algebra of type ``t`` and rank ``n``.

Example
======
```julia
julia> dim, monomials, vectors = PolyBases.MB.basisLieHighestWeight(:A, 2, [1,0])
(3, Vector{UInt8}[[0x00, 0x00, 0x00], [0x01, 0x00, 0x00], [0x00, 0x00, 0x01]], SparseArrays.SparseVector{Int64, Int64}[  [1]  =  1,   [2]  =  1,   [3]  =  1])
```
"""
function basisLieHighestWeight(t::String, n::Int, hw::Vector{Int}; roots=[]) #--- :: Tuple{Int64,Vector{Vector{UInt8}},Vector{SRow{ZZRingElem}}}
  L, CH = lieAlgebra(t, n)
  ops = CH[1] # positive root vectors
  # .. reorder..
  wts = weightsForOperators(L, CH[3], ops)
  wts = (v -> Int.(v)).(wts)

  mats = tensorMatricesForOperators(L, hw, ops)
  hwv = sparse_row(ZZ, [(1, 1)])

  monomials = compute(hwv, mats, wts)
  ZZx, x = polynomial_ring(ZZ, length(monomials[1]))
  monomials = [
    finish(push_term!(MPolyBuildCtx(ZZx), ZZ(1), Int.(mon))) for mon in monomials
  ]
  monomials = Set(monomials)
  return monomials
end

nullMon(m) = zeros(UInt8, m)

function compute(v0, mats, wts::Vector{Vector{Int}})
  m = length(mats)
  monomials = [nullMon(m)]
  lastPos = 0
  id(mon) = sum((1 << (sum(mon[1:i]) + i - 1) for i in 1:(m - 1)); init=1)
  e = [UInt8.(1:m .== i) for i in 1:m]
  maxid(deg) = id(deg .* e[1])

  blacklists = [falses(maxid(0))]
  blacklists[end][id(monomials[1])] = 1
  newMons(deg) = (push!(blacklists, falses(maxid(deg))))
  checkMon(mon) = blacklists[end - 1][id(mon)]
  setMon(mon) = (blacklists[end][id(mon)] = true)

  vectors = [v0]
  weights = [0 * wts[1]]
  space = Dict(weights[1] => SparseVectorSpaceBasis([], []))
  addAndReduce!(space[weights[1]], v0)

  deg = 0
  while length(monomials) > lastPos
    startPos = lastPos + 1
    newPos = length(monomials)
    deg = deg + 1
    newMons(deg)
    # iteration in degrevlex ordering 
    for i in m:-1:1, di in deg:-1:1
      for p in startPos:newPos
        if !all(monomials[p][(i + 1):m] .== 0)
          continue
        end

        if monomials[p][i] + 1 > di
          startPos = p + 1
          continue
        end
        if monomials[p][i] + 1 < di
          break
        end

        mon = monomials[p] + e[i]

        if any(i != j && mon[j] > 0 && !checkMon(mon - e[j]) for j in 1:m)
          continue
        end

        wt = weights[p] + wts[i]
        if !haskey(space, wt)
          space[wt] = SparseVectorSpaceBasis([], [])
        end

        vec = vectors[p] * transpose(mats[i])
        vec = addAndReduce!(space[wt], vec)
        if vec == 0
          continue
        end

        setMon(mon)
        push!(monomials, mon)
        push!(weights, wt)
        push!(vectors, vec)
      end
    end
    lastPos = newPos
  end

  return monomials
end

end
