# I would like to credit Dario Mathiä who produced an initial version of the following code

export wreath_macs, fake_deg

Oscar.canonical_unit(x::QQAbFieldElem{AbsSimpleNumFieldElem}) = x #meanwhile canonical_unit gets fixed

# Tools

#gives indices of a subarray
function subarray_indices(sub,arr)
  return [findfirst(==(s), arr) for s in sub]
end

# Computes the b-invariant
function b1(lambda::Vector{Int})
  return sum((i-1)*l_i for (i, l_i) in enumerate(lambda); init=0)
end

# Computes the b-invariant of a multipartition
function b_inv(lbb::Multipartition)
  return b1(map(lambda -> sum(lambda)[1], lbb)) + length(lbb)*sum(i -> b1(data(lbb[i])),1:length(lbb))
end

function beta_number_to_partition(beta::Vector{Int})
  lb=Int[]
  lbeta=length(beta)
  for j in 1:lbeta
    nbholes=beta[j]-beta[1]-(j-1)
    if nbholes >= 1
    pushfirst!(lb, nbholes)
    end
  end
  return partition(lb)
end

# Residue notation is in terms of rows and columns

function core(coroot::Vector{Int},r::Int)
  beta=Int[]
  m=minimum(coroot)
  M=maximum(coroot)
  for k in m:M
    for j in 1:r
       if k <= coroot[j]
          push!(beta, k*r+(j-1))
       end
    end
  end
  beta=sort!(beta)
  return beta_number_to_partition(beta)
end

#Inspired by Sagemath's algorithm to obtain a partition from core and quotient

function tau_om(lbb::Multipartition, wperm::PermGroupElem, coroot::Vector{Int})
   r=length(lbb)
   w0=perm([r-i for i in 0:(r-1)])
   wperm=inv(wperm)
   lbb_perm=multipartition(permuted(lbb.mp,wperm*w0))
   gamma=core([-x for x in coroot],r)
   lg=length(gamma)
   k=r*maximum(length(lbb_perm[i]) for i in 1:r) + lg
   v=[[gamma[i]-(i-1) for i in 1:lg]; [-i for i in lg:(k-1)]]
   w=[[x for x in v if mod((x-i),r) == 0] for i in 1:r]
   new_w=[]
   for i in 1:r
     lw=length(w[i])
     lq=length(lbb_perm[i])
     new_w=[new_w;[w[i][j] + r*lbb_perm[i][j] for j in 1:lq]]
     new_w=[new_w;[w[i][j] for j in lq+1:lw]]
   end
   sort!(new_w,rev=true)
   new_w=[new_w[i]+(i-1) for i in 1:length(new_w)]
   filter!(x-> x!=0,new_w)
   return partition(new_w)
end

#Coroots are given in the canonical basis (e_i) of Z^I.
#Note also that \alpha^{\vee}_i=e_{i-1}-e_i.
#The coroot lattice is thus the lattice of null sum elements in Z^I.

#Reoders the CharTable from Chevie
function reorder(charTirr::Matrix{QQAbFieldElem{AbsSimpleNumFieldElem}}, r::Int, n::Int, Modules::Vector{Multipartition{Int}}, K::QQAbField{AbsSimpleNumField})
  mps=multipartitions(n,r)
  new_ord=Int[]
    for lbb in mps
      i=findfirst(x-> x==lbb,Modules)
      push!(new_ord,i)
    end
  return charTirr[new_ord,:]
end

#Tools from representation theory

# NEED TO FIX THAT (issue is with the fraction_field and canonical_unit...)
# Computes the fake degree of a multipartition cf. Ste89 Thm 5.3 and Prop. 3.3.2 Haiman cdm
function fake_deg(lbb::Multipartition, K::QQAbField{AbsSimpleNumField}, Q::AbstractAlgebra.Generic.FracField{AbstractAlgebra.Generic.MPoly{QQAbFieldElem{AbsSimpleNumFieldElem}}}, var::AbstractAlgebra.Generic.MPoly{QQAbFieldElem{AbsSimpleNumFieldElem}})
  r=size(lbb)[1]
  n=sum(lbb)[1]
  res=Q(1)
  for p in 1:n
    res=res*(1-var^(r*p))
  end
  for k in 1:r
    for i in 1:length(lbb[k])
      for j in 1:lbb[k][i]
        hookfactor=(Q(1)-var^(r*(1+lbb[k][i]+length(findall(row->row >= j,lbb[k]))-j-i)))
        res=res//hookfactor
      end
    end
  end
  res=res*var^b_inv(lbb)
  return res
end

# computes C_Delta defined in PhD relation (1.45)
function C_Delta(r::Int, n::Int, K::QQAbField{AbsSimpleNumField}, Q::AbstractAlgebra.Generic.FracField{AbstractAlgebra.Generic.MPoly{QQAbFieldElem{AbsSimpleNumFieldElem}}}, var::AbstractAlgebra.Generic.MPoly{QQAbFieldElem{AbsSimpleNumFieldElem}})
  charTable=character_table_complex_reflection_group(r,1,n)
  Modules=class_parameters(charTable)
  charTirr=[charTable[i,j] for i in 1:nrows(charTable), j in 1:ncols(charTable)]
  Modules=[multipartition(map(x->partition(x),lbb)) for lbb in Modules]
  l=length(Modules)
  charTirr=reorder(charTirr,r,n,Modules,K)
  charTirrT=solve_init(matrix(K,transpose(charTirr)))
  mps=collect(multipartitions(n,r))
  rows=zero_matrix(Q,l,0)
  for i in 1:l
    col=zero_matrix(Q,l,1)
    for j in 1:l
      v=matrix(K,l,1,map(k->charTirr[i,k]*charTirr[j,k],1:l))
      x=solve(charTirrT,v,side=:right)
      col=col+fake_deg(mps[j],K,Q,var)*x
    end
    rows=hcat(rows,col)
  end
  return rows
end

# computes the character of the simple representations of the Cherednik algebra.
function C_L(r::Int, n::Int, K::QQAbField{AbsSimpleNumField}, Q::AbstractAlgebra.Generic.FracField{AbstractAlgebra.Generic.MPoly{QQAbFieldElem{AbsSimpleNumFieldElem}}}, var::AbstractAlgebra.Generic.MPoly{QQAbFieldElem{AbsSimpleNumFieldElem}})
  C_D=C_Delta(r,n,K,Q,var)
  mps=multipartitions(n,r)
  d = [var^b_inv(mp)//fake_deg(mp,K,Q,var) for mp in mps]
  diag=diagonal_matrix(Q,d)
  return diag * C_D
end

function bigger_ord(lbb::Multipartition, wperm::PermGroupElem, coroot::Vector{Int})
  r=length(lbb)
  n=sum(lbb)
  mps=multipartitions(n,r)
  lbbquot=tau_om(lbb,wperm,coroot)
  res=map(x->x.mp,Iterators.filter(x-> dominates(tau_om(x,wperm,coroot),lbbquot),mps))
  return res
end

function smaller_ord(lbb::Multipartition, wperm::PermGroupElem, coroot::Vector{Int})
  r=length(lbb)
  n=sum(lbb)
  mps=multipartitions(n,r)
  lbbquot=tau_om(lbb,wperm,coroot)
  res=map(x->x.mp,Iterators.filter(x-> dominates(lbbquot,tau_om(x,wperm,coroot)),mps))
  return res
end

function wreath_macs(n::Int, r::Int, wperm::PermGroupElem, coroot::Vector{Int})
  K,_ = abelian_closure(QQ)
  R, (q,t) = polynomial_ring(K, [:q,:t])
  Q = fraction_field(R)

  mps=multipartitions(n,r)
  l=length(mps)

  c_L=C_L(r,n,K,Q,t)
  c_L_q=map_entries(f->numerator(f)(0,q)//denominator(f)(0,q),c_L)
  c_L_tinv=map_entries(f->numerator(f)(0,1//t)//denominator(f)(0,1//t),c_L)

  rows=[]
  for lbb in mps
    smallers=smaller_ord(lbb, wperm, coroot) #tinv
    biggers=bigger_ord(lbb, wperm, coroot) #q
    smaller_indices=sort!(subarray_indices(smallers,mps))
    bigger_indices=sort!(subarray_indices(biggers,mps))
    sub_smaller_tinv=vcat(map(i->c_L_q[i:i,:],smaller_indices)...)
    sub_bigger_q=vcat(map(i->c_L_tinv[i:i,:],bigger_indices)...)
    M=vcat(sub_smaller_tinv,sub_bigger_q)
    B=kernel(M, side=:left)
    rows=vcat(rows,B[1,1:length(smaller_indices)]*sub_smaller_tinv)
  end
  c_L_qt = matrix(Q,l,l,rows)

  triv=[partition([n])]
  for i in 2:r
    push!(triv,partition([]))
  end
  triv=multipartition(triv)
  index_triv=findfirst(x-> x==triv,mps)
  d=[1//c_L_qt[i,index_triv] for i in 1:l]
  diag=diagonal_matrix(Q,d)
  c_L_qt_H=diag*c_L_qt
  return c_L_qt_H
end
