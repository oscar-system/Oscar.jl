module BlockSys

using Oscar

struct BlockSystems
  n::Int
  l::Int
  cur::Vector{Vector{Int}}
  function BlockSystems(n::Int, l::Int)
    @assert n % l == 0
    return new(n, l, [collect((i-1)*l+1:i*l) for i=1:divexact(n, l)])
  end
end

function Base.iterate(B::BlockSystems)
  return B.cur, deepcopy(B.cur)
end

function Base.iterate(B::BlockSystems, st::Array{Vector{Int}})
  if B.l==1||B.l==B.n
    return nothing
  end
  i = length(B.cur)-1
  while true
    j = B.l
    while true
      if st[i][j] < B.n - B.l + j
        st[i][j] += 1
        free = Set(1:B.n)
        for l=1:i-1
          setdiff!(free, st[l])
        end
        if !(st[i][j] in free) 
          continue
        end
        if length(intersect(free, Set(st[i][j]+1:B.n)))<B.l-j
          continue
        end
        setdiff!(free, st[i][1:j])
        while j < B.l
          j += 1
          I = intersect(free, Set(st[i][j-1]:B.n))
          if isempty(I)
            break
          end
          st[i][j] = minimum(I)
          pop!(free, st[i][j])
        end
        i += 1
        while i <= length(st)
          for j=1:B.l
            st[i][j] = minimum(free)
            pop!(free, st[i][j])
          end
          i += 1
        end
        return deepcopy(st), st
      end
      j -= 1
      if j == 1
        i -= 1
        i == 0 && return nothing
        break
      end
    end
  end
end
Base.IteratorSize(::BlockSystems) = Base.SizeUnknown()
<<<<<<< Updated upstream


=======
################################################################################################################################
# collect_subfields (gives a list and its blocksystems(debbuging) for valid subfield structures)
###############################################################################################################################
function collect_subfields(C::Oscar.GaloisGrp.GaloisCtx{Hecke.qAdicRootCtx},filterinvar=false)
  FieldList,BL = [],[]
  K, _ = number_field(C.f)
  S = Oscar.GaloisGrp.SubfieldLattice(K, C)
  n = degree(C.f)
  Rts = roots(C,7)
  cyc = gap_perm([findfirst(y -> y == x, Rts) for x = map(frobenius,Rts)])  # get the cyc from froeb.
  for d in divisors(n) # generate all possible blocksys.
    for bls in BlockSystems(n,d)
      if filterinvar   # Here seems to be a bug
        if !isinvar(bls, cyc)
          @show bls, cyc
        else
          push!(BL, bls)
        end
      else 
        push!(BL,bls) 
      end 
    end 
  end
##################################################################################################
  VB = [] # found valid blocksys
##################################################################################################
  for bls in BL
      SF = Nothing()
      !istriv(bls) || continue
      # pseudotest / other with VB
      #TODO
      SF = Oscar.GaloisGrp.subfield(S, bls)
      if typeof(SF) == Tuple{AnticNumberField, NfToNfMor} # i.e not equal nothing -> update SF
        push!(FieldList,(bls,SF))
        push!(VB,bls)
      end
      SF = Nothing()
  end 
  #debug 
  for (a,b) in FieldList 
    if b !=  Oscar.GaloisGrp.subfield(S, a) 
      @warn "Wrong matching >>> ($b, $a)"
    end 
  end 
  return FieldList
end
function isinvar(bls,cyc)
  #return true if no contradiction with cyc
  for bl in bls
    Set(cyc.(bl)) in Set.(bls) || return false 
  end 
  return true 
end
function istriv(bls)
  length(bls)!=1 || return true
  iszero(length.(bls).-1) ? (return true) : (return false)  
end

end 
################################################################################################################################
# End Module
################################################################################################################################


################################################################################################################################
# basic tree with unlimited children per node, inpush updated paths to all leaves
# ids are updated during puhses, nodes are stored by id in a Dict
# acces via node-id. 
# update leafs in the end #TODO improve this feature inpush
################################################################################################################################
mutable struct Node
  parent::Union{Int64,Missing}
  children::Vector{Int64}
  nchildren::Int64
  value::Any   #later d,n and Vector here
  path::Vector{Int64}  #<- TODO new type path with attributes length and nodeids etc...
end
mutable struct Tree  
  nodes::Dict{Int64, Node} # for performance reasons maybe Array here (unnecessary to call hash function on our practical size)
  nnodes::Int64
  leafs::Vector{Int64}
end 
mutable struct compressed_Blocksystem
  T::Tree # dataflowtree with k and used Cyc
  d::Int64 # blocksize
  Cyc::Vector{Vector{Int64}} # disjoint cycles
  N::Vector{Int64}
end 
function Tree(value)
  return Tree(Dict(zip([1],[Node(missing,Vector{Int64}(),0,value,[1])])),1,Vector{Int64}())
end 
function root(T::Tree)
  return T.nodes[1]
end
function pushchild!(tree::Tree, parentid::Int , value::Any)
  1 <= parentid <= tree.nnodes || throw(BoundsError(tree,parentid))
  tree.nnodes += 1
  push!(tree.nodes[parentid].children,tree.nnodes) # add child id to parents child array
  tree.nodes[parentid].nchildren+=1     # add child counter
  push!(tree.nodes,(tree.nnodes => Node(parentid,Vector{Int64}(),0,value,vcat(tree.nodes[parentid].path,tree.nnodes)))) # push new child to the Dictionary & Update path
  return tree
end
function isroot(n::Node)
  typeof(n.parent) == Missing || return false 
  return true 
end 
function isleaf(n::Node)
  n.children == Vector{Int64}() || return false
  return true 
end 
function leafs(T::Tree)
  leafs = []
  for n_i in 1:T.nnodes
    !isleaf(T.nodes[n_i]) || push!(leafs,n_i) 
  end 
  return leafs
end
function all_root2leaf_paths(T::Tree)
  P = []
  for n_i in 1:T.nnodes
      !isleaf(T.nodes[n_i]) || push!(P,T.nodes[n_i].path) 
  end 
  return P
end
################################################################################################################################

#Also good idea:
#
# find all pos. blocks cont 1: update in a BST (by Orbits of perm.) During iter we can avoid doublecheck and recalc orbits due to the tree struct.
#TODO 
# Runtime tests
# Pseudotests & new Date Structure 
# 
function block_klueners(C::Oscar.GaloisGrp.GaloisCtx{Hecke.qAdicRootCtx},d) # change signature and function name (after i fixed my new BS type)
  # if d = 1 manualle or else k= 1 = d and dk - n1 < 0 possible
  Rts = roots(C,7)
  cyc = gap_perm([findfirst(y -> y == x, Rts) for x = map(frobenius,Rts)])
  n = degree(C.f)
  divexact(n,d) # catch errors fast
  Cyc = map(collect,orbits(gset(sub(symmetric_group(n),[cyc])[1]))) 
  l = length(Cyc)
  Z_values = Dict(zip([i for i in 1:l],length.(Cyc)))
  T = Tree(0)
  O = kdivisors(vcat(length.(Cyc),[d]))
  #O = divisors(d)
  A = Set([i for i=1:l])
  T = rec_klueners(T,d,A,Z_values,cyc,1,O)
  T.leafs =  leafs(T)
  return compressed_Blocksystem(T,d,Cyc,length.(Cyc))
  # returns the blocksys structure Tree 
end
function kdivisors(U)
  O = Set()
  for g in U # experimental here
    O = union(O,Set(divisors(g)))
  end 
  return O
end 
function rec_klueners(T,d,A,Z_values,cyc,parentid,g) 
  l = length(Z_values)#
  for k in g   # note there exists k_i for each conjugated block 
    for set in specialsets(A,Z_values,k*d,k) # find all subsets. st. dk-minA = sum others
      pushchild!(T,parentid,(k => set))
      A_new = setdiff(A,set)
      isempty(A_new) ||  rec_klueners(T,d,A_new,Z_values,cyc,T.nnodes,g) # append to last pushed child (with id T.nnodes)
    end 
  end
  return T 
end 
function specialsets(A,Zv,v,k)
  #return subsets A s.t. kd - min(A) = Sum n_a    , k | n_a inkl min.
  if isempty(A)
    return []
  end 
  m = minimum(A)
  M = Set([m])
  A_new = setdiff(A,M)
  B = Set()   
  for i in A_new
    if Zv[i]%k == 0 
      push!(B,i)
    end 
  end   # dies kann man noch optimieren (readout Z_valuse.%k == 0 und erstelle A )
  L = []
  for C in Hecke.subsets(B)
    D = union(C,M)
    #if !(min in A)  # fixed one in recursion
    #  continue
    #end
    if v == sum([Zv[i] for i in D])
      push!(L,D)
    end 
  end 
  return L
end 



function read_path_to_leaf(CB::compressed_Blocksystem,id)# where id is the leaf id
  #return all blocksystems along the unique path to the leaf with given id
  d = CB.d
  #CB.Cyc
  B = []
  for ids in CB.T.nodes[id].path[2:end] # path without the root
    push!(B,all_choice_blocks(CB,ids))
  end 
  BL = [[b] for b in B[1]]
  BL2 = []
  for i = 2:length(B)
    for bl in BL
      for bli in B[i]
        push!(BL2,vcat(bl,[bli]))
      end 
    end 
    BL = BL2
  end 
  return BL 
end 

################################################################################################################################
### TODO PART 
function all_choice_blocks(CB,ids)
  # Cyc = CB.Cyc
  # k,usedC = CB.T.nodes[ids].values
  # for used_cycle in CB.T.nodes[ids].values[2]
  n_blocks = length(CB.T.nodes[ids].values[2]) # je n/k aus jedem block 
  # TODO
  E = []
  for j = 1:length(CB.T.nodes[ids].values[2])
    choose = divexavt(CB.T.nodes[ids].values[2][k],CB.T.nodes[ids].values[1]) # na_ / k , how many elements from each cycle we can choose to get a valid block
    if j == 1 
      for set in Hecke.subsets(CB.Cyc[j],choose)
        push!(E,set)
      end 
    end 
  end 
# TODO
end 
#=
Also TODO:
- intersect of blocksys,
- other pseudotests,
- ...
=#
################################################################################################################################
# TEST PART
function timetest(TESTFIELDS)
  #where Testfields is a List of defining polynomial
  Zx, x = ZZ["x"];
  for polynom in TESTFIELDS
    d = discriminant(polynom[2])
    p = 17 
    while d%p == 0
      p = next_prime(p)
    end # avoid p | disc(f)
    C = Oscar.GaloisGrp.GaloisCtx(polynom[2], p);
    K, _ = number_field(C.f);
    S = Oscar.GaloisGrp.SubfieldLattice(K, C);
    L = collect_subfields(C,false)
  end 
end 
function timetest2(ll)
  for b in ll
    g = galois_group(b[2])[1];
    #tg = transitive_identification(g);
    #si = signature(b[2])[1];
    #fil = open("bal$(tg)st$(si)", "a")
    #println(fil, "[",b[1],",",b[2],"],")
    #close(fil)
  end 
>>>>>>> Stashed changes
end
