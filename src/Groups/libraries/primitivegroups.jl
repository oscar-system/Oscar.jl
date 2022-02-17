export
    all_primitive_groups,
    number_primitive_groups,
    primitive_group


###################################################################
# Primitive groups, block system
###################################################################
"""
    number_primitive_groups(n::Int)

Return the number of primitive permutation groups of degree `n`,
up to permutation isomorphism.
"""
number_primitive_groups(n::Int) = GAP.Globals.NrPrimitiveGroups(n)

"""
    primitive_group(deg::Int, i::Int)

Return the `i`-th group in the catalogue of primitive groups of degree `n`
in GAP's Primitive Groups Library. The output is a group of type `PermGroup`.
"""
function primitive_group(deg::Int, n::Int)
   N = number_primitive_groups(deg)
   @assert n <= N "There are only $N primitive groups of degree $deg."
   return PermGroup(GAP.Globals.PrimitiveGroup(deg,n), deg)
end

"""
    all_primitive_groups(L...)

Return the list of all primitive groups (up to permutation isomorphism)
satisfying the conditions in `L`.
Here, `L` is a vector whose arguments are organized as
`L` = [ `func1`, `arg1`, `func2`, `arg2`, ... ],
and the function returns all the groups `G` satisfying the conditions
`func1(G) == arg1`, `func2(G) == arg2`, etc.
An argument can be omitted if it corresponds to the boolean value `true`.

# Examples
```jldoctest
julia> all_primitive_groups(degree, 4, isabelian)
PermGroup[]
```
returns the list of all abelian primitive groups acting on a set of order 4.

The type of the groups is `PermGroup`.
"""
function all_primitive_groups(L...)
   valid, temp = CheckValidType(L; isapg=true)
   @assert valid "Wrong type inserted"
   isargument = false                     # says if the inserted value is the argument of the previous value

   L1 = Vector(undef, length(L)+temp)
   pos = 1
   for i in 1:length(L)
      if typeof(L[i]) <: Function
         if isargument
            L1[pos] = true
            pos += 1
         end
         L1[pos] = find_index_function(L[i], true)[2]
         isargument = true
      else
         L1[pos] = GAP.julia_to_gap(L[i])
         isargument = false
      end
   pos+=1
   end
   if isargument L1[length(L1)]=true end

   K = GAP.Globals.AllPrimitiveGroups(L1...)
   return [PermGroup(x) for x in K]          # GAP.julia_to_gap(K) does not work
end
