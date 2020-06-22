export
    moved_points,
    number_moved_points

moved_points(x::PermGroupElem) = [y for y in GAP.gap_to_julia(GAP.Globals.MovedPoints(x.X))]
moved_points(G::PermGroup) = [y for y in GAP.gap_to_julia(GAP.Globals.MovedPoints(G.X))]
number_moved_points(x::Union{PermGroupElem,PermGroup}) = length(moved_points(x))

include("perfectgroups.jl")
include("primitivegroups.jl")
include("smallgroups.jl")
include("transitivegroups.jl")

#############################################################################

# useful methods for all functions of type all_"property"_groups

#############################################################################

# return the output of the function f and the corresponding GAP function
function find_index_function(f, isapg)
   if f == isabelian return Bool, GAP.Globals.IsAbelian
   elseif f == isalmostsimple return Bool, GAP.Globals.IsAlmostSimpleGroup
   elseif f == iscyclic return Bool, GAP.Globals.IsCyclic
   elseif f == isnilpotent return Bool, GAP.Globals.IsNilpotentGroup
   elseif f == isperfect return Bool, GAP.Globals.IsPerfectGroup
   elseif f == issimple return Bool, GAP.Globals.IsSimpleGroup
   elseif f == issolvable return Bool, GAP.Globals.IsSolvableGroup
   elseif f == issupersolvable return Bool, GAP.Globals.IsSupersolvableGroup
   elseif f == istransitive return Bool, GAP.Globals.IsTransitive
   elseif f == number_moved_points || f == degree
      if isapg
         return Union{Int, AbstractVector{Int}}, GAP.Globals.NrMovedPoints
      else
         throw(ArgumentError("Function not supported"))
      end
   elseif f == order return Union{Int, AbstractVector{Int}}, GAP.Globals.Size
   else throw(ArgumentError("Function not supported"))
   end
end

# check whether the input of all_small_group is valid (see below)
function CheckValidType(L; isapg=false)    # isapg = Is a Permutation Group
   isargument = false
   expected_type = Any
   arguments_to_add = 0                          # times the argument of a Boolean function is not specified
   
   for i in 1:length(L)
      if isargument
         if expected_type == Bool
            if typeof(L[i]) <: Function
               expected_type = find_index_function(L[i], isapg)[1]
               arguments_to_add += 1
            elseif typeof(L[i])==Bool
               isargument = false
            else
               return false
            end
         else
            if typeof(L[i]) <: expected_type
               isargument = false
            else
               return false
            end
         end
      else
         if !(typeof(L[i]) <: Function)
            return false
         else
            isargument = true
            expected_type = find_index_function(L[i], isapg)[1]
         end
      end
   end

   if isargument
      if expected_type == Bool
         arguments_to_add += 1
         return true, arguments_to_add
      else
         return false
      end
   else
      return true, arguments_to_add
   end
end

