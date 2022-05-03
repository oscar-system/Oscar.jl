include("atlasgroups.jl")
include("perfectgroups.jl")
include("primitivegroups.jl")
include("smallgroups.jl")
include("transitivegroups.jl")

#############################################################################

# useful methods for all functions of type all_"property"_groups

#############################################################################

const _IntOrIntVec = Union{IntegerUnion, AbstractVector{<:IntegerUnion}}

const _group_filter_attrs = Dict{Any,Tuple{Type, GapObj, Any}}()
const _permgroup_filter_attrs = Dict{Any,Tuple{Type, GapObj, Any}}()
#const _ctbl_filter_attrs = Dict{Any,Tuple{Type, GapObj, Any}}()

function _add_bool_attr(attrs::Dict{Any,Tuple{Type, GapObj, Any}}, k, v::GapObj)
  attrs[k] = (Bool, v, true)
  attrs[!k] = (Bool, v, false)
end

function __init_group_libraries()
  props = [
    isabelian => GAP.Globals.IsAbelian,
    isalmostsimple => GAP.Globals.IsAlmostSimpleGroup,
    iscyclic => GAP.Globals.IsCyclic,
    is_duplicate_table => GAP.Globals.IsDuplicateTable,  # TODO: should be moved into a separate list
    isnilpotent => GAP.Globals.IsNilpotentGroup,
    isperfect => GAP.Globals.IsPerfectGroup,
    is_quasisimple => GAP.Globals.IsQuasisimple,
    issimple => GAP.Globals.IsSimpleGroup,
    is_sporadic_simple => GAP.Globals.IsSporadicSimple,
    issolvable => GAP.Globals.IsSolvableGroup,
    issupersolvable => GAP.Globals.IsSupersolvableGroup,
  ]

  empty!(_group_filter_attrs)
  for (k, v) in props
    _add_bool_attr(_group_filter_attrs, k, v)
  end

  _group_filter_attrs[order] = (_IntOrIntVec, GAP.Globals.Size, nothing)
  _group_filter_attrs[number_conjugacy_classes] = (_IntOrIntVec, GAP.Globals.NrConjugacyClasses, nothing)

#  _ctbl_filter_attrs = copy(_group_filter_attrs)
#  _add_bool_attr(_ctbl_filter_attrs, is_duplicate_table, GAP.Globals.IsDuplicateTable)

  copy!(_permgroup_filter_attrs, _group_filter_attrs)
  _add_bool_attr(_permgroup_filter_attrs, istransitive, GAP.Globals.IsTransitive)
  _add_bool_attr(_permgroup_filter_attrs, isprimitive, GAP.Globals.IsPrimitive)
  _permgroup_filter_attrs[number_moved_points] = (_IntOrIntVec, GAP.Globals.NrMovedPoints, nothing)
  _permgroup_filter_attrs[degree] = (_IntOrIntVec, GAP.Globals.NrMovedPoints, nothing)
  _permgroup_filter_attrs[transitivity] = (_IntOrIntVec, GAP.Globals.Transitivity, nothing)

end

# return the output of the function f and the corresponding GAP function
# TODO: use @nospecialize???
function find_index_function(f, permgroups::Bool)

   if permgroups
     haskey(_permgroup_filter_attrs, f) && return _permgroup_filter_attrs[f]
   else
     haskey(_group_filter_attrs, f) && return _group_filter_attrs[f]
   end
   throw(ArgumentError("Function not supported"))
end

# check whether the input of all_small_group is valid (see below)
function translate_group_library_args(args::Tuple; permgroups::Bool=false)
   gapargs = []
   i = 1
   while i <= length(args)
      arg = args[i]
      if arg isa Pair
         # handle e.g. `isabelian => false`
         func = arg[1]
         data = arg[2]
         i += 1
         expected_type, gapfunc, _ = find_index_function(func, permgroups)
         if data isa expected_type
            push!(gapargs, gapfunc, GAP.Obj(data))
         else
            throw(ArgumentError("bad argument $(arg[2]) for function $(func)"))
         end
      elseif arg isa Function
         # handle e.g. `isabelian` or `isabelian, false`
         func = arg
         expected_type, gapfunc, default = find_index_function(func, permgroups)
         i += 1
         if i <= length(args) && args[i] isa expected_type
            push!(gapargs, gapfunc, GAP.Obj(args[i]))
            i += 1
         elseif default !== nothing
            # no argument given: default to `true`
            push!(gapargs, gapfunc, default)
         else
            if i <= length(args)
               throw(ArgumentError("bad argument $(args[i]) for function $(func)"))
            else
               throw(ArgumentError("missing argument for function $(func)"))
            end
         end
      else
         throw(ArgumentError("expected a function or a pair, got $arg"))
      end
   end

   return gapargs
end
