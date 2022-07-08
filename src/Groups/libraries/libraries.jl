
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
const _ctbl_filter_attrs = Dict{Any,Tuple{Type, GapObj, Any}}()
const _atlas_group_filter_attrs = Dict{Any,Tuple{Type, GapObj, Any}}()

function _add_bool_attr(attrs::Dict{Any,Tuple{Type, GapObj, Any}}, k, v::GapObj)
  attrs[k] = (Bool, v, true)
  attrs[!k] = (Bool, v, false)
end

function __init_group_libraries()
  props = [
    is_abelian => GAP.Globals.IsAbelian,
    is_almostsimple => GAP.Globals.IsAlmostSimpleGroup,
    is_cyclic => GAP.Globals.IsCyclic,
    is_nilpotent => GAP.Globals.IsNilpotentGroup,
    is_perfect => GAP.Globals.IsPerfectGroup,
    is_quasisimple => GAP.Globals.IsQuasisimple,
    is_simple => GAP.Globals.IsSimpleGroup,
    is_sporadic_simple => GAP.Globals.IsSporadicSimple,
    is_solvable => GAP.Globals.IsSolvableGroup,
    is_supersolvable => GAP.Globals.IsSupersolvableGroup,
  ]

  empty!(_group_filter_attrs)
  for (k, v) in props
    _add_bool_attr(_group_filter_attrs, k, v)
  end

  _group_filter_attrs[order] = (_IntOrIntVec, GAP.Globals.Size, nothing)
  _group_filter_attrs[number_conjugacy_classes] = (_IntOrIntVec, GAP.Globals.NrConjugacyClasses, nothing)

  copy!(_permgroup_filter_attrs, _group_filter_attrs)
  _add_bool_attr(_permgroup_filter_attrs, is_transitive, GAP.Globals.IsTransitive)
  _add_bool_attr(_permgroup_filter_attrs, is_primitive, GAP.Globals.IsPrimitive)
  _permgroup_filter_attrs[number_moved_points] = (_IntOrIntVec, GAP.Globals.NrMovedPoints, nothing)
  _permgroup_filter_attrs[degree] = (_IntOrIntVec, GAP.Globals.NrMovedPoints, nothing)
  _permgroup_filter_attrs[transitivity] = (_IntOrIntVec, GAP.Globals.Transitivity, nothing)

  copy!(_ctbl_filter_attrs, _group_filter_attrs)
  _add_bool_attr(_ctbl_filter_attrs, is_duplicate_table, GAP.Globals.IsDuplicateTable)

  empty!(_atlas_group_filter_attrs)
  _atlas_group_filter_attrs[degree] = (_IntOrIntVec, GAP.Globals.NrMovedPoints, nothing)
  _atlas_group_filter_attrs[rank_action] = (_IntOrIntVec, GAP.Globals.RankAction, nothing)
  _atlas_group_filter_attrs[transitivity] = (_IntOrIntVec, GAP.Globals.Transitivity, nothing)
  _add_bool_attr(_atlas_group_filter_attrs, is_transitive, GAP.Globals.IsTransitive)
  _add_bool_attr(_atlas_group_filter_attrs, is_primitive, GAP.Globals.IsPrimitive)
  _atlas_group_filter_attrs[base_ring] = (AbstractAlgebra.Ring, GAP.Globals.Ring, nothing)
  _atlas_group_filter_attrs[character] = (Oscar.GroupClassFunction, GAP.Globals.Character, nothing)
  _atlas_group_filter_attrs[characteristic] = (_IntOrIntVec, GAP.Globals.Characteristic, nothing)
  _atlas_group_filter_attrs[dim] = (_IntOrIntVec, GAP.Globals.Dimension, nothing)
end

# return the output of the function f and the corresponding GAP function
# TODO: use @nospecialize???
function find_index_function(f, filter_attrs::Dict)
   haskey(filter_attrs, f) && return filter_attrs[f]
   throw(ArgumentError("Function not supported"))
end

# check whether the input of all_small_group is valid (see below)
function translate_group_library_args(args::Tuple; filter_attrs::Dict = _group_filter_attrs)
   gapargs = []
   for arg in args
      if arg isa Pair
         # handle e.g. `is_abelian => false`
         func = arg[1]
         data = arg[2]
         expected_type, gapfunc, _ = find_index_function(func, filter_attrs)
         data isa expected_type || throw(ArgumentError("bad argument $(data) for function $(func)"))
         push!(gapargs, gapfunc, GAP.Obj(data))
      elseif arg isa Function
         # handle e.g. `is_abelian` or `! is_abelian`
         func = arg
         expected_type, gapfunc, default = find_index_function(func, filter_attrs)
         default !== nothing || throw(ArgumentError("missing argument for function $(func)"))
         push!(gapargs, gapfunc, default)
      else
         throw(ArgumentError("expected a function or a pair, got $arg"))
      end
   end

   return gapargs
end
