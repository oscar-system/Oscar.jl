@attributes mutable struct GroupRecognitionTree{T <: GAPGroup}
  input_group::T
  gap_tree::GapObj

  function GroupRecognitionTree{T}(G::GAPGroup, gap_tree::GapObj) where T
     @req G isa PermGroup || G isa MatrixGroup "only matrix and permutation groups are supported"
     res = new{T}(G, gap_tree)
     return res
  end
end

GapObj(tree::GroupRecognitionTree) = tree.gap_tree


"""
    is_ready(tree::GroupRecognitionTree)

Return `true` if the recognition procedure for the group of `tree`
was successful, and `false` otherwise.

# Examples
```jldoctest
julia> rec = recognize(GL(4, 2));  is_ready(rec)
true
```
"""
is_ready(tree::GroupRecognitionTree) = GAPWrap.IsReady(GapObj(tree))


# does the tree describe a leaf?
is_leaf(tree::GroupRecognitionTree) = GAPWrap.IsLeaf(GapObj(tree))

# the group from which `tree` was constructed
input_group(tree::GroupRecognitionTree) = tree.input_group

# the group stored in `tree`, provided that recognition did not fail.
# According to the Recog documentation, the `GapObj`s of its elements store
# memory information, and then it makes sense to distinguish between this
# group and the input group.
# However, this claim is apparently not correct.
function group(tree::GroupRecognitionTree)
  @req is_ready(tree) "recognition failed, no group stored"
  g = GAPWrap.Grp(GapObj(tree))
  if GAPWrap.IsPermGroup(g)
    # preserve the degree
    return permutation_group(g, degree(input_group(tree)))
  else
    # preserve the base ring
    return _as_subgroup_bare(input_group(tree), g)
  end
end


"""
    nice_gens(tree::GroupRecognitionTree)

Return the vector of generators of the group of `tree` w.r.t. which
the straight line programs for group elements computed by
[`straight_line_program`](@ref) are written.

# Examples
```jldoctest
julia> rec = recognize(GL(4, 2));  is_ready(rec)
true

julia> x = rand(group(rec));

julia> slp = straight_line_program(rec, x);

julia> evaluate(slp, nice_gens(rec)) == x
true
```
"""
@attr Vector{<:GAPGroupElem} function nice_gens(tree::GroupRecognitionTree)
  G = group(tree)
  return [group_element(G, x) for x in GAPWrap.NiceGens(GapObj(tree))]
end

function Base.show(io::IO, tree::GroupRecognitionTree)
  if is_terse(io)
    print(io, "Recognition tree of a group")
  else
    print(io, "Recognition tree of ")
    io = pretty(io)
    print(terse(io), Lowercase(), group(tree))
  end
end

function Base.show(io::IO, mode::MIME"text/plain", tree::GroupRecognitionTree)
  io = pretty(io)
  _show_recursive(io, mode, GapObj(tree); toplevel = true)
end

function _show_recursive(io::IO, mode::MIME"text/plain", gtree::GapObj; toplevel::Bool = true)

  if GAPWrap.IsReady(gtree)
    init = "Recognition tree:"
  else
    init = "Failed recognition tree:"
  end
  if GAP.hasbangproperty(gtree, :projective) &&
     GAP.getbangproperty(gtree, :projective)
    init = "$init (projective)"
  end
  if GAPWrap.Hasfhmethsel(gtree)
    ms = GAPWrap.fhmethsel(gtree)
    if GAPWrap.IsRecord(ms)
      if hasproperty(ms, :successMethod)
        init = "$init $(string(getproperty(ms, :successMethod)))"
      else
        init = "$init NO STAMP"
      end
    elseif GAPWrap.IsString(ms)
      init = "$init $String(ms)"
    end
    if GAP.hasbangproperty(gtree, :comment)
      init = "$init Comment=$(string(GAP.getbangproperty(gtree, :comment)))"
    end
  end
  if GAPWrap.HasIsRecogInfoForSimpleGroup(gtree) && GAPWrap.IsRecogInfoForSimpleGroup(gtree)
    init = "$init Simple"
  elseif GAPWrap.HasIsRecogInfoForAlmostSimpleGroup(gtree) && GAPWrap.IsRecogInfoForAlmostSimpleGroup(gtree)
    init = "$init AlmostSimple"
  end
  if GAPWrap.HasSize(gtree)
    init = "$init Size=$(string(GAPWrap.Size(gtree)))"
  end
  if GAPWrap.HasGrp(gtree) && GAPWrap.IsMatrixGroup(GAPWrap.Grp(gtree))
    init = "$init Dim=$(string(GAP.getbangproperty(gtree, :dimension))) Field=$(string(GAPWrap.Size(GAP.getbangproperty(gtree, :field))))"
  end
  if toplevel
    print(io, init)
  else
    print(io, Lowercase(), init)
  end

  if ! GAPWrap.IsLeaf(gtree)
    println(io, Indent())
    print(io, "F: ")
    if GAPWrap.HasImageRecogNode(gtree)
      _show_recursive(io, mode, GAPWrap.ImageRecogNode(gtree); toplevel = false)
      println(io, "")
    else
      println(io, "has no image")
    end
    print(io, "K: ")
    if GAPWrap.HasKernelRecogNode(gtree)
      if GAPWrap.KernelRecogNode(gtree) === GAP.Globals.fail
        print(io, "trivial kernel")
      else
        _show_recursive(io, mode, GAPWrap.KernelRecogNode(gtree); toplevel = false)
      end
    else
      print(io, "has no kernel")
    end
    print(io, Dedent())
  end
end


"""
    recognize(G::Union{PermGroup, MatrixGroup})

Return a `GroupRecognitionTree` object that describes the structure of `G`
in a recursive way.
If the recognition was successful (see [`is_ready`](@ref)) then
the result provides a membership test that is usually more efficient than
the membership test without the recognition information.

# Examples
```jldoctest
julia> recognize(symmetric_group(5))
Recognition tree: MovesOnlySmallPoints Size=120

julia> g = general_linear_group(4, 9);

julia> s = sub(g, [rand(g), rand(g)])[1];

julia> rec = recognize(s);  is_ready(rec)
true

julia> rand(s) in rec
true
```
"""
function recognize(G::Union{PermGroup, MatrixGroup})
  res = GAPWrap.RecognizeGroup(GapObj(G))
  T = typeof(G)
  return GroupRecognitionTree{T}(G, res)
end

# membership test:
# Note that we cannot simply delegate the question to GAP.
# - A matrix group element is regarded as an element of a matrix group
#   only if the `base_ring` values are equal,
#   whereas GAP uses natural embeddings of base rings.
# - For a permutation and a permutation group, we do not require the degrees
#   to be equal.
Base.in(g::PermGroupElem, tree::GroupRecognitionTree{PermGroup}) = GapObj(g) in GapObj(tree)

function Base.in(g::MatrixGroupElem{T, S}, tree::GroupRecognitionTree{MatrixGroup{T, S}}) where {T, S}
  base_ring(parent(g)) != base_ring(group(tree)) && return false
  return GapObj(g) in GapObj(tree)
end

# For convenience, we support also membership tests for matrices.
function Base.in(g::MatElem{T}, tree::GroupRecognitionTree{MatrixGroup{T, S}}) where {T, S}
  G = group(tree)
  base_ring(g) != base_ring(G) && return false
  return map_entries(_ring_iso(G), g) in GapObj(tree)
end


"""
    straight_line_program(tree::GroupRecognitionTree, g::GAPGroupElem)

Return a straight line program for the element `g` of the group of `tree`.
The inputs of this program correspond to `nice_gens(tree)`,
see [`nice_gens`](@ref) for an example.
"""
function straight_line_program(tree::GroupRecognitionTree, g::GAPGroupElem)
  return straight_line_program(GAPWrap.SLPforElement(GapObj(tree), GapObj(g)))
end
