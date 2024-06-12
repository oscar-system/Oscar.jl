@attributes mutable struct GroupRecognitionNode{T <: GAPGroup}
  input_group::T
  gap_node::GapObj

  function GroupRecognitionNode{T}(G::GAPGroup, gap_node::GapObj) where T
     res = new{T}(G, gap_node)
     return res
  end
end

GapObj(node::GroupRecognitionNode) = node.gap_node

# did the recognition not fail?
is_ready(node::GroupRecognitionNode) = GAP.Globals.IsReady(node.gap_node)

# does the node describe a leaf?
is_leaf(node::GroupRecognitionNode) = GAP.Globals.IsLeaf(node.gap_node)

# the group from which `node` was constructed
input_group(node::GroupRecognitionNode) = node.input_group

# the group stored in `node`, provided that recognition did not fail;
# the `GapObj`s of its elements store memory information
function group(node::GroupRecognitionNode)
  @req is_ready(node) "recognition failed, no group stored"
  return _oscar_group(GAP.Globals.Grp(node.gap_node))
end

function Base.show(io::IO, node::GroupRecognitionNode)
  if is_terse(io)
    print(io, "Recognition node of a group")
  else
    print(io, "Recognition node of ")
    io = pretty(io)
    print(terse(io), Lowercase(), group(node))
  end
end

function Base.show(io::IO, mode::MIME"text/plain", node::GroupRecognitionNode)
  io = pretty(io)
  _show_recursive(io, mode, GapObj(node); toplevel = true)
end

function _show_recursive(io::IO, mode::MIME"text/plain", gnode::GapObj; toplevel::Bool = true)

  if GAP.Globals.IsReady(gnode)
    init = "Recognition node:"
  else
    init = "Failed recognition node:"
  end
  if GAP.hasbangproperty(gnode, :projective) &&
     GAP.getbangproperty(gnode, :projective)
    init = "$init (projective)"
  end
  if GAP.Globals.Hasfhmethsel(gnode)
    ms = GAP.Globals.fhmethsel(gnode)
    if GAP.Globals.IsRecord(ms)
      if hasproperty(ms, :successMethod)
        init = "$init $(String(getproperty(ms, :successMethod)))"
      else
        init = "$init NO STAMP"
      end
    elseif GAP.Globals.IsString(ms)
      init = "$init $String(ms)"
    end
    if GAP.hasbangproperty(gnode, :comment)
      init = "$init Comment=$(String(GAP.getbangproperty(gnode, :comment)))"
    end
  end
  if GAP.Globals.HasIsRecogInfoForSimpleGroup(gnode) && GAP.Globals.IsRecogInfoForSimpleGroup(gnode)
    init = "$init Simple"
  elseif GAP.Globals.HasIsRecogInfoForAlmostSimpleGroup(gnode) && GAP.Globals.IsRecogInfoForAlmostSimpleGroup(gnode)
    init = "$init AlmostSimple"
  end
  if GAP.Globals.HasSize(gnode)
    init = "$init Size=$(GAP.Globals.Size(gnode))"
  end
  if GAP.Globals.HasGrp(gnode) && GAP.Globals.IsMatrixGroup(GAP.Globals.Grp(gnode))
    init = "$init Dim=$(GAP.getbangproperty(gnode, :dimension)) Field=$(GAP.Globals.Size(GAP.getbangproperty(gnode, :field)))"
  end
  if toplevel
    println(io, init)
  else
    println(io, Lowercase(), init)
  end

  print(io, Indent())
  if ! GAP.Globals.IsLeaf(gnode)
    println(io, "F:")
    print(io, Indent())
    if GAP.Globals.HasImageRecogNode(gnode)
      _show_recursive(io, mode, GAP.Globals.ImageRecogNode(gnode); toplevel = false)
    else
      println(io, "has no image")
    end
    print(io, Dedent())
    println(io, "K:")
    print(io, Indent())
    if GAP.Globals.HasKernelRecogNode(gnode)
      if GAP.Globals.KernelRecogNode(gnode) === GAP.Globals.fail
        println(io, "trivial kernel")
      else
        _show_recursive(io, mode, GAP.Globals.KernelRecogNode(gnode); toplevel = false)
      end
    else
      println(io, "has no kernel")
    end
    print(io, Dedent())
  end
  print(io, Dedent())
end

# call GAP's recognition algorithm
function recognize(G::Union{PermGroup, MatrixGroup})
  res = GAP.Globals.RecognizeGroup(GapObj(G))
  T = typeof(G)
  return GroupRecognitionNode{T}(G, res)
end
