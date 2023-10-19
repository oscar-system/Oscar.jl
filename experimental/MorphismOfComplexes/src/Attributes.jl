domain(moc::MorphismOfComplexes) = moc.domain
codomain(moc::MorphismOfComplexes) = moc.codomain

function degree(moc::MorphismOfComplexes)
  return moc.degree
end

function has_right_bound(moc::MorphismOfComplexes)
  return isdefined(moc, :right_bound)
end

function has_left_bound(moc::MorphismOfComplexes)
  return isdefined(moc, :left_bound)
end

function left_bound(moc::MorphismOfComplexes)
  has_left_bound(moc) || error("morphism of complexes has no known bound to the left")
  return moc.left_bound
end

function right_bound(moc::MorphismOfComplexes)
  has_right_bound(moc) || error("morphism of complexes has no known bound to the right")
  return moc.right_bound
end

function extends_left(moc::MorphismOfComplexes)
  return moc.extends_left
end

function extends_right(moc::MorphismOfComplexes)
  return moc.extends_right
end

function Base.range(moc::MorphismOfComplexes)
  (has_right_bound(moc) && has_left_bound(moc)) || error("no bounds known")
  return left_bound(moc):right_bound(moc)
end

# Some missing methods for other things
function has_right_bound(c::ComplexOfMorphisms) 
  is_complete(c) && return true
  if typ(c) == :chain
    return !isdefined(c, :fill)
  end
  # must be a cochain complex which does not know about filling to the right.
  return true
end

function has_left_bound(c::ComplexOfMorphisms) 
  is_complete(c) && return true
  if typ(c) == :cochain
    return !isdefined(c, :fill)
  end
  # must be a chain complex which does not know about filling to the left.
  return true
end

function left_bound(c::ComplexOfMorphisms)
  has_left_bound(c) || error("complex has no known bound to the left")
  return (typ(c) == :cochain ? first(range(c)) : last(range(c)))
end

function right_bound(c::ComplexOfMorphisms) 
  has_right_bound(c) || error("complex has no known bound to the right")
  return (typ(c) == :chain ? first(range(c)) : last(range(c)))
end
