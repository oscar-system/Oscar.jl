# Use the @wrap macro to set up optimized versions of a bunch of frequently
# used GAP functions. So instead of writing e.g. GAP.Globals.DenominatorRat(x)
# you'd write GAPWrap.DenominatorRat(x). The former always performs a variable
# lookup in the GAP kernel and is not type stable. The latter accesses the
# underlying GAP function object directly and also has information about the
# return type.
#
# Note that the macro GAP.@wrap has a similar purpose as @gapwrap and @gappatribute
# have, but works on a much lower level on purpose. We may actually phase out
# use of @gapwrap and @gappatribute in the future.
module GAPWrap

using GAP

GAP.@wrap Coefficients(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DenominatorRat(x::Any)::GapInt
GAP.@wrap IsDoneIterator(x::Any)::Bool
GAP.@wrap IsList(x::Any)::Bool
GAP.@wrap IN(x::Any, y::Any)::Bool
GAP.@wrap NextIterator(x::GapObj)::Any
GAP.@wrap NumberColumns(x::GapObj)::GapInt
GAP.@wrap NumberRows(x::GapObj)::GapInt
GAP.@wrap NumeratorRat(x::Any)::GapInt

end
