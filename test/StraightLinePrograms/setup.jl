using Oscar.StraightLinePrograms

using Oscar.AbstractAlgebra: AbstractAlgebra, Perm, SymmetricGroup, order, @perm_str

using Oscar.StraightLinePrograms: Const, Exp, Gen, Minus, Plus, LazyRec,
    Times, UniMinus, Call, pushconst!, pushop!,
    Line, Arg, constants, lines, evaluate!

using Oscar.StraightLinePrograms: pushline!, GAPStraightLine
using Oscar.StraightLinePrograms: AtlasSLProgram, AtlasLine

replstr(c) = sprint(show, MIME"text/plain"(), c)

Oscar.StraightLinePrograms.order(x::AbstractAlgebra.GroupElem) = AbstractAlgebra.order(x)
