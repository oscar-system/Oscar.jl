using Test, StraightLinePrograms

using AbstractAlgebra: AbstractAlgebra, Perm, SymmetricGroup, order, @perm_str

using StraightLinePrograms: Const, Exp, Gen, Minus, Plus, LazyRec,
    Times, UniMinus, Call, pushconst!, pushop!,
    Line, Arg, constants, lines, evaluate, evaluate!

using StraightLinePrograms: pushline!, GAPStraightLine
using StraightLinePrograms: AtlasSLProgram, AtlasLine

const SL = StraightLinePrograms

replstr(c) = sprint((io, x) -> show(io, "text/plain", x), c)

StraightLinePrograms.order(x::AbstractAlgebra.GroupElem) = AbstractAlgebra.order(x)
