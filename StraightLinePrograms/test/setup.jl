using Test, StraightLinePrograms, AbstractAlgebra

using StraightLinePrograms: Const, Exp, Gen, Minus, Plus, Lazy,
    Times, UniMinus, pushconst!, pushop!,
    Line, Arg, constants, lines, evaluate, evaluate!

using StraightLinePrograms: pushline!, GAPStraightLine
using StraightLinePrograms: AtlasSLProgram, AtlasLine

const SL = StraightLinePrograms

replstr(c) = sprint((io, x) -> show(io, "text/plain", x), c)
