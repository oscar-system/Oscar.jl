module MatrixGroups

using GAP
using Oscar
  
function __init__()
    GAP.Globals.Reread(GAP.GapObj(joinpath(Oscar.oscardir, "experimental", "Matrix", "matrix.g")))
end

function _wrap_for_gap(m::MatrixElem)
    return GAP.Globals.MakeJuliaMatrixRep(m)
end

end #module MatrixGroups
