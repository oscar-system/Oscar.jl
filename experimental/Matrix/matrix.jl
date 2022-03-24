module MatrixGroupExperimental

using GAP
using Oscar
  
function __init__()
    GAP.Globals.Reread(GAP.GapObj(joinpath(Oscar.oscardir, "experimental", "Matrix", "matrix.g")))
end

end #module MatrixGroupExperimental


