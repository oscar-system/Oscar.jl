function symplectic_basis_of_homology(O::Origami)
    return GAP.Globals.SymplecticBasisOfHomology(GapObj(O))
end

function has_spin_structure(O::Origami)
    return GAP.Globals.HasSpinStructure(GapObj(O))::Bool
end

function spin_parity(O::Origami)
    return GAP.Globals.SpinParity(GapObj(O))::Int
end
