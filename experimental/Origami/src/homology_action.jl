function homology(o::Origami)
    return GAP.gap_to_julia(GAP.Globals.HomologyOrigami(GapObj(o)))
end

function non_taut_part_of_homology(o::Origami, homology_basis)
    return GAP.gap_to_julia(
        GAP.Globals.NonTautPartOfHomologyOrigami(GapObj(o), GapObj(homology_basis)))
end

function action_of_matrix_on_homology(o::Origami, a::Matrix)
    # TODO does not work yet because of error in GAP package, change
    # IsMatrixObj to IsMatrix in homologyaction.gd
    return GAP.gap_to_julia(
        GAP.Globals.ActionOfMatrixOnHom(GapObj(o), GapObj(a)))
end

function action_of_matrix_on_non_taut(o::Origami, a::Matrix)
    # TODO does not work yet because of error in GAP package, change
    # IsMatrixObj to IsMatrix in homologyaction.gd
    return GAP.gap_to_julia(
        GAP.Globals.ActionOfMatrixOnNonTaut(GapObj(o), GapObj(a)))
end

function shadow_veech_group(o::Origami)
    shadowVeech = GAP.Globals.ShadowVeechGroup(GapObj(o))::GapObj
    generators = GAP.Globals.GeneratorsOfGroup(shadowVeech)::GapObj
    julia_matrices = [matrix(ZZ, gen) for gen in generators]
    return matrix_group(julia_matrices)
end

function homology_to_string(h)
    # TODO
end
