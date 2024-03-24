function homology(o::Origami)
    return GAP.gap_to_julia(GAP.Globals.HomologyOrigami(GapObj(o)))
end

function non_taut_part_of_homology(o::Origami, homology_basis)
    return GAP.gap_to_julia(
        GAP.Globals.NonTautPartOfHomologyOrigami(GapObj(o), GAP.julia_to_gap(homology_basis)))
end

function action_of_matrix_on_homology(o::Origami, a::Matrix)
    
end

function action_of_matrix_on_non_taut(o::Origami, a::Matrix)

end

function shadow_veech_group(o::Origami)
    shadowVeech = GAP.Globals.ShadowVeechGroup(GapObj(o))
    generators = GAP.Globals.GeneratorsOfGroup(shadowVeech)
    return matrix_group(GAP.gap_to_julia(generators))
    # need way to get julia vector from gap
end

function homology_to_string(h)

end