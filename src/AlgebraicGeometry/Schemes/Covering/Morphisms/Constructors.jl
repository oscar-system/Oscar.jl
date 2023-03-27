########################################################################
# Constructors for CoveringMorphism                                    #
########################################################################

function identity_map(C::Covering)
  map_dict = IdDict{AbsSpec, AbsSpecMor}()
  for U in patches(C)
    map_dict[U] = identity_map(U)
  end
  return CoveringMorphism(C, C, map_dict, check=false)
end


