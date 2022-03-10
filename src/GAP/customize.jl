# We want to switch off GAP's info statements by default
# when GAP is used by Oscar.
# More generally, `GAP_info_messages_off` can set all info levels either
# to zero or to those values that were valid when Oscar was started.
const __GAP_info_levels_default = Pair{GAP.GapObj, Int}[]

function __GAP_info_messages_off(off::Bool = true)
  if length(__GAP_info_levels_default) == 0
    # Initialize the info about default levels.
    for c in Vector{GAP.GapObj}(GAP.Globals.INFO_CLASSES)
      push!(__GAP_info_levels_default, c => GAP.Globals.InfoLevel(c))
    end
  end
  # Set the info levels.
  for pair in __GAP_info_levels_default
    GAP.Globals.SetInfoLevel(pair[1], off ? 0 : pair[2])
  end
end
