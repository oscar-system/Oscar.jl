################################################################################
# Ladder Game
#
#
#
################################################################################

module LadderGame

using Oscar
import Oscar: IntegerUnion

import Oscar: pretty, Lowercase, Indent, Dedent, terse, @show_name, @show_special

include("LadderGameTypes.jl")
include("main.jl")

end
