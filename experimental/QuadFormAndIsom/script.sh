#!/bin/sh
#

# some settings that avoid weirdness in sed when it tries to
# adapt to your locale (e.g. if your system uses German as system language)
export LANG=C
export LC_CTYPE=C
export LC_ALL=C

# Files to modify (default uses all files known to git,
# but obviously you can modify it)
FILES=$(git ls-files)

# on macOS, you may need to change the following
SED_I="sed -i"
#SED_I="gsed -i"
#SED_I="sed -i ''"


# AbstractAlgebra constructors
$SED_I -e "s;\bHecke.absolute_simple_field\b;absolute_simple_field;g" $FILES



echo DONE
