# JToric

[![Build Status](https://github.com/HereAround/JToric.jl/workflows/CI/badge.svg)](https://github.com/HereAround/JToric.jl/actions)
[![Coverage](https://codecov.io/gh/HereAround/JToric.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/HereAround/JToric.jl)


This package will only be functional if you configure the Gap-package CddInterface appropriately. The following steps aim to establish such a configuration.

1: In Julia

    using GAP
    gap_location = GAP.GAPROOT

2: exit Julia, then

    cd $(gap_location)/pkg
    git clone https://github.com/homalg-project/CddInterface.git
    git reset --hard b806a6f93788c38b0323f642bdb220287b9fc41e
    ./install  $(gap_location)
