using Pkg
Pkg.add(name="MixedSubdivisions", version="1.1"; io=devnull)

import LibGit2
repo = LibGit2.clone("https://github.com/isaacholt100/generic_root_count", "generic_root_count");
LibGit2.checkout!(repo, "ac17f42d0897d72e378310826b1c47db4d65df36");
