import LoadFlint
using LoadFlint.GMP_jll
using LoadFlint.MPFR_jll
using LoadFlint.FLINT_jll
using Pkg.Artifacts

# TODO: use ARGS to specify custom build dir?
function jll_artifact_dir(the_jll::Module)
    artifacts_toml = joinpath(dirname(dirname(Base.pathof(the_jll))), "StdlibArtifacts.toml")

    # If this file exists, it's a stdlib JLL and we must download the artifact ourselves
    if isfile(artifacts_toml)
        # the artifact name is always equal to the module name minus the "_jll" suffix
        name = replace(string(nameof(the_jll)), "_jll" => "")
        meta = artifact_meta(name, artifacts_toml)
        hash = Base.SHA1(meta["git-tree-sha1"])
        if !artifact_exists(hash)
            dl_info = first(meta["download"])
            download_artifact(hash, dl_info["url"], dl_info["sha256"])
        end
        return artifact_path(hash)
    end

    # Otherwise, we can just use the artifact directory given to us by GMP_jll
    return the_jll.find_artifact_dir()
end

const gmp_prefix = jll_artifact_dir(GMP_jll)
const mpfr_prefix = jll_artifact_dir(MPFR_jll)
const flint_prefix = jll_artifact_dir(FLINT_jll)


cd("antic")
run(`./configure
    --with-gmp=$(gmp_prefix)
    --with-mpfr=$(mpfr_prefix)
    --with-flint=$(flint_prefix)
    --prefix=/tmp/antic
`)

run(`make install -j$(Sys.CPU_THREADS)`)

println(
"""


Add the following lines to your .julia/artifacts/Overrides.toml

[e21ec000-9f72-519e-ba6d-10061e575a27]
Antic = "/tmp/antic"

and restart julia

""")
