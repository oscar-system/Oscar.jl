import LoadFlint
using LoadFlint.GMP_jll
using LoadFlint.MPFR_jll
using LoadFlint.FLINT_jll
using Pkg.Artifacts

# TODO: use ARGS to specify custom build dir?
function gmp_artifact_dir()
    artifacts_toml = joinpath(dirname(dirname(Base.pathof(GMP_jll))), "StdlibArtifacts.toml")

    # If this file exists, it's a stdlib JLL and we must download the artifact ourselves
    if isfile(artifacts_toml)
        meta = artifact_meta("GMP", artifacts_toml)
        hash = Base.SHA1(meta["git-tree-sha1"])
        if !artifact_exists(hash)
            dl_info = first(meta["download"])
            download_artifact(hash, dl_info["url"], dl_info["sha256"])
        end
        return artifact_path(hash)
    end

    # Otherwise, we can just use the artifact directory given to us by GMP_jll
    return GMP_jll.find_artifact_dir()
end

function mpfr_artifact_dir()
    artifacts_toml = joinpath(dirname(dirname(Base.pathof(MPFR_jll))), "StdlibArtifacts.toml")

    # If this file exists, it's a stdlib JLL and we must download the artifact ourselves
    if isfile(artifacts_toml)
        meta = artifact_meta("MPFR", artifacts_toml)
        hash = Base.SHA1(meta["git-tree-sha1"])
        if !artifact_exists(hash)
            dl_info = first(meta["download"])
            download_artifact(hash, dl_info["url"], dl_info["sha256"])
        end
        return artifact_path(hash)
    end

    # Otherwise, we can just use the artifact directory given to us by GMP_jll
    return MPFR_jll.find_artifact_dir()
end

function flint_artifact_dir()
    artifacts_toml = joinpath(dirname(dirname(Base.pathof(FLINT_jll))), "StdlibArtifacts.toml")

    # If this file exists, it's a stdlib JLL and we must download the artifact ourselves
    if isfile(artifacts_toml)
        meta = artifact_meta("FLINT", artifacts_toml)
        hash = Base.SHA1(meta["git-tree-sha1"])
        if !artifact_exists(hash)
            dl_info = first(meta["download"])
            download_artifact(hash, dl_info["url"], dl_info["sha256"])
        end
        return artifact_path(hash)
    end

    # Otherwise, we can just use the artifact directory given to us by GMP_jll
    return FLINT_jll.find_artifact_dir()
end


const gmp_prefix = gmp_artifact_dir()
const mpfr_prefix = mpfr_artifact_dir()
const flint_prefix = flint_artifact_dir()


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
