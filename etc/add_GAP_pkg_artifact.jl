#
# This script is used to update Artifacts.toml.
# Any changes to unrelated parts of Artifacts.toml (like comments) *should be reverted* after running this script.
#
# Usage variants:
#   julia --project=etc etc/add_GAP_pkg_artifact.jl name https://github.com/.../releases/download/.../....tar.gz
#

using Downloads: download
import Pkg
using Pkg.Artifacts
import TOML
import SHA

function sha256sum(tarball_path)
    return open(tarball_path, "r") do io
        return bytes2hex(SHA.sha256(io))
    end
end

function add_artifact_for_package(pkgname::String, pkgurl::String, artifact_toml = joinpath(@__DIR__, "..", "Artifacts.toml"))
  pkg_artifact_name = "GAP_pkg_$(lowercase(pkgname))"
  
  tarball_path = download(pkgurl)
  tarball_hash = sha256sum(tarball_path)
  
  pkg_artifact_hash = create_artifact() do artifact_dir
      Pkg.PlatformEngines.unpack(tarball_path, artifact_dir)
  end

  rm(tarball_path)

  bind_artifact!(
    artifact_toml,
    pkg_artifact_name,
    pkg_artifact_hash;
    download_info=[(pkgurl, tarball_hash)],
    force=true,
  )
end

@assert length(ARGS) == 2
pkgname = ARGS[1]
pkgurl = ARGS[2]
add_artifact_for_package(pkgname, pkgurl)
