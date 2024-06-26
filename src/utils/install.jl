# Install a GAP package via the URL of a GIT repository.
# If the package was already installed then update the installation,
# otherwise clone the repository.
# (The code is a modified copy of GAP's `InstallPackageFromGit`.)
function GAP_Packages_install(url::String; interactive::Bool = false, branch::String = "")
  # point PackageManager to GAP.jl's pkg dir
  GAP.Globals.PKGMAN_CustomPackageDir = GapObj(GAP.Packages.DEFAULT_PKGDIR[])
  # avoid info messages from PackageManager
  GAP.Globals.SetInfoLevel(GAP.Globals.InfoPackageManager, 0)

  url = GapObj(url)
  name = GAP.Globals.PKGMAN_NameOfGitRepo(url)
  name === GAP.Globals.fail && return false
  dir = GAP.Globals.Filename(GAP.Globals.Directory(GAP.Globals.PKGMAN_PackageDir()), name)

  # check for existing repository
  allinfo = GAP.Globals.PackageInfo(name)
  info = GAP.Globals.Filtered(allinfo,
             x -> GAP.Globals.StartsWith(x.InstallationPath, GAP.Globals.PKGMAN_PackageDir()))
  dirs = GAP.Globals.List(info, i -> i.InstallationPath)
  repo = GAP.Globals.Filename(GAP.Globals.List(dirs, GAP.Globals.Directory), GapObj(".git"));
  if repo !== GAP.Globals.fail
    # the package is already installed, update it
    return GAP.Globals.UpdatePackage(name, interactive)
  end

  ! GAP.Globals.PKGMAN_IsValidTargetDir(dir) && return false
  if branch == ""
    exec = GAP.Globals.PKGMAN_Exec(GapObj("."), GapObj("git"), GapObj("clone"), url, dir)
  else
    exec = GAP.Globals.PKGMAN_Exec(GapObj("."), GapObj("git"), GapObj("clone"), url, dir, GapObj("-b"), GapObj(branch))
  end

  (exec === GAP.Globals.fail || exec.code != 0) && return false
  GAP.Globals.PKGMAN_RefreshPackageInfo()

  # check for PackageInfo.g
  info = GAP.Globals.Filename(GAP.Globals.Directory(dir), GapObj("PackageInfo.g"));
  if ! GAP.Globals.IsReadableFile(info)
    if GAP.Globals.ValueOption("debug") != true
      GAP.Globals.PKGMAN_RemoveDir(dir)
    end
    return false
  end

  # install dependencies
  if GAP.Globals.PKGMAN_InstallDependencies(dir) != true
    if GAP.Globals.ValueOption("debug") != true
      GAP.Globals.PKGMAN_RemoveDir(dir)
    end
    return false
  end

  # compile, make doc, and check
#TODO: build the documentation once PackageManager sets the right current directory, 
# return GAP.Globals.PKGMAN_CheckPackage(dir)
  return true
end
