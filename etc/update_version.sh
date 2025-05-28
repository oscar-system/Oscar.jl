#!/bin/sh
set -e

# print error in red and exit
error() {
    printf '\033[31mERROR: %s\033[0m\n' "$*"
    exit 1
}

# error if version is missing
# TODO: validate format
if test -z "$1"; then
    echo "Error, no version given"
    exit 1
fi

# release version
relvers=$1

# release date, default is today
reldate=$(date '+%d/%m/%Y')
reldate_iso=$(date '+%Y-%m-%d')

git update-index --refresh > /dev/null || error "uncommitted changes detected"
git diff-index --quiet HEAD -- || error "uncommitted changes detected"

echo "Setting version to $relvers, released $reldate_iso"

# update version in several files
perl -pi -e 's;version = "[^"]+";version = "'$relvers'";' Project.toml
perl -pi -e "s;Version [^ },]+;Version $relvers;" README.md
perl -pi -e "s;{V}ersion [^ },]+;{V}ersion $relvers;" README.md
perl -pi -e 's;Date := "[^"]+",;Date := "'$reldate'",;' gap/OscarInterface/PackageInfo.g
perl -pi -e 's;Version := "[^"]+",;Version := "'$relvers'",;' gap/OscarInterface/PackageInfo.g


# commit it
git commit -m "Version $relvers" Project.toml README.md gap/OscarInterface/PackageInfo.g
