#!/usr/bin/env python3
#############################################################################
# Usage:
#     ./release_notes.py [VERSION]
#
# For example
#     ./release_notes.py 4.13.1
#
# This assumes that the tags named v4.13.1, 4.13dev (?) and v4.13.0 (???) already exists.
#
# A version ending in .0 is consider MAJOR, any other MINOR
# Don't use this with versions like 4.13.0-beta1

import json
import os
import subprocess
import sys
from datetime import datetime
from typing import Any, Dict, List


ownpath = os.path.abspath(sys.argv[0])
dirpath = os.path.dirname(ownpath)
repopath = os.path.dirname(os.path.dirname(os.path.dirname(ownpath)))
newfile = f"{dirpath}/new.md"
finalfile = f"{repopath}/CHANGELOG.md"

def usage(name: str) -> None:
    print(f"Usage: `{name} [NEWVERSION] [OLDVERSION]`")
    sys.exit(1)


def is_existing_tag(tag: str) -> bool:
    print(tag)
    res = subprocess.run(
        [
            "gh",
            "release",
            "list",
            "--json=name",
            "-q",
            f""".[] | select(.name | contains("{tag.strip()}"))"""
        ],
        shell=False,
        capture_output=True
    )
    return res.stdout.decode() != ""


def find_previous_version(version: str) -> str:
    major, minor, patchlevel = map(int, version.split("."))
    if major != 1:
        error("unexpected OSCAR version, not starting with '1.'")
    if patchlevel != 0:
        patchlevel -= 1
        return f"{major}.{minor}.{patchlevel}"
    minor -= 1
    patchlevel = 0
    while True:
        v = f"{major}.{minor}.{patchlevel}"
        if not is_existing_tag("v" + v):
            break
        patchlevel += 1
    if patchlevel == 0:
        error("could not determine previous version")
    patchlevel -= 1
    return f"{major}.{minor}.{patchlevel}"

def notice(s):
    print(s)

def error(s):
    print(s)
    exit()

def warning(s):
    print('===================================================')
    print(s)
    print('===================================================')

# the following is a list of pairs [LABEL, DESCRIPTION]; the first entry is the name of a GitHub label
# (be careful to match them precisely), the second is a headline for a section the release notes; any PR with
# the given label is put into the corresponding section; each PR is put into only one section, the first one
# one from this list it fits in.
# See also <https://github.com/gap-system/gap/issues/4257>.
prioritylist = [
    ["release notes: highlight", "Highlights"],

    ["renaming", "Renamings"],

    ["topic: algebraic geometry",   "Algebraic Geometry"],
    ["topic: combinatorics",        "Combinatorics"],
    ["topic: commutative algebra",  "Commutative Algebra"],
    ["topic: FTheoryTools",         "F-Theory Tools"],
    ["topic: groups",               "Groups"],
    ["topic: lie theory",           "Lie Theory"],
    ["topic: number theory",        "Number Theory"],
    ["topic: polyhedral geometry",  "Polyhedral Geometry"],
    ["topic: toric geometry",       "Toric Geometry"],
    ["topic: tropical geometry",    "Tropical Geometry"],

    ["serialization",               "Changes related to serializing data in the MRDI file format"],

    ["enhancement",                 "New features or extended functionality"],
    ["experimental",                "Only changes experimental parts of OSCAR"],
    ["optimization",                "Performance improvements or improved testing"],
    ["bug: crash",                  "Fixed bugs that could lead to crashes"],
    ["bug",                         "Other fixed bugs"],
    ["documentation",               "Improvements or additions to documentation"],

    ["package: AbstractAlgebra",    "Changes related to the package AbstractAlgebra"],
    ["package: AlgebraicSolving",   "Changes related to the package AlgebraicSolving"],
    ["package: GAP",                "Changes related to the package GAP"],
    ["package: Hecke",              "Changes related to the package Hecke"],
    ["package: Nemo",               "Changes related to the package Nemo"],
    ["package: Polymake",           "Changes related to the package Polymake"],
    ["package: Singular",           "Changes related to the package Singular"],

]


def get_tag_date(tag: str) -> str:
    if is_existing_tag(tag):
        res = subprocess.run(
            [
                "gh",
                "release",
                "view",
                f"{tag}",
                "--json=createdAt"
            ],
            shell=False,
            capture_output=True
        )
        res = json.loads(res.stdout.decode())
    else:
        error("tag does not exist!")
    return res['createdAt'][0:10]


def get_pr_list(date: str, extra: str) -> List[Dict[str, Any]]:
    query = f'merged:>={date} -label:"release notes: not needed" -label:"release notes: added" base:master {extra}'
    print("query: ", query)
    res = subprocess.run(
        [
            "gh",
            "pr",
            "list",
            "--search",
            query,
            "--json",
            "number,title,closedAt,labels,mergedAt",
            "--limit",
            "200",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    jsonList = json.loads(res.stdout.strip())
    jsonList = sorted(jsonList, key=lambda d: d['number']) # sort by ascending PR number
    return jsonList


def pr_to_md(pr: Dict[str, Any]) -> str:
    """Returns markdown string for the PR entry"""
    k = pr["number"]
    title = pr["title"]
    return f"- [#{k}](https://github.com/oscar-system/Oscar.jl/pull/{k}) {title}\n"


def has_label(pr: Dict[str, Any], label: str) -> bool:
    return any(x["name"] == label for x in pr["labels"])


def changes_overview(
    prs: List[Dict[str, Any]], startdate: str, new_version: str
) -> None:
    """Writes files with information for release notes."""

    date = datetime.now().strftime("%Y-%m-%d")
    release_url = f"https://github.com/oscar-system/Oscar.jl/releases/tag/v{new_version}"

    # Could also introduce some consistency checks here for wrong combinations of labels
    notice("Writing release notes into file " + newfile)
    with open(newfile, "w", encoding="utf-8") as relnotes_file:
        prs_with_use_title = [
            pr for pr in prs if has_label(pr, "release notes: use title")
        ]
        # Write out all PRs with 'use title'
        relnotes_file.write(
            f"""# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
tries to adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [{new_version}]({release_url}) - {date}

The following gives an overview of the changes compared to the previous release. This list is not
complete, many more internal or minor changes were made, but we tried to only list those changes
which we think might affect some users directly.

"""
        )
        totalPRs = len(prs)
        print(f"Total number of PRs: {totalPRs}")
        countedPRs = 0
        for priorityobject in prioritylist:
            matches = [
                pr for pr in prs_with_use_title if has_label(pr, priorityobject[0])
            ]
            print("PRs with label '" + priorityobject[0] + "': ", len(matches))
            countedPRs = countedPRs + len(matches)
            if len(matches) == 0:
                continue
            relnotes_file.write("### " + priorityobject[1] + "\n\n")
            for pr in matches:
                relnotes_file.write(pr_to_md(pr))
                prs_with_use_title.remove(pr)
            relnotes_file.write("\n")
        print(f"Remaining PRs: {totalPRs - countedPRs}")
        # The remaining PRs have no "kind" or "topic" label from the priority list
        # (may have other "kind" or "topic" label outside the priority list).
        # Check their list in the release notes, and adjust labels if appropriate.
        if len(prs_with_use_title) > 0:
            relnotes_file.write("### Other changes\n\n")
            for pr in prs_with_use_title:
                relnotes_file.write(pr_to_md(pr))
            relnotes_file.write("\n")

        # Report PRs that have to be updated before inclusion into release notes.
        prs_to_be_added = [pr for pr in prs if has_label(pr, "release notes: to be added")]
        if len(prs_to_be_added) > 0:
            relnotes_file.write("### **TODO** release notes: to be added" + "\n\n")
            relnotes_file.write(
                "If there are any PRs listed below, check their title and labels.\n"
            )
            relnotes_file.write(
                'When done, change their label to "release notes: use title".\n\n'
            )
            for pr in prs_to_be_added:
                relnotes_file.write(pr_to_md(pr))
            relnotes_file.write("\n")

        # remove PRs already handled earlier
        prs = [pr for pr in prs if not has_label(pr, "release notes: to be added")]
        prs = [pr for pr in prs if not has_label(pr, "release notes: added")]
        prs = [pr for pr in prs if not has_label(pr, "release notes: use title")]

        # Report PRs that have neither "to be added" nor "added" or "use title" label
        if len(prs) > 0:
            relnotes_file.write("### **TODO** Uncategorized PR" + "\n\n")
            relnotes_file.write(
                "If there are any PRs listed below, either apply the same steps\n"
            )
            relnotes_file.write(
                'as above, or change their label to "release notes: not needed".\n\n'
            )
            for pr in prs:
                relnotes_file.write(pr_to_md(pr))
            relnotes_file.write('\n')

        # now read back the rest of changelog.md into newfile
        with open(finalfile, 'r') as oldchangelog:
            oldchangelog.seek(262)
            for line in oldchangelog.readlines():
                relnotes_file.write(line)
        # finally copy over this new file to changelog.md
        os.rename(newfile, finalfile)


def main(new_version: str, old_version: str = "") -> None:
    major, minor, patchlevel = map(int, new_version.split("."))
    extra = ""
    release_type = 0 # 0 by default, 1 for point release, 2 for patch release
    if major != 1:
        error("unexpected OSCAR version, not starting with '1.'")
    if patchlevel == 0:
        # "major" OSCAR release which changes just the minor version
        release_type = 1
        previous_minor = minor - 1
        basetag = f"v{major}.{minor}-dev"
        # *exclude* PRs backported to previous stable-1.X branch
        extra = f'-label:"backport {major}.{previous_minor}.x"'
    else:
        # "minor" OSCAR release which changes just the patchlevel
        release_type = 2
        previous_patchlevel = patchlevel - 1
        basetag = f"v{major}.{minor}.{previous_patchlevel}"
        # *include* PRs backported to current stable-4.X branch
        extra = f'label:"backport {major}.{minor}.x"'

    if old_version != "":
        basetag = f"v{old_version}"
    if release_type == 2:
        startdate = get_tag_date(basetag)
    else:
        # Find the date when the last version bump happened
        # step 1 : find the PR number of the version bump
        l = subprocess.run([
            "gh",
            "pr",
            "list",
            f'--search="Version 1.{minor}.0-dev"',
            "--state=merged",
            "--json=id,title,number,closedAt",
        ], shell=False, capture_output=True)
        j = list(json.loads(l.stdout.decode()))
        startdate = [k['closedAt'] for k in j if f"1.{minor}.0-dev" in k['title'].lower()][0][0:10]
    print("Base tag is", basetag)
    print("Base tag was closed at ", startdate)

    print("Downloading filtered PR list")
    prs = get_pr_list(startdate, extra)
    # print(json.dumps(prs, sort_keys=True, indent=4))

    # reset changelog file to state tracked in git
    
    subprocess.run(f'git checkout -- {finalfile}'.split(), check=True)

    changes_overview(prs, startdate, new_version)


if __name__ == "__main__":
    # the argument is the new version
    if len(sys.argv) == 1:
        itag = subprocess.run(
            [
                "gh",
                "release",
                "list",
                "--json=name,isLatest",
                "-q",
                ".[] | select(.isLatest == true)"
            ],
            shell=False,
            capture_output=True
        )
        itag = json.loads(itag.stdout.decode())["name"][1:]
        itag = itag.split('.')
        itag[-1] = str(int(itag[-1])+1)
        itag = ".".join(itag)
        main(itag)
    elif len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    elif len(sys.argv) >3:
        usage(sys.argv[0])
    else:
        main(sys.argv[1])
