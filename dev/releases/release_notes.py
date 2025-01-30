#!/usr/bin/env python3
#############################################################################
# Usage:
#     ./release_notes.py VERSION
#
# For example
#     ./release_notes.py 4.13.1
#
# This assumes that the tags named v4.13.1, 4.13dev (?) and v4.13.0 (???) already exists.
#
# A version ending in .0 is consider MAJOR, any other MINOR
# Don't use this with versions like 4.13.0-beta1

import gzip
import json
import os
import subprocess
import sys
from datetime import datetime
from typing import Any, Dict, List, TextIO

import requests


def usage(name: str) -> None:
    print(f"Usage: `{name} NEWVERSION`")
    sys.exit(1)


def is_existing_tag(tag: str) -> bool:
    res = subprocess.run(
        ["git", "show-ref", "--quiet", "--verify", "refs/tags/" + tag], check=False
    )
    return res.returncode == 0


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
    ["enhancement", "New features or extended functionality"],
    ["experimental", "Only changes experimental parts of OSCAR"],
    ["optimization", "Performance improvements or improved testing"],
    ["package: AbstractAlgebra", "Changes related to the package AbstractAlgebra"],
    ["package: AlgebraicSolving", "Changes related to the package AlgebraicSolving"],
    ["package: GAP", "Changes related to the package GAP"],
    ["package: Hecke", "Changes related to the package Hecke"],
    ["package: Nemo", "Changes related to the package Nemo"],
    ["package: Polymake", "Changes related to the package Polymake"],
    ["package: Singular", "Changes related to the package Singular"],
    ["renaming", "Items being renamed ?"],
    ["serialization", "Changes related to serializing data in the MRDI file format ?"],
    ["topic: algebraic geometry", "Changes related to Algebraic Geometry"],
    ["topic: combinatorics", "Changes related to Combinatorics"],
    ["topic: FTheoryTools","Changes related to F-Theory Tools"],
    ["topic: groups","Changes related to Groups"],
    ["topic: LieTheory","Changes related to Lie Theory"],
    ["topic: number theory","Changes related to Number Theory"],
    ["topic: polyhedral geometry","Changes related to Polyhedral Geometry"],
    ["topic: rings","Changes related to Rings"],
    ["topic: schemes","Changes related to Schemes"],
    ["topic: toric schemes","Changes related to Toric Schemes"],
    ["topic: toric varieties","Changes related to Toric Varieties"],
    ["topic: tropical geometry","Changes related to Tropical Geometry"],
    ["bug: crash", "Fixed bugs that could lead to crashes"],
    ["bug", "Other fixed bugs"],
    ["documentation", "Improvements or additions to documentation"],
]


def get_tag_date(tag: str) -> str:
    if is_existing_tag(tag):
        res = subprocess.run(
            ["git", "for-each-ref", "--format=%(creatordate:short)", "refs/tags/" + tag],
            check=True,
            capture_output=True,
            text=True,
        )
    else:
        error("tag does not exist!")
    return res.stdout.strip()


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
    return json.loads(res.stdout.strip())


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

    month = datetime.now().month
    year = datetime.now().year
    day = datetime.now().day

    # Could also introduce some consistency checks here for wrong combinations of labels
    newfile = './new.md'
    finalfile = '../../CHANGELOG.md'
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

## [{new_version}] - {year}-{month}-{day}

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
        relnotes_file.write("### " + "release notes: to be added" + "\n\n")
        relnotes_file.write(
            "If there are any PRs listed below, check their title and labels.\n"
        )
        relnotes_file.write(
            'When done, change their label to "release notes: use title".\n\n'
        )

        for pr in prs:
            if has_label(pr, "release notes: to be added"):
                relnotes_file.write(pr_to_md(pr))

        prs = [pr for pr in prs if not has_label(pr, "release notes: to be added")]

        relnotes_file.write("\n")

        # Report PRs that have neither "to be added" nor "added" or "use title" label
        relnotes_file.write("### Uncategorized PR" + "\n\n")
        relnotes_file.write(
            "If there are any PRs listed below, either apply the same steps\n"
        )
        relnotes_file.write(
            'as above, or change their label to "release notes: not needed".\n\n'
        )

        for pr in prs:
            # we need to use both old "release notes: added" label and
            # the newly introduced in "release notes: use title" label
            # since both label may appear in GAP 4.12.0 changes overview
            if not (
                has_label(pr, "release notes: added")
                or has_label(pr, "release notes: use title")
            ):
                relnotes_file.write(pr_to_md(pr))
        # now read back the rest of changelog.md into newfile
        relnotes_file.write('\n')
        with open(finalfile, 'r') as oldchangelog:
            oldchangelog.seek(262)
            for line in oldchangelog.readlines():
                relnotes_file.write(line)
        # finally copy over this new file to changelog.md
        os.rename(newfile, finalfile)


def main(new_version: str) -> None:
    major, minor, patchlevel = map(int, new_version.split("."))
    if major != 1:
        error("unexpected OSCAR version, not starting with '1.'")
    if patchlevel == 0:
        # "major" OSCAR release which changes just the minor version
        previous_minor = minor - 1
        basetag = f"v{major}.{minor}dev"
        # *exclude* PRs backported to previous stable-1.X branch
        extra = f'-label:"backport-to-{major}.{previous_minor}-DONE"'
    else:
        # "minor" OSCAR release which changes just the patchlevel
        previous_patchlevel = patchlevel - 1
        basetag = f"v{major}.{minor}.{previous_patchlevel}"
        # *include* PRs backported to current stable-4.X branch
        extra = f'label:"backport-to-{major}.{minor}-DONE"'
        extra = ''

    print("Base tag is", basetag)

    startdate = get_tag_date(basetag)
    print("Base tag was created ", startdate)

    print("Downloading filtered PR list")
    prs = get_pr_list(startdate, extra)
    # print(json.dumps(prs, sort_keys=True, indent=4))

    # reset changelog file to state tracked in git
    
    subprocess.run('git checkout -- ../../CHANGELOG.md'.split(), check=True)

    changes_overview(prs, startdate, new_version)


if __name__ == "__main__":
    # the argument is the new version
    if len(sys.argv) != 2:
        usage(sys.argv[0])

    main(sys.argv[1])
