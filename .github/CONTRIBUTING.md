# Contributing to WithinHostExperiment

Bug reports and pull requests are welcome at
<https://github.com/CuiweiG/WithinHostExperiment/issues>. Please open
pull requests against `main`.

The rest of this page records how changes move between GitHub and
Bioconductor. The two sides use different names for the development
branch, which is the easiest place to make a mistake.

## Which branch a change goes to

| Change | GitHub branch | Bioconductor branch | Version |
|---|---|---|---|
| New features, changed defaults, refactoring, documentation for the next release | `main` | `devel` | 1.1.z (odd y) |
| Bug fixes and documentation corrections for users of Bioconductor 3.24 | `RELEASE_3_24`, after `main` | `RELEASE_3_24`, after `devel` | 1.0.z (even y) |

- GitHub `main` has the same code as Bioconductor `devel` but only our
  own commits; see Pushing.
- `RELEASE_3_24` takes bug fixes and documentation corrections only.
  New features and changed defaults wait on `main` for Bioconductor
  3.25 (spring 2027).
- A bug fix goes into `main` first and is then cherry-picked to
  `RELEASE_3_24`. Each branch gets its own version bump.
- When Bioconductor 3.25 is released, `RELEASE_3_25` replaces
  `RELEASE_3_24` as the branch for bug fixes, and `RELEASE_3_24` is
  frozen.

## Version numbers

Bioconductor uses x.y.z version numbers
(<https://contributions.bioconductor.org/versionnum.html>):

- y is odd on `devel` (1.1.z) and even on the release branch (1.0.z).
  The Bioconductor team changes y on `devel` at each release; make the
  same change on `main` in a commit of our own. Do not change x or y
  otherwise.
- Bump z in every commit that should reach users: 1.1.0 to 1.1.1 on
  `devel`, 1.0.0 to 1.0.1 on `RELEASE_3_24`. A pushed commit without a
  version bump is built, but not propagated: `BiocManager::install()`
  keeps serving the previous version.
- Commits that only touch files left out of the package by
  `.Rbuildignore` (for example `.github/`) do not need a bump.
- Builds run once a day and take about 24 hours, so a push appears on
  the build report (<https://bioconductor.org/checkResults/>) one or two
  days later. It reaches users once a build with the new version
  finishes without errors.

## Remotes

This is the setup Bioconductor recommends
(<https://contributions.bioconductor.org/git-version-control.html>):

```sh
git remote add upstream git@git.bioconductor.org:packages/WithinHostExperiment.git
git remote -v
# origin    https://github.com/CuiweiG/WithinHostExperiment.git (fetch)
# origin    https://github.com/CuiweiG/WithinHostExperiment.git (push)
# upstream  git@git.bioconductor.org:packages/WithinHostExperiment.git (fetch)
# upstream  git@git.bioconductor.org:packages/WithinHostExperiment.git (push)
```

`ssh -T git@git.bioconductor.org info` lists what the registered SSH key
may read (R) and write (W); `packages/WithinHostExperiment` should show
`R W`.

## Pushing

`main` keeps only our own commits. Bioconductor `devel` also has the
Bioconductor team's version bumps, so never merge `upstream/devel` into
`main`. Commit on `main`, then cherry-pick to `devel`:

```sh
git fetch upstream
git checkout -B bioc-devel upstream/devel
git cherry-pick <commits from main>
git push upstream bioc-devel:devel
git push origin main
```

Release (the same branch on both sides):

```sh
git fetch upstream
git checkout -b RELEASE_3_24 upstream/RELEASE_3_24   # first time only
git cherry-pick <commit on main>
# bump z in DESCRIPTION (1.0.0 to 1.0.1) and commit
git push upstream RELEASE_3_24
git push origin RELEASE_3_24
```

To push exactly the commit that was checked, push by hash:

```sh
git push upstream <hash on bioc-devel>:refs/heads/devel
git push origin <hash on main>:refs/heads/main
```

Never force-push. Bioconductor does not accept rewritten history, and a
GitHub ruleset blocks force-pushes and deletion on `main` and on
`RELEASE_*` branches.

## Continuous integration

The workflows in `.github/workflows/` run on `main` with R-devel and on
`RELEASE_*` branches with R 4.6 and Bioconductor 3.24: `R CMD check` on
Linux, macOS and Windows, `BiocCheckGitClone()` and `BiocCheck()`, test
coverage, and a check of `citation()`, links and spelling. The pkgdown
site is built from `main` only. The versions are set at the top of each
workflow; update them when Bioconductor 3.25 is released.
