## Submission

This is a resubmission. Per Uwe Ligges' feedback on the first
submission, "+ file LICENSE" has been removed from the License field
and the LICENSE file itself has been deleted; copyright-holder
attribution now lives entirely in Authors@R (Kris V. Parag as
ctb/cph for the ported EpiFilter files, Raj Subedi as aut/cre/cph for
everything else), as instructed.

## Test environments

* local Windows 11, R 4.6.1
* GitHub Actions: windows-latest (release), macos-latest (release),
  ubuntu-latest (devel, release, oldrel-1)
* win-builder (release; devel upload repeatedly failed with a
  server-side FTP error unrelated to this package - R-devel is
  independently verified via R-hub's linux platform below)
* R-hub (linux, macos-arm64, windows, donttest, nosuggests,
  ubuntu-release)

## R CMD check results

0 errors | 0 warnings | 1 note (see below)

## Downstream dependencies

This is a new package with no reverse dependencies.

## Notes for CRAN reviewers

* Under R-hub's "nosuggests" check platform, vignette re-building fails
  because rmarkdown (a Suggests package) is deliberately absent in that
  environment - this is expected for any package using a knitr::rmarkdown
  vignette engine (rmarkdown is required to knit .Rmd files at all,
  independent of any code inside the vignette) and is consistent with
  what appears on CRAN's own check pages for other published packages
  using the same vignette engine.
* The "Possibly misspelled words in DESCRIPTION" NOTE (Bracher, EpiFilter,
  Parag, al, et) are author surnames, the package/algorithm name, and
  citation shorthand ("et al.") - all correct as written.
* R/epiFilter.R, R/epiSmoother.R, and the core loop in R/recursPredict.R
  are ports of the EpiFilter algorithm (Parag KV, 2021, PLOS
  Computational Biology 17(9): e1009347), used and redistributed under
  its original GPL-3 license. This is documented in Authors@R (Kris V.
  Parag, role = "ctb, cph") and in each file's roxygen header.
* The generation-interval discretization convention in compute_lambda()
  follows GrowthPredict (Chowell et al. 2024); this is a shared
  mathematical convention, not reused code, and is cited accordingly
  in the function's documentation.
