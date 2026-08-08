## Submission

This is a resubmission addressing Konstanze Lauseker's feedback on the
prior submission:
* All acronyms in the Description field are now spelled out (mean
  absolute error, mean squared error, root mean squared error, in
  place of MAE/MSE/RMSE).
* plot.rtforecast.Rd now has a \value tag documenting its (invisible)
  return value.

While fixing these, markdown-style roxygen (`Roxygen: list(markdown =
TRUE)`) was found to have never actually been enabled in DESCRIPTION,
despite every roxygen comment in the package using markdown syntax -
so none of it was being converted to proper \code{}/\link{} Rd markup
anywhere in the package (it was rendering as literal backticks and
square brackets). This is now fixed and regenerated across all man
pages; in the process, two dangling cross-references to
functions/filenames that don't exist under those names
(epiFilter()/epiSmoother() instead of the actual epi_filter()/
epi_smoother(), and a leftover reference to a write_rtforecast()
function that was never implemented) were also found and corrected.

## Test environments

* local Windows 11, R 4.6.1
* GitHub Actions: windows-latest (release), macos-latest (release),
  ubuntu-latest (devel, release, oldrel-1)
* win-builder (release; devel upload repeatedly failed with a
  server-side FTP error unrelated to this package - R-devel is
  independently verified via R-hub's linux platform and via CRAN's
  own incoming pretest, both of which passed)
* R-hub (linux, macos-arm64, windows, donttest, nosuggests,
  ubuntu-release)
* CRAN's own incoming pretest (Windows and Debian, R-devel)

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
