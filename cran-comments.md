## Submission

This is a documentation-only update (0.1.2). No code changes.

* `?rt_forecast` now states the correct automatic `maxI` minimum (50; the
  help page previously said 2000, which did not match the code) and
  describes every element of the returned object, including
  `forecast_quantiles` and `call_args`.
* The vignette explains why the in-sample check is not forecast
  validation and links to a companion tutorial.

It has been about three months since the 0.1.1 release.

## Test environments

* local macOS (aarch64), R 4.6.1
* GitHub Actions: windows-latest, macos-latest, ubuntu-latest
  (devel, release, oldrel-1)
* win-builder (R-devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Downstream dependencies

None (no reverse dependencies).
