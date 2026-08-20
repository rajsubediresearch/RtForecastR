## Submission

This is a patch release (0.1.1) fixing a correctness bug found after
0.1.0's acceptance: `rt_forecast()`'s `forecast_quantiles` element
(intended for weighted interval score calculations via `wis()`) was
computed from a disconnected in-sample calculation that did not
correspond to the genuine out-of-sample forecast reported in
`fit$forecast`. It now derives from the same forecast distribution, so
`forecast_quantiles` is guaranteed consistent with `fit$forecast`. A
regression test now covers this. See NEWS.md for details.

No user-facing API changes; no new dependencies.

## Test environments

* local Windows 11, R 4.6.1
* R CMD check --as-cran: 0 errors, 0 warnings, 0 notes

## Downstream dependencies

None (no reverse dependencies).
