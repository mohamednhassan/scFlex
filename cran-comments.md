## R CMD check results

0 errors | 0 warnings | 3 notes

* This is a new submission.

* The check reported:
  "unable to verify current time".
  This appears to be specific to the local checking environment.

* HTML validation was skipped because the external `tidy` command was not
  available in the local checking environment.

## Test environments

* Ubuntu 24.04.3 LTS, R 4.4.1
* GitHub Actions R-CMD-check
* R-hub Windows
* R-hub Ubuntu release

Additional R-hub checks on macOS ARM R-devel and a Clang 21 environment did
not reach the scFlex package check because dependency installation failed
before scFlex was checked.

## Downstream dependencies

There are currently no known downstream dependencies.