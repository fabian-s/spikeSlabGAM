## Submission

Maintenance update. The registered `.C()` entry point `sampler()` is now
declared with its true argument types and a `void` return type, matching its
definition (previously `void *` arguments in the registration vs. typed
pointers and an `SEXP` return in the definition). This mismatch was flagged by
clang's `-fsanitize=kcfi` (see https://bugs.r-project.org/show_bug.cgi?id=19347).
Also fixes author names in the CITATION file and replaces `T`/`F` with
`TRUE`/`FALSE` in the R code.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

- local: Linux Mint 22.1 (Ubuntu 24.04), R 4.6.2, `--as-cran`
  (also installed with `-flto -Wlto-type-mismatch`: no type mismatches)
- GitHub Actions: ubuntu-latest (R-devel, release, oldrel-1),
  macos-latest (release), windows-latest (release)
- R-hub: linux, macos-arm64, windows (all R-devel), ubuntu-next, valgrind
  (OK, no memory errors), clang-asan, clang-ubsan (no sanitizer errors; 1 NOTE
  about the non-portable flag `-Wp,-D_FORTIFY_SOURCE=3`, which comes from the
  R-hub containers' own build configuration)
- win-builder: R-devel

## Reverse dependencies

Checked 2 reverse dependencies (countSTAR, SeBR) with r-devel/recheck:
no new problems.
