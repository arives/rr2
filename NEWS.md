# rr2 1.1.0

* Fixed several bugs when working with different models.

# rr2 1.0.3

* Use `phyr`'s new function styles (underscore instead of dot in function names).
* Add methods for `phyr::pglmm()` outputs.
* Use `alphaWarn` information returned by `phylolm::phyloglm()` to avoid calculate R^2^s on models with large alpha values.

# rr2 1.0.2

* Added a `NEWS.md` file to track changes to the package.
* Added support for `nlme::gls()`.