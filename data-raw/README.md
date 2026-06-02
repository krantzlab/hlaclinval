# Single source of truth for the web calculators

`gen_fixtures.R` runs the `hlaclinval` package over a grid of inputs and writes
`hlaclinval_fixtures.json`. `verify_fixtures.mjs` re-implements the calculator math
in JavaScript (which must stay identical to the Quarto widgets) and checks it
against those fixtures.

Workflow whenever the math changes (in the package OR the widgets):

```sh
Rscript gen_fixtures.R hlaclinval_fixtures.json   # R package -> fixtures
node    verify_fixtures.mjs hlaclinval_fixtures.json   # widgets vs fixtures; nonzero exit on divergence
```

Wire `verify_fixtures.mjs` into CI (GitHub Actions) so any drift between the
package and the browser tools fails the build. Keep the JS functions in
`verify_fixtures.mjs` byte-for-byte equal to the math embedded in the calculator
HTML/OJS.
