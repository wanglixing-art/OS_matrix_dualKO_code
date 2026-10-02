# Software environment

`R_sessionInfo.txt` was generated from a live `installed.packages()` call on the machine
that produced the results, so the versions listed there are the ones the manuscript's
numbers actually came from. **Treat it as authoritative**; this README is only narrative.

## Base

| Tool | Version |
|---|---|
| R | 4.6.1 (2026-06-24 ucrt), x86_64-w64-mingw32 |
| OS | Windows 11 |
| Node.js | 22.x — used only by `01_data_download/P0_build_clinical.js` |
| Python | 3.13 — used only by the DOCX/figure assembly scripts (`python-docx`, `pandas`, `matplotlib`) |

## Installing the R dependencies

```r
source("environment/install_packages.R")
```

## Windows: you need a compiler

Three dependency chains build from source, so a toolchain is mandatory:

- `scTenifoldKnk` / `scTenifoldNet` (compiled C++ core)
- `CellChat`
- `infercnv` (via `rjags` / JAGS for the HMM)

**Rtools version must match the R major version** — R 4.5.x and R 4.6.x both use
**Rtools45** (there is no Rtools46). Verify before installing:

```r
pkgbuild::has_build_tools(debug = TRUE)   # must return TRUE
```

If R is installed outside `C:\Program Files\R`, R cannot auto-detect Rtools and you must
write `etc\Renviron.site` yourself with `RTOOLS45_HOME`, an augmented `PATH`, and the
`CC` / `CXX` / `FC` prefixes pointing at the toolchain binaries. On a default Rtools45
install the compilers live under
`<rtools>\x86_64-w64-mingw32.static.posix\bin\` — **not** under `usr\bin`, which only
holds the MSYS2 core utilities (`make`, coreutils). Pointing `PATH` at `usr\bin` is the
most common cause of "Rtools is required to build R packages" errors that persist even
after installing Rtools.

`rjags` additionally needs JAGS itself (a separate installer, not an R package), and
`make` on Windows is easiest to obtain through Rtools' MSYS2. The installer helpers used
during development are kept in `../misc/` for reference.

## Reproducibility caveat

`scTenifoldKnk` builds a single-cell gene regulatory network with **randomised** tensor
decomposition. Set a seed before running `04_dual_knockout/P5_dualKO.R` if you need
run-to-run comparability; network construction is also the slowest step in the pipeline,
so allow substantial wall time.
