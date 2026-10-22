# DEA Machine Learning for Stata

Stata software for the article:

> Lee, C. Forthcoming. *dea_ml: Machine-learning frontier estimators for efficiency analysis in Stata*.
> Accepted for publication in the [Stata Journal](https://www.stata-journal.com/).

Twelve frontier estimators can be used directly or through `dea_ml`.
This repository maintains the article's software and subsequent author updates.
Its initial computational baseline is the curated, verified R8 submission dated
September 16, 2026. Help examples and contact details have been updated separately.
The acceptance notice refers to the article; it does not extend the article's
validation claims to later software revisions.

## Start using the commands

Requires **Stata 16.1 or newer**. Download or clone this repository and keep
**all files in `ado/` together**. Add that folder to Stata's search path:

```stata
adopath ++ "/your/path/DEA_Machine-Learning/ado"
which dea_ml
help dea_ml
```

Replace the path with the location on your computer. Use `which command_name`
to check which copy Stata loads when another installation is present.
This version uses a local-folder installation; it does not provide `net install`
metadata. While the repository is private, access requires permission.

Each help file contains reproducible synthetic examples. Copy the full Examples
block into the Do-file Editor, omit leading dots, and run from `preserve` through
`restore`. The original dataset is restored on successful completion; estimation
results and settings may change. Small samples, short boosting runs and small
tuning grids demonstrate syntax rather than an empirical research design.
No example dataset or separate do-file download is required.

## Choose a command

Open a help link below to read syntax and examples on GitHub. In Stata,
`help command_name` displays the formatted help.

| Family | Command | Purpose | Files |
| --- | --- | --- | --- |
| Interface | `dea_ml` | Dispatch to any of the twelve estimators; optional frontier graph | [code](ado/dea_ml.ado) · [help](ado/dea_ml.sthlp) |
| Trees and forests | `eatree` | Efficiency analysis trees | [code](ado/eatree.ado) · [help](ado/eatree.sthlp) |
| Trees and forests | `ceat` | Convex efficiency analysis trees | [code](ado/ceat.ado) · [help](ado/ceat.sthlp) |
| Trees and forests | `rfeat` | Random forests for efficiency analysis | [code](ado/rfeat.ado) · [help](ado/rfeat.sthlp) |
| Boosting | `beateat` | Boosted efficiency analysis trees | [code](ado/beateat.ado) · [help](ado/beateat.sthlp) |
| Boosting | `lsbmafs` | Least-squares boosting with multivariate adaptive frontier splines | [code](ado/lsbmafs.ado) · [help](ado/lsbmafs.sthlp) |
| Support-vector frontiers | `svfeat` | Support-vector frontier with EAT efficiency calculation | [code](ado/svfeat.ado) · [help](ado/svfeat.sthlp) |
| Support-vector frontiers | `csvf` | Convexified support-vector frontier | [code](ado/csvf.ado) · [help](ado/csvf.sthlp) |
| Support-vector frontiers | `svfspline` | Support-vector frontier with kernel splines and prediction | [code](ado/svfspline.ado) · [help](ado/svfspline.sthlp) |
| By-production | `bpeat` | By-production efficiency analysis trees | [code](ado/bpeat.ado) · [help](ado/bpeat.sthlp) |
| By-production | `bpceat` | Convex by-production efficiency analysis trees | [code](ado/bpceat.ado) · [help](ado/bpceat.sthlp) |
| Dynamic | `dyneat` | Dynamic efficiency analysis trees | [code](ado/dyneat.ado) · [help](ado/dyneat.sthlp) |
| Dynamic | `dynceat` | Convex dynamic efficiency analysis trees | [code](ado/dynceat.ado) · [help](ado/dynceat.sthlp) |

Scores depend on the estimator: static radial models report `e(mean_phi)`,
by-production models `e(mean_b)`, and dynamic models `e(mean_beta)`. RFEAT can
produce radial scores below one. Directional-distance scores and radial ratios
have different interpretations; consult the method-specific help before making
comparisons.

## Required helpers and optional Python backend

Keep these helpers beside the commands:

- [`_dea_ml_frontier_graph.ado`](ado/_dea_ml_frontier_graph.ado) supports `dea_ml, frontiergraph`.
- [`svfspline_p.ado`](ado/svfspline_p.ado) supports `predict` after `svfspline`.
- [`svfspline_highs_bridge.py`](ado/svfspline_highs_bridge.py) supports the optional `solver(highs)` route.

All default estimators run natively in Stata. Python is needed only for an
explicit `svfspline, solver(highs) pythonexec("/path/to/python")` request.
That interpreter needs NumPy, CVXPY and HiGHS (`highspy`). This route supports
fixed-hyperparameter fits; see [SVF-Spline help](ado/svfspline.sthlp) for its
restrictions and numerical diagnostics. R is not needed to use this package.

Some native SVF-Spline specifications return `e(numerical_warning)==1` even
when scores are produced. Inspect the diagnostics documented in the help; a
successful command return alone does not establish numerical adequacy.

## Versions and feedback

Command versions are recorded in `.ado` headers. Help revisions have a separate
dated marker; documentation-only changes do not change the computational version.
Use a Git commit to identify the exact combination of files reviewed or tested.

Report problems through [GitHub Issues](https://github.com/renegade68/DEA_Machine-Learning/issues).
Include the command name, Git commit, Stata version, `which command_name` output,
a small reproducible example, and expected versus actual results. For SVF-Spline,
include the returned numerical diagnostics and solver. Use synthetic data when
reporting a problem with confidential research data.

## Attribution and license

The software is licensed under the [Academic Free License 3.0](LICENSE), following
the author's DEA software licensing policy. See [the software notice](LICENSE_CODE.txt).
Preserve author and method references, cite the article for this software, and
cite the method papers documented in each command's help. The full journal
citation and DOI will be added when available.

This code repository contains only the selected commands, help and runtime
helpers. Manuscripts, publishing agreements, research data, validation archives,
ZIP submissions and experimental versions remain outside the repository.
The software license does not cover article manuscripts or publisher PDFs.

Contact: Choonjoo Lee — bloom.rampike@gmail.com
