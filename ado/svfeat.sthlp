{smcl}
{* Help revision: 2026-10-22.1; command version unchanged.}{...}
{* *! svfeat  v0.2.1  08sep2026}{...}
{viewerjumpto "Syntax" "svfeat##syntax"}{...}
{viewerjumpto "Description" "svfeat##description"}{...}
{viewerjumpto "Options" "svfeat##options"}{...}
{viewerjumpto "Stored results" "svfeat##results"}{...}
{viewerjumpto "Examples" "svfeat##examples"}{...}
{viewerjumpto "References" "svfeat##references"}{...}

{title:Title}
{p2colset 5 18 20 2}{...}
{p2col:{cmd:svfeat} {hline 2}}Support Vector Frontier on an equal-width grid (SVF){p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}
{p 8 16 2}
{cmd:svfeat} {ifin}{cmd:,} {opt x(varlist)} {opt y(varlist)}
[{opt gen(name)} {opt d(#)} {opt c(#)} {opt eps(#)} {opt tune}
 {opt cgrid(numlist)} {opt epsgrid(numlist)} {opt dgrid(numlist)}
 {opt folds(#)} {opt seed(#)} {opt nu(#)}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt x(varlist)}}nonnegative input variables{p_end}
{synopt:{opt y(varlist)}}nonnegative outputs; at least one positive per observation{p_end}
{syntab:Optional}
{synopt:{opt gen(name)}}new variable to store efficiency score phi (>= 1){p_end}
{synopt:{opt d(#)}}equal-width intervals per input; default {cmd:d(3)}{p_end}
{synopt:{opt nu(#)}}deprecated compatibility option; ignored{p_end}
{synopt:{opt c(#)}}slack penalty (tightness of frontier); default {cmd:c(1.0)}{p_end}
{synopt:{opt eps(#)}}epsilon-insensitivity band; default {cmd:eps(0.0)}{p_end}
{synopt:{opt seed(#)}}random seed for tuning folds; fixed fits preserve the RNG state{p_end}
{synopt:{opt tune}}select C, epsilon, and d by leakage-free cross-validation{p_end}
{synopt:{opt cgrid(numlist)}}candidate C values{p_end}
{synopt:{opt epsgrid(numlist)}}candidate epsilon values{p_end}
{synopt:{opt dgrid(numlist)}}candidate grid interval counts{p_end}
{synopt:{opt folds(#)}}number of CV folds; default 5{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:svfeat} implements model (9) of Valero-Carreras, Aparicio, and Guerrero
(2021). It estimates a non-decreasing step production frontier using an
equal-width rectangular grid. A point activates its own grid cell and every
componentwise dominated cell.{p_end}

{pstd}
For multiple outputs, {cmd:svfeat} follows the multi-output extension of
Valero-Carreras, Aparicio, and Guerrero (2022): all outputs share the same input
grid, while each output has its own coefficient vector and frontier. Conditional
on the common grid and hyperparameters, the objective and constraints are
separable across output dimensions.{p_end}

{pstd}
The frontier at observation i is estimated as:{p_end}
{phang2}f(x_i) = w' phi_G(x_i){p_end}
{pstd}
The free grid coefficients are represented as positive and negative parts. The
linear program minimizes:{p_end}
{phang2}||w||_1 + C * sum(xi){p_end}
{pstd}
subject to y_i <= f(x_i) <= y_i + eps + xi_i, xi_i >= 0, and adjacent
cumulative grid values being nondecreasing in every input dimension.{p_end}

{pstd}
The resulting efficiency score is phi = f(x_i) / y_i >= 1, where phi = 1 denotes a
fully efficient DMU on the estimated frontier. The parameter {opt c()} controls
the tradeoff between the L1 structural term and empirical error.{p_end}

{pstd}
Requires Stata 16 or later (uses the LinearProgram class).{p_end}

{marker options}{...}
{title:Options}

{phang}{opt x(varlist)} specifies the input variables. Required.{p_end}

{phang}{opt y(varlist)} specifies the output variables. Required.{p_end}

{phang}{opt gen(name)} stores the efficiency score phi (>= 1) in the new variable
{it:name}. If omitted, results are stored in {cmd:e()} only.{p_end}

{phang}{opt d(#)} sets the number of equal-width intervals per input. The grid
contains (d+1)^p coefficient cells for p inputs. Default is 3.{p_end}

{phang}{opt nu(#)} is retained for command-line compatibility but ignored.
Model (9) fixes the coefficient of the L1 norm at one.{p_end}

{phang}{opt c(#)} sets the slack penalty controlling how tightly the frontier
fits the data. Larger values produce a tighter (less forgiving) frontier.
Default is 1.0.{p_end}

{phang}{opt eps(#)} sets the width of the epsilon-insensitivity band. Observations
within {it:#} units of the frontier incur no penalty. Default is 0.0.{p_end}

{phang}{opt seed(#)} sets the random seed for fold assignment only with
{cmd:tune}. Fixed fits preserve Stata's RNG state. Default is 12345.{p_end}

{phang}{opt tune} selects C, epsilon, and d by held-out MSE. Every fold's input
ranges, grid transformation, and coefficients are constructed from its training
subset only; the selected candidate is refitted on all estimation rows.{p_end}

{phang}{opt cgrid(numlist)}, {opt epsgrid(numlist)}, and {opt dgrid(numlist)}
specify candidate values. Defaults are {cmd:.1 .5 1 2 5},
{cmd:0 .001 .01 .1 .2}, and {cmd:1 2 3}. Candidate scan order supplies the
deterministic tie rule: the command minimizes mean held-out MSE and retains
the first candidate in d--C--epsilon scan order when losses differ by no more
than 1e-14.{p_end}

{phang}{opt folds(#)} sets the number of folds. Default is 5. {opt seed(#)}
controls fold assignment only when {cmd:tune} is specified. The number of
folds may not exceed the number of complete estimation observations.{p_end}

{pstd}Because the number of grid cells is (d+1)^p, the command rejects candidate
grids exceeding 2,500 cells. This guard prevents accidental exponential memory
allocation and should not be interpreted as a methodological recommendation.{p_end}

{marker results}{...}
{title:Stored results}

{pstd}{cmd:svfeat} stores the following in {cmd:e()}:{p_end}
{synoptset 18 tabbed}{...}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(sample)}}estimation-sample indicator{p_end}
{synopt:{cmd:e(mean_phi)}}mean efficiency score across all DMUs{p_end}
{synopt:{cmd:e(n_failed)}}number of DMUs for which the LP failed{p_end}
{synopt:{cmd:e(solver)}}LP solver identifier{p_end}
{synopt:{cmd:e(n_svec)}}number of nonzero grid coefficients{p_end}
{synopt:{cmd:e(d)}}equal-width interval count{p_end}
{synopt:{cmd:e(c_pen)}}slack penalty used{p_end}
{synopt:{cmd:e(eps)}}epsilon-insensitivity band used{p_end}
{synopt:{cmd:e(objective)}}sum of output-specific model (9) objective values{p_end}
{synopt:{cmd:e(max_violation)}}maximum original inequality violation{p_end}
{synopt:{cmd:e(coefficients)}}grid coefficient matrix{p_end}
{synopt:{cmd:e(frontier)}}fitted frontier matrix{p_end}
{synopt:{cmd:e(grid_inputs)}}virtual input grid{p_end}
{synopt:{cmd:e(grid_frontier)}}fitted outputs at the virtual grid points{p_end}
{synopt:{cmd:e(tuned)}}1 if CV was requested{p_end}
{synopt:{cmd:e(folds)}}number of CV folds, or 0 for a fixed fit{p_end}
{synopt:{cmd:e(cv_failures)}}number of candidates with at least one failed fold{p_end}
{synopt:{cmd:e(cv_results)}}candidate C, epsilon, d, and mean held-out MSE; returned with {cmd:tune}{p_end}
{synopt:{cmd:e(cv_fold_results)}}candidate-by-fold held-out MSE; returned with {cmd:tune}{p_end}
{synopt:{cmd:e(cv_predictions)}}held-out predictions and identifiers; returned with {cmd:tune}{p_end}
{synopt:{cmd:e(fold_id)}}fold assignment by estimation-sample row; returned with {cmd:tune}{p_end}
{synopt:{cmd:e(x)}}names of input variables{p_end}
{synopt:{cmd:e(y)}}names of output variables{p_end}
{synopt:{cmd:e(cmd)}}"svfeat"{p_end}
{synopt:{cmd:e(cmdline)}}executed sub-command invocation{p_end}

{marker examples}{...}
{title:Examples}

{pstd}Run the complete example block in the Do-file Editor, omitting leading dots.
The examples create their own data; {cmd:preserve} and {cmd:restore} retain your
original dataset on successful completion. Settings and estimation results may change.
Small samples and tuning grids demonstrate syntax rather than a research design.{p_end}
{* DEA_ML_HELP_EXAMPLES_BEGIN}{...}
{phang2}{cmd:. preserve}{p_end}

{pstd}Create reproducible example data:{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 80}{p_end}
{phang2}{cmd:. set seed 42}{p_end}
{phang2}{cmd:. generate double x1 = 1 + 4*runiform()}{p_end}
{phang2}{cmd:. generate double x2 = 2 + 3*runiform()}{p_end}
{phang2}{cmd:. generate double x = x1}{p_end}
{phang2}{cmd:. generate double y = x1^0.4*x2^0.3*exp(-0.2*runiform())}{p_end}

{pstd}Basic SVF using the default grid:{p_end}
{phang2}{cmd:. svfeat, x(x) y(y) gen(phi_basic)}{p_end}

{pstd}Use two equal-width intervals and a larger error penalty:{p_end}
{phang2}{cmd:. svfeat, x(x1 x2) y(y) gen(phi_grid) d(2) c(5.0)}{p_end}

{pstd}Add an epsilon-insensitivity band:{p_end}
{phang2}{cmd:. svfeat, x(x) y(y) gen(phi_eps) eps(0.05)}{p_end}

{pstd}Tune C, epsilon, and d and inspect the selection table:{p_end}
{phang2}{cmd:. svfeat, x(x) y(y) gen(phi_tuned) tune cgrid(.5 1 2) epsgrid(0 .05) dgrid(1 2) folds(5) seed(2026)}{p_end}
{phang2}{cmd:. matrix list e(cv_results)}{p_end}

{pstd}Display number of nonzero grid coefficients and mean efficiency:{p_end}
{phang2}{cmd:. svfeat, x(x) y(y) gen(phi_display)}{p_end}
{phang2}{cmd:. display "Nonzero grid weights: " e(n_svec)}{p_end}
{phang2}{cmd:. display "Mean phi: " e(mean_phi)}{p_end}

{phang2}{cmd:. restore}{p_end}
{* DEA_ML_HELP_EXAMPLES_END}{...}


{marker references}{...}
{title:References}

{phang}Valero-Carreras, D., Aparicio, J., and Guerrero, N. M. (2021). Support vector
frontiers: A new approach for estimating production functions through support vector
machines. {it:Omega} 104:102490.
{browse "https://doi.org/10.1016/j.omega.2021.102490"}
{p_end}

{phang}Valero-Carreras, D., Aparicio, J., and Guerrero, N. M. (2022). Multi-output
support vector frontiers. {it:Computers & Operations Research} 143:105765.
{browse "https://doi.org/10.1016/j.cor.2022.105765"}
{p_end}

{marker author}{...}
{title:Author}

{pstd}Choonjoo Lee{break}
Korea National Defense University{break}
Nonsan-si, Republic of Korea{break}
bloom.rampike@gmail.com{p_end}

{title:Publication and license}

{pstd}Lee, C. Forthcoming. {it:dea_ml: Machine-learning frontier estimators for efficiency analysis in Stata}.
Accepted for publication in the {browse "https://www.stata-journal.com/":Stata Journal}.{p_end}
{pstd}Author-maintained software under the Academic Free License 3.0 (AFL-3.0).
See LICENSE and LICENSE_CODE.txt in the repository. The acceptance notice refers
to the article; later software revisions are maintained separately.{p_end}
