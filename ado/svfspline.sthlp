{smcl}
{* Help revision: 2026-10-22.1; command version unchanged.}{...}
{* *! svfspline  v0.3.0  13sep2026}{...}
{viewerjumpto "Syntax" "svfspline##syntax"}{...}
{viewerjumpto "Description" "svfspline##description"}{...}
{viewerjumpto "Options" "svfspline##options"}{...}
{viewerjumpto "Stored results" "svfspline##results"}{...}
{viewerjumpto "Postestimation prediction" "svfspline##predict"}{...}
{viewerjumpto "Examples" "svfspline##examples"}{...}
{viewerjumpto "Specification and reproducibility" "svfspline##reproducibility"}{...}
{viewerjumpto "References" "svfspline##references"}{...}

{title:Title}
{p2colset 5 18 20 2}{...}
{p2col:{cmd:svfspline} {hline 2}}Support Vector Frontier with Splines (SVF-Spline){p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}
{p 8 16 2}
{cmd:svfspline} [{it:if}] [{it:in}]{cmd:,} {opt x(varlist)} {opt y(varlist)}
[{opt gen(name)} {opt nknots(#)} {opt c(#)} {opt eps(#)} {opt tune}
 {opt cgrid(numlist)} {opt epsgrid(numlist)} {opt dgrid(numlist)}
 {opt folds(#)} {opt seed(#)} {opt lambda(#)} {opt solver(native|highs)}
 {opt pythonexec(string)}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt x(varlist)}}nonnegative input variables{p_end}
{synopt:{opt y(varlist)}}strictly positive output variables{p_end}
{syntab:Optional}
{synopt:{opt gen(name)}}new variable to store efficiency score phi (>= 1){p_end}
{synopt:{opt nknots(#)}}equal-width intervals d per input; creates d+1 nodes; default {cmd:nknots(3)}{p_end}
{synopt:{opt c(#)}}slack penalty; default {cmd:c(10.0)}{p_end}
{synopt:{opt eps(#)}}epsilon-insensitivity band; default {cmd:eps(0.0)}{p_end}
{synopt:{opt lambda(#)}}deprecated compatibility option; ignored{p_end}
{synopt:{opt tune}}select C, epsilon, and d by leakage-free cross-validation{p_end}
{synopt:{opt cgrid(numlist)}}candidate C values for {cmd:tune}{p_end}
{synopt:{opt epsgrid(numlist)}}candidate epsilon values for {cmd:tune}{p_end}
{synopt:{opt dgrid(numlist)}}candidate interval counts d for {cmd:tune}{p_end}
{synopt:{opt folds(#)}}number of CV folds; default 5{p_end}
{synopt:{opt seed(#)}}CV fold seed; default {cmd:seed(12345)}{p_end}
{synopt:{opt solver(native|highs)}}numerical backend; default {cmd:native}{p_end}
{synopt:{opt pythonexec(string)}}Python executable for optional {cmd:highs} backend{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:svfspline} estimates the continuous, piecewise-linear, concave,
non-decreasing SVF-SP production frontier in Guerrero et al. (2024). It uses an
additive truncated-linear spline transformation and squared-L2 regularization:
the single-output formulation is their model (14), and the multiple-output
formulation is their separable model (19), with a common input grid.
The frontier is specified as:{p_end}
{phang2}f(x) = sum_j [beta0_j + beta_j*x_j + sum_k beta_jk*max(0,x_j-t_jk)]{p_end}
{pstd}
where cumulative segment slopes are nonnegative and hinge coefficients beta_jk
are nonpositive. These restrictions enforce monotonicity and concavity. Each
input range is divided into {opt nknots()} equal-width intervals, producing
{opt nknots()}+1 endpoint nodes.{p_end}

{pstd}
The convex program minimizes:{p_end}
{phang2}beta'beta + C * sum(xi){p_end}
{pstd}
subject to y_i <= f(x_i) <= y_i + eps + xi_i, xi_i >= 0, and the shape
constraints above. The command uses Mata {cmd:optimize()} with a decreasing
logarithmic barrier and rejects solutions whose original constraint residuals
exceed 1e-6.  It separately reports approximate KKT stationarity and a barrier
duality-gap bound.  Values above 1e-2 and 1e-3, respectively, trigger a visible
diagnostic warning and should be checked with an independent QP solver.{p_end}

{pstd}
The efficiency score is phi = f(x_i) / y_i >= 1, where phi = 1 indicates full
efficiency on the estimated frontier.{p_end}

{pstd}
Requires Stata 16 or later. Stata 16--18 use the documented
{cmd:optimize()} log-barrier implementation.{p_end}

{marker options}{...}
{title:Options}

{phang}{opt x(varlist)} specifies nonnegative input variables. Required.{p_end}

{phang}{opt y(varlist)} specifies strictly positive output variables. Required.{p_end}

{phang}{opt gen(name)} stores the efficiency score phi (>= 1) in the new variable
{it:name}. If omitted, results are stored in {cmd:e()} only.{p_end}

{phang}{opt nknots(#)} sets d, the number of equal-width intervals in each input
range. There are d+1 nodes and p(d+3) coefficients for p inputs. Default is 3.{p_end}

{phang}{opt c(#)} sets the slack penalty controlling the tightness of the frontier
fit. Larger values penalize violations more heavily. Default is 10.0.{p_end}

{phang}{opt eps(#)} sets the width of the epsilon-insensitivity band. Observations
within {it:#} units of the frontier incur no penalty. Default is 0.0.{p_end}

{phang}{opt lambda(#)} is retained temporarily for command-line compatibility
but is ignored. Neither model (14) nor its multiple-output extension (19) has a
separate lambda hyperparameter.{p_end}

{phang}{opt solver(native|highs)} selects the numerical backend. {cmd:native}
is the default Mata log-barrier implementation. {cmd:highs} is an optional,
fixed-hyperparameter sparse-QP backend and cannot be combined with {cmd:tune}.

{phang}{opt pythonexec(string)} gives the full path to Python with {cmd:cvxpy}
and {cmd:highspy} installed. It is required with {cmd:solver(highs)}. The
command creates no scores unless HiGHS reports {cmd:optimal} and independently
recomputed primal and KKT residuals are no larger than 1e-7.{p_end}

{phang}{opt tune} selects C, epsilon, and d using held-out mean squared error.
For every fold, nodes and coefficients are estimated using the training subset
only. The selected model is then refitted on all estimation observations.{p_end}

{phang}{opt cgrid(numlist)}, {opt epsgrid(numlist)}, and {opt dgrid(numlist)}
specify CV candidates. Defaults are respectively
{cmd:.1 .5 1 2 5}, {cmd:0 .001 .01 .1 .2}, and {cmd:1 2 3}.
The command minimizes mean held-out MSE and retains the first candidate in
d--C--epsilon scan order when losses differ by no more than 1e-14.{p_end}

{phang}{opt folds(#)} sets the number of CV folds. Default is 5. It may not
exceed the number of complete estimation observations.{p_end}

{phang}{opt seed(#)} sets the random fold assignment under {cmd:tune}. It does
not affect a fixed-hyperparameter fit. Default is 12345.{p_end}

{marker results}{...}
{title:Stored results}

{pstd}{cmd:svfspline} stores the following in {cmd:e()}:{p_end}
{synoptset 18 tabbed}{...}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(sample)}}estimation-sample indicator{p_end}
{synopt:{cmd:e(mean_phi)}}mean efficiency score across all DMUs{p_end}
{synopt:{cmd:e(n_failed)}}number of output-specific optimization failures{p_end}
{synopt:{cmd:e(n_basis)}}total number of spline basis functions used{p_end}
{synopt:{cmd:e(nknots)}}equal-width interval count d used{p_end}
{synopt:{cmd:e(c_pen)}}slack penalty used{p_end}
{synopt:{cmd:e(eps)}}epsilon-insensitivity band used{p_end}
{synopt:{cmd:e(objective)}}sum of output-specific objective values{p_end}
{synopt:{cmd:e(max_violation)}}maximum original inequality violation{p_end}
{synopt:{cmd:e(max_stationarity)}}maximum normalized approximate KKT stationarity residual{p_end}
{synopt:{cmd:e(max_barrier_gap)}}maximum relative barrier duality-gap bound{p_end}
{synopt:{cmd:e(all_stages_converged)}}1 if every intermediate barrier stage reported convergence{p_end}
{synopt:{cmd:e(barrier_stage_warning)}}1 if a native intermediate barrier stage reached its iteration cap; an informational review flag, 0 for HiGHS{p_end}
{synopt:{cmd:e(numerical_warning)}}1 if the final primal, KKT-stationarity, or barrier-gap diagnostic exceeds its recommended threshold; inspect results before use{p_end}
{synopt:{cmd:e(max_complementarity)}}maximum normalized KKT complementarity residual; HiGHS backend{p_end}
{synopt:{cmd:e(solver_iterations)}}HiGHS QP iterations; HiGHS backend{p_end}
{synopt:{cmd:e(solver)}}solver identifier{p_end}
{synopt:{cmd:e(coefficients)}}spline coefficient matrix{p_end}
{synopt:{cmd:e(frontier)}}fitted frontier matrix{p_end}
{synopt:{cmd:e(xranges)}}training-sample input minima and maxima defining spline nodes{p_end}
{synopt:{cmd:e(n_inputs)}}number of input variables{p_end}
{synopt:{cmd:e(n_outputs)}}number of output variables{p_end}
{synopt:{cmd:e(tuned)}}1 if cross-validation was requested{p_end}
{synopt:{cmd:e(folds)}}number of CV folds, or 0 for a fixed fit{p_end}
{synopt:{cmd:e(cv_failures)}}number of candidates with at least one failed fold{p_end}
{synopt:{cmd:e(cv_results)}}candidate C, epsilon, d, and mean validation MSE; returned with {cmd:tune}{p_end}
{synopt:{cmd:e(cv_fold_results)}}candidate-by-fold validation MSE; returned with {cmd:tune}{p_end}
{synopt:{cmd:e(cv_predictions)}}held-out predictions with candidate/fold/observation identifiers; returned with {cmd:tune}{p_end}
{synopt:{cmd:e(fold_id)}}fold assignment by estimation-sample row; returned with {cmd:tune}{p_end}
{synopt:{cmd:e(x)}}names of input variables{p_end}
{synopt:{cmd:e(y)}}names of output variables{p_end}
{synopt:{cmd:e(cmd)}}"svfspline"{p_end}
{synopt:{cmd:e(cmdline)}}executed sub-command invocation{p_end}

{marker predict}{...}
{title:Postestimation prediction}

{pstd}
After estimation, {cmd:predict newvar} evaluates the fitted frontier for the
first output.  With multiple outputs, use {cmd:predict newvar, output(#)}.
Prediction reuses the nodes defined by {cmd:e(xranges)}; it does not estimate
new nodes from the prediction observations.  This permits a frontier fitted on
a training sample to be evaluated on a held-out or fixed test panel.{p_end}

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
{phang2}{cmd:. set obs 30}{p_end}
{phang2}{cmd:. set seed 42}{p_end}
{phang2}{cmd:. generate double x1 = 1 + 4*runiform()}{p_end}
{phang2}{cmd:. generate double x2 = 2 + 3*runiform()}{p_end}
{phang2}{cmd:. generate double x = x1}{p_end}
{phang2}{cmd:. generate double y = x1^0.4*x2^0.3*exp(-0.2*runiform())}{p_end}

{pstd}Basic SVF-Spline with default settings:{p_end}
{phang2}{cmd:. svfspline, x(x) y(y) gen(phi_basic)}{p_end}
{pstd}After a native fit, inspect both {cmd:e(numerical_warning)} and
{cmd:e(barrier_stage_warning)}. A score is produced only when the original
shape and fit constraints are feasible. The latter flag records an intermediate
iteration cap even if the final diagnostics are acceptable; the former is the
review warning for an unacceptable final primal, KKT, or barrier-gap diagnostic.
For a fixed specification requiring an independent QP check, use the optional
{cmd:solver(highs)} route described above.{p_end}

{pstd}With {cmd:tune} and no grids, five-fold cross-validation fits 75
candidate specifications (375 training fits). This can take tens of minutes
or longer. The runnable example below uses a smaller instructional grid.{p_end}

{pstd}Use a small instructional candidate grid:{p_end}
{phang2}{cmd:. svfspline, x(x) y(y) gen(phi_small_tune) tune cgrid(1 2) epsgrid(.05) dgrid(1 2) folds(3) seed(2026)}{p_end}

{pstd}Inspect CV results:{p_end}
{phang2}{cmd:. matrix list e(cv_results)}{p_end}
{phang2}{cmd:. matrix list e(cv_fold_results)}{p_end}

{pstd}Display basis size and mean efficiency:{p_end}
{phang2}{cmd:. svfspline, x(x) y(y) gen(phi_display)}{p_end}
{phang2}{cmd:. display "Basis functions: " e(n_basis)}{p_end}
{phang2}{cmd:. display "Mean phi: " e(mean_phi)}{p_end}

{pstd}Fit on a training subsample and predict for all observations:{p_end}
{phang2}{cmd:. generate byte train = mod(_n,5) != 0}{p_end}
{phang2}{cmd:. svfspline if train, x(x1 x2) y(y) nknots(2) c(1) eps(.05)}{p_end}
{phang2}{cmd:. predict double frontier_oos}{p_end}

{phang2}{cmd:. restore}{p_end}
{* DEA_ML_HELP_EXAMPLES_END}{...}


{marker reproducibility}{...}
{title:Specification and reproducibility note}

{pstd}
The squared-L2 objective implemented here is also used by the authors' public
Python implementation.  The displayed norm in equations (14) and (19) of
Guerrero et al. (2024) is typographically ambiguous because its exponent is
not printed. This command follows the executable squared-L2 specification.
Their {it:q}=1 denotes the order of the spline basis, not an L1 regularizer.
The directional-distance LP in their equation (26) is a postestimation
efficiency calculation on the fitted virtual-point technology; it is not an
LP reformulation of the SVF-Spline fitting problem.{p_end}

{pstd}The optional HiGHS route solves the same fixed-hyperparameter convex QP.
It is opt-in while its multi-output and cross-validation performance contracts
are expanded and independently tested. For this route,
{cmd:e(max_barrier_gap)} is missing and {cmd:e(max_complementarity)} is returned.{p_end}

{pstd}
This implementation deliberately recomputes nodes and coefficients on
the training rows of each CV fold and retains every input block in prediction.
Those safeguards are stated explicitly because the pinned public Python source
reviewed at the pinned commit contains CV data leakage and a first-input-only
prediction path for multiple-input data.  These defects are not reproduced.{p_end}

{pstd}
The reported efficiency is the conventional frontier ratio f(x_i)/y_i.  The
epsilon option affects estimation of f; epsilon is not subtracted from the
reported score.{p_end}

{marker references}{...}
{title:References}

{phang}Guerrero, N. M., Moragues, J., Aparicio, J., and Valero-Carreras, D.
(2024). Support Vector Frontiers with kernel splines.
{it:Omega} 128:103130.
{browse "https://doi.org/10.1016/j.omega.2024.103130"}
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
