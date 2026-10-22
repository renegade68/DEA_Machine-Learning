{smcl}
{* Help revision: 2026-10-22.1; command version unchanged.}{...}
{* *! lsbmafs  v0.1.1  08sep2026}{...}
{viewerjumpto "Syntax" "lsbmafs##syntax"}{...}
{viewerjumpto "Description" "lsbmafs##description"}{...}
{viewerjumpto "Options" "lsbmafs##options"}{...}
{viewerjumpto "Stored results" "lsbmafs##results"}{...}
{viewerjumpto "Examples" "lsbmafs##examples"}{...}
{viewerjumpto "References" "lsbmafs##references"}{...}

{title:Title}
{p2colset 5 18 20 2}{...}
{p2col:{cmd:lsbmafs} {hline 2}}Least Squares Boosting with MARS (LSB-MAFS){p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}
{p 8 16 2}
{cmd:lsbmafs} {ifin}{cmd:,} {opt x(varlist)} {opt y(varlist)}
[{opt gen(name)} {opt maxterms(#)} {opt iter(#)} {opt lr(#)} {opt numstop(#)} {opt seed(#)}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt x(varlist)}}nonnegative input variables{p_end}
{synopt:{opt y(varlist)}}one strictly positive output variable{p_end}
{syntab:Optional}
{synopt:{opt gen(name)}}new variable to store efficiency score phi (>= 1){p_end}
{synopt:{opt maxterms(#)}}maximum MARS basis terms per weak learner; default {cmd:maxterms(5)}{p_end}
{synopt:{opt iter(#)}}boosting iterations; default {cmd:iter(100)}{p_end}
{synopt:{opt lr(#)}}learning rate (shrinkage) in (0,1]; default {cmd:lr(0.1)}{p_end}
{synopt:{opt numstop(#)}}optional MAFS minimum-span override; by default computed from sample size{p_end}
{synopt:{opt seed(#)}}accepted for backward compatibility; the estimator is deterministic and preserves the RNG state{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:lsbmafs} estimates a production frontier using gradient boosting (L2 loss) with
Multivariate Adaptive Regression Splines (MARS) as weak learners. Following the
published adapted forward MARS (MAFS) algorithm, the command implements the rough
piecewise-linear LSB-MAFS frontier without the paper's optional post-fit smoothing
transformation. Each selected knot contributes the
reflected pair max(0, x_j-t) and max(0, t-x_j). Their coefficients are estimated by
linear programming under envelopment, monotonicity, and concavity restrictions.{p_end}

{pstd}
The boosting algorithm mirrors {helpb beateat} but replaces the EAT weak learner with
MARS:{p_end}
{phang2}1. Initialize F_0 = max(y).{p_end}
{phang2}2. For m = 1, ..., {it:iter}:{p_end}
{phang2}{space 4}a. Compute residuals r = Y - F (r <= 0 initially since F >= Y).{p_end}
{phang2}{space 4}b. Fit a constrained MAFS weak learner directly to r, adding
reflected hinge pairs by forward selection.{p_end}
{phang2}{space 4}c. Update F = F + lr * h.{p_end}
{phang2}3. Compute phi_i = F(x_i)/y_i for each DMU.{p_end}

{pstd}
MARS hinge functions provide a continuous, piecewise-linear approximation to the
frontier, in contrast to the step-function approximation of EAT. Because MARS
captures non-linear monotone relationships through additive spline terms, it may
offer better approximation with fewer weak learners than {helpb beateat} in some
applications.{p_end}

{pstd}
The theoretical score is phi = F(x_i)/y_i. The command uses a lower bound of one
only as a numerical safeguard; phi = 1 indicates full efficiency.{p_end}

{marker options}{...}
{title:Options}

{phang}{opt x(varlist)} specifies the input variables. Required.{p_end}

{phang}{opt y(varlist)} specifies exactly one output variable, as required by the
published LSB-MAFS estimator and its reference implementation. Required.{p_end}

{phang}{opt gen(name)} stores the efficiency score phi (>= 1) in the new variable
{it:name}. If omitted, results are stored in {cmd:e()} only.{p_end}

{phang}{opt maxterms(#)} sets the maximum number of MARS basis terms (hinge
functions) per weak learner. Larger values allow more flexible but potentially
overfit weak learners. Default is 5.{p_end}

{phang}{opt iter(#)} sets the number of boosting iterations (MARS weak learners to
add). Default is 100.{p_end}

{phang}{opt lr(#)} sets the learning rate (shrinkage factor) applied to each MARS
weak learner's contribution. Must be in (0,1]. Smaller values require more
iterations but may improve generalization. Default is 0.1.{p_end}

{phang}{opt numstop(#)} overrides the minimum span used to separate candidate knots
in the adapted forward MARS search. If omitted, the span is computed from the sample
size using the rule in the reference implementation. Usually this option should be
left at its automatic setting.{p_end}

{phang}{opt seed(#)} is accepted for backward compatibility. The published
implementation is deterministic, and this option preserves Stata's RNG state.{p_end}

{pstd}Failed candidate LPs are skipped and counted. A failed selected-basis LP
halts estimation before any requested output variable is created.{p_end}

{marker results}{...}
{title:Stored results}

{pstd}{cmd:lsbmafs} stores the following in {cmd:e()}:{p_end}
{synoptset 18 tabbed}{...}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(sample)}}estimation-sample indicator{p_end}
{synopt:{cmd:e(n_lp_failed)}}number of failed candidate or selected-basis LP solves{p_end}
{synopt:{cmd:e(n_final_failed)}}number of failed selected-basis LP solves{p_end}
{synopt:{cmd:e(solver)}}LP solver identifier{p_end}
{synopt:{cmd:e(mean_phi)}}mean efficiency score across all DMUs{p_end}
{synopt:{cmd:e(n_iter)}}number of boosting iterations used{p_end}
{synopt:{cmd:e(maxterms)}}maximum MARS basis terms per weak learner used{p_end}
{synopt:{cmd:e(lr)}}learning rate used{p_end}
{synopt:{cmd:e(frontier)}}fitted output frontier at the estimation observations{p_end}
{synopt:{cmd:e(x)}}names of input variables{p_end}
{synopt:{cmd:e(y)}}names of output variables{p_end}
{synopt:{cmd:e(cmd)}}"lsbmafs"{p_end}
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
{phang2}{cmd:. set obs 30}{p_end}
{phang2}{cmd:. set seed 42}{p_end}
{phang2}{cmd:. generate double x1 = 1 + 4*runiform()}{p_end}
{phang2}{cmd:. generate double x2 = 2 + 3*runiform()}{p_end}
{phang2}{cmd:. generate double x = x1}{p_end}
{phang2}{cmd:. generate double y = x1^0.4*x2^0.3*exp(-0.2*runiform())}{p_end}

{pstd}{it:Note:} The published forward-search specification solves many constrained linear
programs. Large samples, more basis terms, or many boosting iterations can take
substantially longer than the small examples below.{p_end}

{pstd}Basic LSB-MAFS with a short instructional run:{p_end}
{phang2}{cmd:. lsbmafs, x(x) y(y) gen(phi_basic) maxterms(3) iter(5)}{p_end}

{pstd}More MARS basis terms and more iterations with a slower learning rate:{p_end}
{phang2}{cmd:. lsbmafs, x(x1 x2) y(y) gen(phi_terms) maxterms(5) iter(10) lr(0.05)}{p_end}

{pstd}Compare LSB-MAFS and boosted EAT on the same data:{p_end}
{phang2}{cmd:. beateat, x(x) y(y) gen(phi_eat) iter(5) lr(0.1)}{p_end}
{phang2}{cmd:. lsbmafs, x(x) y(y) gen(phi_mars) maxterms(3) iter(5) lr(0.1)}{p_end}
{phang2}{cmd:. summarize phi_eat phi_mars}{p_end}

{pstd}Display mean efficiency after estimation:{p_end}
{phang2}{cmd:. lsbmafs, x(x) y(y) gen(phi_display) maxterms(3) iter(5)}{p_end}
{phang2}{cmd:. display e(mean_phi)}{p_end}

{phang2}{cmd:. restore}{p_end}
{* DEA_ML_HELP_EXAMPLES_END}{...}


{marker references}{...}
{title:References}

{phang}Guillen, M. D., J. Aparicio, J. L. Zofio, and V. J. Espana. 2024.
Improving the predictive accuracy of production frontier models for efficiency
measurement using machine learning: The LSB-MAFS method.
{it:Computers & Operations Research} 171: 106793.
{browse "https://doi.org/10.1016/j.cor.2024.106793"}
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
