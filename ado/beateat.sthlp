{smcl}
{* Help revision: 2026-10-22.1; command version unchanged.}{...}
{* *! beateat  v0.1.1  08sep2026}{...}
{viewerjumpto "Syntax" "beateat##syntax"}{...}
{viewerjumpto "Description" "beateat##description"}{...}
{viewerjumpto "Options" "beateat##options"}{...}
{viewerjumpto "Stored results" "beateat##results"}{...}
{viewerjumpto "Examples" "beateat##examples"}{...}
{viewerjumpto "References" "beateat##references"}{...}

{title:Title}
{p2colset 5 18 20 2}{...}
{p2col:{cmd:beateat} {hline 2}}Boosted EAT (EATBoosting){p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}
{p 8 16 2}
{cmd:beateat} {ifin}{cmd:,} {opt x(varlist)} {opt y(varlist)}
[{opt gen(name)} {opt numleaves(#)} {opt iter(#)} {opt lr(#)} {opt numstop(#)} {opt seed(#)}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt x(varlist)}}nonnegative input variables{p_end}
{synopt:{opt y(varlist)}}nonnegative outputs; at least one positive per observation{p_end}
{syntab:Optional}
{synopt:{opt gen(name)}}new variable to store efficiency score phi (>= 1){p_end}
{synopt:{opt numleaves(#)}}leaves per EAT weak learner; default {cmd:numleaves(5)}{p_end}
{synopt:{opt iter(#)}}boosting iterations; default {cmd:iter(100)}{p_end}
{synopt:{opt lr(#)}}learning rate (shrinkage) in (0,1]; default {cmd:lr(0.3)}{p_end}
{synopt:{opt numstop(#)}}minimum observations in a node before a split is attempted; default {cmd:numstop(5)}{p_end}
{synopt:{opt seed(#)}}accepted for backward compatibility; the estimator is deterministic and preserves the RNG state{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:beateat} estimates a production frontier using gradient boosting with EAT
(Efficiency Analysis Tree) weak learners and the squared-error (L2) loss. The
frontier is progressively refined over {opt iter()} boosting rounds, each adding a
shrunk EAT fit to the residuals.{p_end}

{pstd}
The boosting algorithm proceeds as follows:{p_end}
{phang2}1. Initialize F_0 = max(y) for each output (conservative upper bound).{p_end}
{phang2}2. For m = 1, ..., {it:iter}:{p_end}
{phang2}{space 4}a. Compute residuals r = Y - F (r <= 0 initially since F >= Y).{p_end}
{phang2}{space 4}b. Shift residuals: r_min = min(r); r_pos = r - r_min (non-negative).{p_end}
{phang2}{space 4}c. Fit EAT weak learner h to r_pos (produces non-negative predictions).{p_end}
{phang2}{space 4}d. Update F = F + lr * (h + r_min) (reduces frontier toward data).{p_end}
{phang2}3. Compute phi_i = max(1, min_r{F_r(x_i) / y_ir}) for each DMU.{p_end}

{pstd}The residual shift is a computational translation for the nonnegative-response
EAT routine. Adding r_min back in the update makes it equivalent to fitting the
unshifted residuals in the published EATBoost algorithm. The lower bound of one
in the reported score is a numerical safeguard and is inactive under exact
envelopment.{p_end}

{pstd}
The efficiency score phi >= 1, where phi = 1 indicates full efficiency. Each EAT
weak learner is constrained to {opt numleaves()} terminal nodes and fit using the
{opt numstop()} stopping criterion.{p_end}

{marker options}{...}
{title:Options}

{phang}{opt x(varlist)} specifies the input variables. Required.{p_end}

{phang}{opt y(varlist)} specifies the output variables. Required.{p_end}

{phang}{opt gen(name)} stores the efficiency score phi (>= 1) in the new variable
{it:name}. If omitted, results are stored in {cmd:e()} only.{p_end}

{phang}{opt numleaves(#)} sets the maximum number of terminal leaf nodes for each
EAT weak learner. Default is 5.{p_end}

{phang}{opt iter(#)} sets the number of boosting iterations (weak learners to add).
Default is 100.{p_end}

{phang}{opt lr(#)} sets the learning rate (shrinkage factor) applied to each weak
learner's contribution. Must be in (0,1]. Smaller values require more iterations but
may improve generalization. Default is 0.3.{p_end}

{phang}{opt numstop(#)} sets the minimum number of observations in a node before
a weak-learner split is attempted. Default is 5.{p_end}

{phang}{opt seed(#)} is accepted for backward compatibility. The published
implementation is deterministic, and this option preserves Stata's RNG state.{p_end}

{marker results}{...}
{title:Stored results}

{pstd}{cmd:beateat} stores the following in {cmd:e()}:{p_end}
{synoptset 18 tabbed}{...}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(sample)}}estimation-sample indicator{p_end}
{synopt:{cmd:e(mean_phi)}}mean efficiency score across all DMUs{p_end}
{synopt:{cmd:e(n_iter)}}number of boosting iterations used{p_end}
{synopt:{cmd:e(n_leaves)}}number of leaves per weak learner{p_end}
{synopt:{cmd:e(lr)}}learning rate used{p_end}
{synopt:{cmd:e(frontier)}}fitted output frontier at the estimation observations{p_end}
{synopt:{cmd:e(x)}}names of input variables{p_end}
{synopt:{cmd:e(y)}}names of output variables{p_end}
{synopt:{cmd:e(cmd)}}"beateat"{p_end}
{synopt:{cmd:e(cmdline)}}executed sub-command invocation{p_end}

{marker examples}{...}
{title:Examples}

{pstd}Run the complete example block in the Do-file Editor, omitting leading dots.
The examples create their own data; {cmd:preserve} and {cmd:restore} retain your
original dataset on successful completion. Settings and estimation results may change.
Small samples and tuning grids demonstrate syntax rather than a research design.{p_end}
{* DEA_ML_HELP_EXAMPLES_BEGIN}{...}
{phang2}{cmd:. preserve}{p_end}

{pstd}Create reproducible data supporting all examples below:{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 80}{p_end}
{phang2}{cmd:. set seed 42}{p_end}
{phang2}{cmd:. generate double x1 = 1 + 4*runiform()}{p_end}
{phang2}{cmd:. generate double x2 = 2 + 3*runiform()}{p_end}
{phang2}{cmd:. generate double x = x1}{p_end}
{phang2}{cmd:. generate double y = x1^0.4*x2^0.3*exp(-0.2*runiform())}{p_end}
{phang2}{cmd:. generate double y1 = y}{p_end}
{phang2}{cmd:. generate double y2 = x^0.35*exp(-0.2*runiform())}{p_end}

{pstd}Basic boosted EAT with default settings:{p_end}
{phang2}{cmd:. beateat, x(x) y(y) gen(phi_basic)}{p_end}

{pstd}Increase iterations and reduce the learning rate for a less irregular fitted frontier:{p_end}
{phang2}{cmd:. beateat, x(x1 x2) y(y) gen(phi_smooth) iter(200) lr(0.1) numleaves(3)}{p_end}

{pstd}Multiple outputs:{p_end}
{phang2}{cmd:. beateat, x(x) y(y1 y2) gen(phi_multi) iter(100) lr(0.3)}{p_end}

{pstd}Display mean efficiency after estimation:{p_end}
{phang2}{cmd:. beateat, x(x) y(y) gen(phi_display)}{p_end}
{phang2}{cmd:. display e(mean_phi)}{p_end}

{phang2}{cmd:. restore}{p_end}
{* DEA_ML_HELP_EXAMPLES_END}{...}


{marker references}{...}
{title:References}

{phang}Guillen, M. D., J. Aparicio, and M. Esteve. 2023. Gradient tree boosting
and the estimation of production frontiers.
{it:Expert Systems with Applications} 214: 119134.
{browse "https://doi.org/10.1016/j.eswa.2022.119134"}
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
