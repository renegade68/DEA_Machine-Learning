{smcl}
{* Help revision: 2026-10-22.1; command version unchanged.}{...}
{* *! dea_ml  v0.2.2  08sep2026}{...}
{viewerjumpto "Syntax" "dea_ml##syntax"}{...}
{viewerjumpto "Description" "dea_ml##description"}{...}
{viewerjumpto "Methods" "dea_ml##methods"}{...}
{viewerjumpto "Examples" "dea_ml##examples"}{...}
{viewerjumpto "Author" "dea_ml##author"}{...}

{title:Title}

{p2colset 5 18 20 2}{...}
{p2col:{cmd:dea_ml} {hline 2}}Unified machine-learning production frontier estimator{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:dea_ml} {ifin}{cmd:,}
{opt method(name)}
[{opt frontiergraph} {opt fgroptions(string)}]
[{it:method-specific options}]

{pstd}
The required {opt method()} option selects the estimation method.
All remaining options are passed unchanged to the selected sub-command.
Sub-commands are also available as standalone Stata commands.

{pstd}
{opt frontiergraph} draws the fitted production frontier after estimation.
It is available for static specifications with one output and either one or
two inputs.  With one input it plots observed output and the fitted frontier;
with two inputs it draws fitted output as a contour map over the two-input
plane.  The contours use fitted values at the estimation observations;
shading between irregularly located observations uses Stata's
thin-plate-spline interpolation and is a graphical aid, not model extrapolation.

{pstd}
{opt fgroptions(string)} passes quoted Stata graph options to the frontier
graph, for example {cmd:fgroptions("title(EAT) xscale(log) yscale(log)")}.


{marker description}{...}
{title:Description}

{pstd}
{cmd:dea_ml} is a unified interface for the machine-learning production
frontier estimators of the Aparicio research group.
After estimation, all results are accessible via {cmd:e()} exactly as
if the sub-command had been called directly.

{pstd}
Each sub-command is self-contained and can be called independently,
so users who need only one method do not need to load the full package.
New methods can be added to the package without changing the calling
convention.

{pstd}
All commands require Stata 16.1 or later (Mata's {cmd:LinearProgram}
class is used for LP-based methods).


{marker methods}{...}
{title:Methods}

{synoptset 14 tabbed}{...}
{synopthdr:method()}
{synoptline}
{syntab:EAT-based}
{synopt:{helpb eatree}}Efficiency Analysis Trees (EAT){p_end}
{synopt:{helpb ceat}}Convex EAT (CEAT){p_end}
{synopt:{helpb rfeat}}Random Forest EAT (RFEAT){p_end}
{syntab:By-production}
{synopt:{helpb bpeat}}By-production EAT (BP-EAT){p_end}
{synopt:{helpb bpceat}}By-production Convex EAT (BP-CEAT){p_end}
{syntab:Dynamic}
{synopt:{helpb dyneat}}Dynamic EAT{p_end}
{synopt:{helpb dynceat}}Dynamic Convex EAT{p_end}
{syntab:Boosting}
{synopt:{helpb beateat}}Boosted EAT (EATBoosting){p_end}
{synopt:{helpb lsbmafs}}Least Squares Boosting with MARS (LSB-MAFS){p_end}
{syntab:Support Vector Frontier}
{synopt:{helpb svfeat}}SVF on an equal-width dominated-cell grid{p_end}
{synopt:{helpb csvf}}Convexified Support Vector Frontier{p_end}
{synopt:{helpb svfspline}}SVF with splines{p_end}
{synoptline}

{pstd}
For options specific to each method, see the corresponding help file.


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
{phang2}{cmd:. set obs 30}{p_end}
{phang2}{cmd:. set seed 42}{p_end}
{phang2}{cmd:. generate double x1 = 1 + 4*runiform()}{p_end}
{phang2}{cmd:. generate double x2 = 2 + 3*runiform()}{p_end}
{phang2}{cmd:. generate double x = x1}{p_end}
{phang2}{cmd:. generate double y = x1^0.4*x2^0.3*exp(-0.2*runiform())}{p_end}

{pstd}Convex EAT (VRS):{p_end}
{phang2}{cmd:. dea_ml, method(ceat) x(x1 x2) y(y) gen(phi_ceat) rts(vrs)}{p_end}

{pstd}Three-variable frontier display (two inputs and one output):{p_end}
{phang2}{cmd:. dea_ml, method(ceat) x(x1 x2) y(y) gen(phi_graph) frontiergraph fgroptions("xscale(log) yscale(log)")}{p_end}

{pstd}Feature-randomized RFEAT with 500 trees:{p_end}
{phang2}{cmd:. dea_ml, method(rfeat) x(x1 x2) y(y) gen(phi_rfeat) ntrees(500) s(1) seed(42)}{p_end}

{pstd}Boosted EAT:{p_end}
{phang2}{cmd:. dea_ml, method(beateat) x(x) y(y) gen(phi_beateat) iter(100) lr(0.3)}{p_end}

{pstd}LSB-MAFS frontier:{p_end}
{phang2}{cmd:. dea_ml, method(lsbmafs) x(x1 x2) y(y) gen(phi_lsbmafs) maxterms(3) iter(5)}{p_end}

{pstd}SVF with leakage-free tuning:{p_end}
{phang2}{cmd:. dea_ml, method(svfeat) x(x) y(y) gen(phi_svfeat) tune cgrid(1 2) epsgrid(.05) dgrid(1 2) folds(3) seed(2026)}{p_end}

{pstd}SVF with splines:{p_end}
{phang2}{cmd:. dea_ml, method(svfspline) x(x) y(y) gen(phi_svfspline) tune cgrid(1 2) epsgrid(.05) dgrid(1 2) folds(3) seed(2026)}{p_end}
{phang2}{cmd:. predict double frontier_svfspline}{p_end}

{pstd}
After SVF-Spline estimation, {cmd:predict} evaluates the fitted frontier using
the training-sample spline nodes.  See {helpb svfspline} for multiple-output
prediction and numerical diagnostics.  Method-specific options, including
{cmd:solver(highs)} and {cmd:pythonexec()}, are passed through unchanged; the
HiGHS route is restricted to fixed-hyperparameter SVF-Spline fits.{p_end}

{pstd}Sub-commands can also be called directly (equivalent):{p_end}
{phang2}{cmd:. ceat, x(x) y(y) gen(phi_direct_ceat) rts(vrs)}{p_end}
{phang2}{cmd:. beateat, x(x) y(y) gen(phi_direct_beateat) iter(100) lr(0.3)}{p_end}


{phang2}{cmd:. restore}{p_end}
{* DEA_ML_HELP_EXAMPLES_END}{...}


{title:Stored results}

{pstd}
{cmd:dea_ml} returns whatever the selected sub-command stores in {cmd:e()}.
In particular, {cmd:e(cmd)} always contains the sub-command name
(e.g., {cmd:"ceat"}, {cmd:"beateat"}), not {cmd:"dea_ml"}.
Likewise, {cmd:e(cmdline)} records the executed sub-command invocation rather
than the original dispatcher call containing {cmd:method()}.
The package standardizes a common core of stored results, while allowing
method-specific scalars and matrices to be added when needed.

{pstd}
This design keeps wrapper commands and second-stage extensions model-agnostic:
they can rely on the shared baseline metadata and branch only when a
method-specific statistic is genuinely required.


{marker author}{...}
{title:Author}

{pstd}Choonjoo Lee{break}
Korea National Defense University{break}
Nonsan-si, Republic of Korea{break}
bloom.rampike@gmail.com{p_end}

{pstd}For the method references, see the corresponding standalone command help
files.{p_end}

{title:Publication and license}

{pstd}Lee, C. Forthcoming. {it:dea_ml: Machine-learning frontier estimators for efficiency analysis in Stata}.
Accepted for publication in the {browse "https://www.stata-journal.com/":Stata Journal}.{p_end}
{pstd}Author-maintained software under the Academic Free License 3.0 (AFL-3.0).
See LICENSE and LICENSE_CODE.txt in the repository. The acceptance notice refers
to the article; later software revisions are maintained separately.{p_end}
